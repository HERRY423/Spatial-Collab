"""OAuth resource server with explicit subject/project ACLs; no anonymous public mode.

The deployment supplies an OAuth authorization server. This module validates its
RS256 access tokens and advertises discovery; it never invents an identity provider.
"""
from __future__ import annotations

import argparse
from contextlib import AsyncExitStack, asynccontextmanager
import json
from pathlib import Path
import re
import threading
import time
from urllib.parse import urlsplit

from mcp.server.auth.middleware.auth_context import get_access_token
from mcp.server.auth.provider import AccessToken
from mcp.server.auth.settings import AuthSettings
from mcp.server.transport_security import TransportSecuritySettings
from starlette.applications import Starlette
from starlette.concurrency import run_in_threadpool
from starlette.middleware.trustedhost import TrustedHostMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse
from starlette.routing import Mount, Route

from .server import READ_TOOLS, ToolService, create_server
from .store import SpatialError


def https_url(value):
    p = urlsplit(value)
    if p.scheme != "https" or not p.hostname or p.username or p.password or p.query or p.fragment:
        raise ValueError("OAuth issuer, JWKS and public URLs must be explicit HTTPS URLs without credentials, query or fragment.")
    return value


class JWTVerifier:
    def __init__(self, issuer, audience, jwks_url, subjects, *, key_provider=None):
        import jwt
        self.issuer = https_url(issuer)
        self.audience = https_url(audience)
        self.subjects = subjects
        self.keys = key_provider or jwt.PyJWKClient(https_url(jwks_url), lifespan=300, timeout=5)

    def decode(self, token):
        import jwt
        key = self.keys.get_signing_key_from_jwt(token).key
        claims = jwt.decode(token, key, algorithms=["RS256"], audience=self.audience, issuer=self.issuer,
                            options={"require": ["exp", "iat", "sub", "iss", "aud"]})
        if not isinstance(claims["sub"], str) or claims["sub"] not in self.subjects:
            return None
        if claims["exp"] - claims["iat"] > 3600:
            return None
        scope = claims.get("scope", "")
        if not isinstance(scope, str):
            return None
        scopes = scope.split()
        if "spatial:read" not in scopes:
            return None
        return AccessToken(token=token, client_id=str(claims.get("client_id", claims.get("azp", ""))),
                           subject=claims["sub"], scopes=scopes, expires_at=int(claims["exp"]),
                           resource=self.audience, claims={"iss": self.issuer})

    async def verify_token(self, token):
        if len(token) > 16384:
            return None
        try:
            return await run_in_threadpool(self.decode, token)
        except Exception:
            # Never include token, identity-provider response or claims in a log.
            return None


class AuthorizedService(ToolService):
    def __init__(self, project_root, subjects, *, max_jobs=2):
        super().__init__(project_root)
        self.subjects = subjects
        self.max_jobs = max_jobs
        self.audit_lock = threading.Lock()

    def audit(self, token, name, outcome):
        record = {"timestamp": time.time(), "subject": token.subject if token else None,
                  "tool": name, "outcome": outcome}
        # Attribution is taken from the verified token, separately from researcher
        # names. No arguments, raw expression, credentials or file bytes are logged.
        with self.audit_lock, (self.project.root / "access-audit.jsonl").open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record, ensure_ascii=False)+"\n")

    def call(self, name, arguments=None):
        token = get_access_token()
        role = self.subjects.get(token.subject) if token else None
        if role not in {"read", "write"} or "spatial:read" not in token.scopes:
            self.audit(token, name, "denied")
            raise SpatialError("This identity is not authorized for this project.")
        if name not in READ_TOOLS and (role != "write" or "spatial:write" not in token.scopes):
            self.audit(token, name, "denied_write")
            raise SpatialError("This operation requires project write access and spatial:write scope.")
        if name in {"submit_analysis", "retry_analysis_job", "submit_integration", "retry_integration_job"}:
            from .transfers import projects
            active = 0
            for project in projects(self.project):
                with project._db() as db:
                    for table in ("analysis_jobs", "integration_jobs"):
                        if db.execute("SELECT 1 FROM sqlite_master WHERE name=?", (table,)).fetchone():
                            active += db.execute(f"SELECT count(*) FROM {table} WHERE status IN ('queued','running','cancelling')").fetchone()[0]
            if active >= self.max_jobs:
                raise SpatialError("Project job quota reached; wait for a job to finish or cancel it.")
        self.audit(token, name, "started")
        try:
            result = super().call(name, arguments)
        except Exception:
            self.audit(token, name, "failed")
            raise
        self.audit(token, name, "succeeded")
        return result


class RequestLimits:
    def __init__(self, app, max_bytes=2_000_000, concurrent=8, body_timeout=15):
        import asyncio
        self.app, self.max_bytes = app, max_bytes
        self.semaphore = asyncio.Semaphore(concurrent)
        self.body_timeout = body_timeout

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            return await self.app(scope, receive, send)
        if self.semaphore.locked():
            return await JSONResponse({"error": "Server busy; retry later."}, 429)(scope, receive, send)
        async with self.semaphore:
            # Buffer only bounded request data, before MCP parsing. No unlimited
            # body read for chunked requests or forged Content-Length headers.
            messages, total = [], 0
            import asyncio
            deadline = asyncio.get_running_loop().time() + self.body_timeout
            while True:
                try:
                    message = await asyncio.wait_for(receive(), max(0, deadline - asyncio.get_running_loop().time()))
                except TimeoutError:
                    return await JSONResponse({"error": "Request body timed out."}, 408)(scope, receive, send)
                if message["type"] == "http.disconnect":
                    return
                total += len(message.get("body", b""))
                if total > self.max_bytes:
                    return await JSONResponse({"error": "Request too large."}, 413)(scope, receive, send)
                messages.append(message)
                if not message.get("more_body", False):
                    break
            async def replay():
                return messages.pop(0) if messages else await receive()
            await self.app(scope, replay, send)


def create_gateway(settings, *, key_provider=None):
    origin = https_url(settings["public_origin"]).rstrip("/")
    if urlsplit(origin).path:
        raise ValueError("public_origin must contain only scheme and host.")
    issuer = https_url(settings["issuer"])
    jwks = https_url(settings["jwks_url"])
    routes, servers, roots, specs = [], [], [], {}
    for alias, entry in settings["projects"].items():
        if not re.fullmatch(r"[a-z0-9-]{1,64}", alias):
            raise ValueError("Project aliases must be bounded lower-case names.")
        root = Path(entry["path"]).resolve(strict=True)
        if any(root.is_relative_to(other) or other.is_relative_to(root) for other in roots):
            raise ValueError("Configured project roots must not overlap.")
        roots.append(root)
        subjects = entry["subjects"]
        if not subjects or any(not isinstance(s, str) or not s or r not in {"read", "write"} for s, r in subjects.items()):
            raise ValueError("Each project needs explicit OAuth subject-to-role mappings.")
        resource = origin + f"/projects/{alias}/mcp"
        max_jobs = entry.get("max_jobs", 2)
        if type(max_jobs) is not int or not 1 <= max_jobs <= 100:
            raise ValueError("max_jobs must be an integer from 1 to 100.")
        service = AuthorizedService(root, subjects, max_jobs=max_jobs)
        verifier = JWTVerifier(issuer, resource, jwks, subjects, key_provider=key_provider)
        server = create_server(root, service=service, token_verifier=verifier,
                               auth=AuthSettings(issuer_url=issuer, resource_server_url=resource, required_scopes=["spatial:read"]),
                               transport_security=TransportSecuritySettings(enable_dns_rebinding_protection=True,
                                   allowed_hosts=[urlsplit(origin).netloc], allowed_origins=[origin]))
        routes.append(Mount(f"/projects/{alias}", app=server.streamable_http_app()))
        servers.append(server)
        specs[alias] = {"resource": resource, "authorization_servers": [issuer],
                        "scopes_supported": ["spatial:read", "spatial:write"]}

    async def metadata(request: Request):
        payload = specs.get(request.path_params["alias"])
        return JSONResponse(payload if payload else {"error": "Unknown resource"}, status_code=200 if payload else 404)

    async def health(request):
        return JSONResponse({"status": "ok", "mode": "oauth_resource_server"})

    routes.insert(0, Route("/.well-known/oauth-protected-resource/projects/{alias}/mcp", metadata))
    routes.insert(0, Route("/healthz", health))
    @asynccontextmanager
    async def lifespan(app):
        async with AsyncExitStack() as stack:
            for server in servers:
                await stack.enter_async_context(server.session_manager.run())
            yield
    app = Starlette(routes=routes, lifespan=lifespan)
    app.add_middleware(RequestLimits)
    app.add_middleware(TrustedHostMiddleware, allowed_hosts=[urlsplit(origin).hostname])
    return app


def main(argv=None):
    parser = argparse.ArgumentParser(description="Authenticated Spatial Collab gateway behind a TLS reverse proxy")
    parser.add_argument("--config", required=True)
    parser.add_argument("--port", type=int, default=8790)
    args = parser.parse_args(argv)
    import uvicorn
    settings = json.loads(Path(args.config).read_text(encoding="utf-8"))
    uvicorn.run(create_gateway(settings), host="127.0.0.1", port=args.port, proxy_headers=False, access_log=False)


if __name__ == "__main__":
    main()
