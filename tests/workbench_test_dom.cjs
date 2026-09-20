// Bounded DOM fixture for production workbench JavaScript. Not a browser engine.
const assert = require('node:assert/strict');
class Element {
  constructor(tag = 'div') {
    this.tag = tag; this.children = []; this.listeners = {}; this.listenerLists = {};
    this.style = {}; this.hidden = false; this.checked = false; this.disabled = false;
    this.classList = { add() {}, remove() {} }; this._text = ''; this._value = undefined;
  }
  set value(value) { this._value = String(value); }
  get value() { return this._value ?? (this.tag === 'select' ? this.options[0]?.value || '' : ''); }
  set textContent(value) { this._text = String(value); this.children = []; }
  get textContent() { return this._text + this.children.map(x => x.textContent ?? String(x)).join(' '); }
  set innerHTML(value) { throw new Error('Untrusted values must be rendered as text'); }
  append(...nodes) { for (const node of nodes) { node.parent = this; this.children.push(node); } }
  prepend(...nodes) { for (const node of nodes) node.parent = this; this.children.unshift(...nodes); }
  replaceChildren(...nodes) { this._text = ''; this.children = []; this._value = this.tag === 'select' ? undefined : this._value; this.append(...nodes); }
  remove() { if (this.parent) this.parent.children = this.parent.children.filter(node => node !== this); }
  addEventListener(name, callback) {
    (this.listenerLists[name] ||= []).push(callback);
    this.listeners[name] = async (event = { target: this }) => {
      for (const fn of this.listenerLists[name]) await fn(event);
    };
  }
  querySelector(selector) {
    const match = /^option\[value="([^"]+)"\]$/.exec(selector);
    if (match) return this.options.find(option => option.value === match[1]) || null;
    throw new Error(`Unsupported fixture selector: ${selector}`);
  }
  get options() { return this.children.filter(child => child.tag === 'option'); }
  getContext() { return new Proxy({}, { get: () => () => {} }); }
  getBoundingClientRect() { return { width: 800, height: 500, left: 0, top: 0 }; }
}

function createDOM(html) {
  const elements = new Map();
  for (const match of html.matchAll(/<([a-z]+)([^>]*\bid="([^"]+)"[^>]*)>/g)) {
    const node = new Element(match[1]), attributes = match[2];
    node.id = match[3];
    const value = /\bvalue="([^"]*)"/.exec(attributes);
    if (value) node.value = value[1];
    node.checked = /\bchecked\b/.test(attributes); node.hidden = /\bhidden\b/.test(attributes);
    elements.set(match[3], node);
  }
  for (const match of html.matchAll(/<select[^>]*\bid="([^"]+)"[^>]*>([\s\S]*?)<\/select>/g)) {
    const select = elements.get(match[1]);
    for (const item of match[2].matchAll(/<option[^>]*value="([^"]*)"[^>]*>([^<]*)<\/option>/g)) {
      const option = new Element('option'); option.value = item[1]; option.textContent = item[2]; select.append(option);
    }
  }
  const get = id => { assert(elements.has(id), `Missing real HTML id: ${id}`); return elements.get(id); };
  const all = () => {
    const seen = new Set();
    const visit = node => { if (seen.has(node)) return; seen.add(node); (node.children || []).forEach(visit); };
    elements.forEach(visit); return [...seen];
  };
  const document = { getElementById: get, createElement: tag => new Element(tag),
    createTextNode: text => ({ textContent: text }), querySelector: () => ({ content: 'fixture-csrf' }),
    querySelectorAll: selector => all().filter(node => selector.split(',').map(x => x.trim()).includes(node.tag)) };
  return { get, elements, document, all };
}
module.exports = { Element, createDOM };
