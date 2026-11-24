
/* The renderers supply Turkish copy. Keep the original DOM strings so switching
   language never rebuilds forms, clears results, or changes analysis values. */
const LOCALE_STORAGE_KEY = 'os-language';
const LOCALIZED_ATTRIBUTES = ['aria-label', 'title', 'placeholder', 'alt'];
const UNTRANSLATED = '[translate="no"], script, style, pre, code';
const originalText = new WeakMap();
const originalAttributes = new WeakMap();

const englishMessages = new Map();
const englishPatterns = [];
const messageSpace = value => value.replace(/\s+/g, ' ').trim();
const escapePattern = value => value.replace(/[.*+?^$(){}|[\]\\]/g, '\\$&');

for (const [source, target] of Object.entries(ENGLISH_MESSAGES)) {
  if (!/\{\d+\}/.test(source)) {
    englishMessages.set(source, target);
    continue;
  }
  const keys = [...source.matchAll(/\{(\d+)\}/g)].map(match => match[1]);
  const parts = source.split(/\{\d+\}/);
  const expression = parts.map(escapePattern).join('(.*?)');
  const specificity = parts.reduce((count, part) => count + part.length, 0);
  englishPatterns.push({match:new RegExp('^' + expression + '$', 'u'), keys, target, specificity});
}
englishPatterns.sort((a, b) => b.specificity - a.specificity);

let siteLanguage = 'tr';
try {
  if (localStorage.getItem(LOCALE_STORAGE_KEY) === 'en') siteLanguage = 'en';
} catch (_) {}

function translateMessage(value) {
  if (siteLanguage !== 'en' || !value?.trim()) return value;
  const source = messageSpace(value);
  let translated = englishMessages.get(source);
  if (translated === undefined) {
    for (const {match, keys, target} of englishPatterns) {
      const found = source.match(match);
      if (!found) continue;
      const args = new Map(keys.map((key, index) => [key, found[index + 1]]));
      translated = target.replace(/\{(\d+)\}/g, (_, key) => {
        const arg = args.get(key) ?? '';
        return englishMessages.get(arg) ?? arg;
      });
      break;
    }
  }
  if (translated === undefined && source.includes(' · ')) {
    translated = source.split(' · ').map(translateMessage).join(' · ');
    if (translated === source) return value;
  }
  if (translated === undefined) return value;
  return value.match(/^\s*/)[0] + translated + value.match(/\s*$/)[0];
}

function localizeTextNode(node) {
  if (node.parentElement?.closest(UNTRANSLATED)) return;
  const previous = originalText.get(node);
  const source = previous?.localized === node.data ? previous.source : node.data;
  const localized = translateMessage(source);
  originalText.set(node, {source, localized});
  if (node.data !== localized) node.data = localized;
}

function localizeAttributes(element) {
  if (element.closest(UNTRANSLATED)) return;
  if (element.tagName === 'OPTION' && !element.hasAttribute('value')) {
    element.setAttribute('value', element.value);
  }
  const entries = originalAttributes.get(element) || new Map();
  for (const name of LOCALIZED_ATTRIBUTES) {
    const value = element.getAttribute(name);
    if (value === null) {
      entries.delete(name);
      continue;
    }
    const previous = entries.get(name);
    const source = previous?.localized === value ? previous.source : value;
    const localized = translateMessage(source);
    entries.set(name, {source, localized});
    if (value !== localized) element.setAttribute(name, localized);
  }
  if (entries.size) originalAttributes.set(element, entries);
}

function localizeBranch(node) {
  if (node.nodeType === Node.TEXT_NODE) {
    localizeTextNode(node);
    return;
  }
  if (node.nodeType !== Node.ELEMENT_NODE || node.closest(UNTRANSLATED)) return;
  localizeAttributes(node);
  for (const child of node.childNodes) localizeBranch(child);
}

function applySiteLanguage(language) {
  siteLanguage = language;
  document.documentElement.lang = language;
  for (const button of document.querySelectorAll('[data-site-language]')) {
    button.setAttribute('aria-pressed', String(button.dataset.siteLanguage === language));
  }
  localizeBranch(document.body);
}

function initSiteLanguage() {
  applySiteLanguage(siteLanguage);
  const languageObserver = new MutationObserver(records => {
    for (const record of records) {
      if (record.type === 'childList') {
        for (const node of record.addedNodes) if (node.isConnected) localizeBranch(node);
      } else if (record.target.isConnected) {
        if (record.type === 'attributes') localizeAttributes(record.target);
        else localizeTextNode(record.target);
      }
    }
  });
  languageObserver.observe(document.body, {
    subtree:true, childList:true, characterData:true, attributes:true,
    attributeFilter:LOCALIZED_ATTRIBUTES,
  });
  document.querySelector('.language-picker').addEventListener('click', event => {
    const button = event.target.closest('[data-site-language]');
    if (!button || button.dataset.siteLanguage === siteLanguage) return;
    applySiteLanguage(button.dataset.siteLanguage);
    try { localStorage.setItem(LOCALE_STORAGE_KEY, siteLanguage); } catch (_) {}
  });
}
