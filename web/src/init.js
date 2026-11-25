/* ================= başlat ================= */
(function () {
  const nav = $('#nav'); let g = null, h = '';
  for (const v of VIEWS) { if (v.grp !== g) { g = v.grp; h += `<div class="grp">${g}</div>`; }
    const ct = v.id === 'katalog' ? CAT.length : '';
    h += `<a href="/${v.id}" data-v="${v.id}">${ic(v.i)}<span>${v.t}</span>${ct !== '' ? `<span class="ct">${ct}</span>` : ''}</a>`; }
  nav.innerHTML = h;
  let lastW = innerWidth; window.addEventListener('resize', () => { clearTimeout(window._rz); window._rz = setTimeout(() => { if (Math.abs(innerWidth - lastW) > 40) { lastW = innerWidth; rerender(); } }, 220); });
  matchMedia('(prefers-color-scheme: dark)').addEventListener?.('change', () => rerender());
  initSiteLanguage();
  route();
  refreshConnection();
  connectLocal();
})();
