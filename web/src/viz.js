/* ================= görselleştirme ================= */
const W_OF = el => Math.max(320, Math.floor((el && el.clientWidth) || 900));

/* --- şablon genel izi: cetvel, GC, varyant, çift glifleri --- */
function overviewTrack(width, opts) {
  const {tpl=TPL, off=OFF, pairs=[], sel=null, showVars=true, showGC=true, region=null, compact=false, win=null} = opts;
  const L = tpl.length, lo = win ? Math.max(0, win[0]) : 0, hi = win ? Math.min(L, win[1]) : L;
  const m = {l:58, r:14}, iw = width - m.l - m.r; const x = i => m.l + ((Math.max(lo, Math.min(hi, i)) - lo) / (hi - lo)) * iw;
  let y = 8; const out = [];
  const span = hi - lo, maxT = Math.max(3, Math.floor(iw / 92)); const step = [10,20,25,50,100,200,250,500,1000,2000,5000,10000].find(s => span / s <= maxT) || 20000;
  out.push(`<line x1="${m.l}" x2="${width-m.r}" y1="${y+14}" y2="${y+14}" stroke="${css('--line2')}"/>`);
  for (let i = Math.ceil(lo / step) * step; i <= hi; i += step) { const X = x(i);
    out.push(`<line x1="${X}" x2="${X}" y1="${y+10}" y2="${y+18}" stroke="${css('--line2')}"/><text x="${X}" y="${y+6}" class="axm" text-anchor="middle">${off ? (off + i + 1).toLocaleString('tr-TR') : (i + 1).toLocaleString('tr-TR')}</text>`); }
  out.push(`<text x="4" y="${y+17}" class="ax">${off ? CHR : 'konum'}</text>`); y += 26;
  if (region) out.push(`<rect x="${x(region[0])}" y="${y-4}" width="${Math.max(2, x(region[1]+1)-x(region[0]))}" height="2000" fill="${css('--accent')}" opacity=".06"/>`);
  if (showGC) {
    const w = span > 600 ? 25 : 15, pts = []; const stp = Math.max(1, Math.floor(span / 500));
    for (let i = Math.max(0, lo); i + w <= Math.min(L, hi + w); i += stp) { if (i + w > L) break; let n = 0; for (let k = i; k < i + w; k++) if (tpl[k] === 'G' || tpl[k] === 'C') n++; pts.push([x(i + w/2), y + 30 - (n/w) * 30]); }
    out.push(`<text x="4" y="${y+19}" class="ax">GC %</text><line x1="${m.l}" x2="${width-m.r}" y1="${y+15}" y2="${y+15}" stroke="${css('--grid')}" stroke-dasharray="2 3"/>`);
    if (pts.length > 1) out.push(`<path d="M${pts.map(p => p[0].toFixed(1)+' '+p[1].toFixed(1)).join('L')}L${pts.at(-1)[0]} ${y+30}L${pts[0][0]} ${y+30}Z" fill="${css('--accent')}" opacity=".14"/><path d="M${pts.map(p => p[0].toFixed(1)+' '+p[1].toFixed(1)).join('L')}" fill="none" stroke="${css('--accent')}" stroke-width="1.3"/>`);
    y += 40;
  }
  if (showVars && off === OFF && tpl === TPL) {
    const H = 34; out.push(`<text x="4" y="${y+22}" class="ax">gnomAD</text><line x1="${m.l}" x2="${width-m.r}" y1="${y+H}" y2="${y+H}" stroke="${css('--line')}"/>`);
    for (const v of VARS) { if (v.idx < lo || v.idx > hi) continue; const X = x(v.idx + .5), hh = Math.max(4, Math.min(H-2, (Math.log10(Math.max(v.af||1e-6, 1e-6)) + 6.2) / 6 * H));
      out.push(`<line x1="${X}" x2="${X}" y1="${y+H}" y2="${y+H-hh}" stroke="${varColor(v.af)}" stroke-width="${v.af>=0.001?1.8:1}" opacity="${v.af>=0.001?1:.55}"/>${v.af>=0.001?`<circle cx="${X}" cy="${y+H-hh}" r="2.8" fill="${varColor(v.af)}"/>`:''}<rect x="${X-2.5}" y="${y}" width="5" height="${H}" fill="transparent" data-tip="${esc(varTip(v))}"/>`); }
    y += H + 10;
  }
  const rowH = compact ? 14 : 17;
  pairs.forEach((p, k) => {
    const isSel = sel && p.key === sel; const yy = y + k * rowH;
    if (isSel) out.push(`<rect x="0" y="${yy-2}" width="${width}" height="${rowH}" fill="${css('--accent-soft')}" opacity=".8"/>`);
    out.push(`<text x="4" y="${yy+9}" class="axm" ${isSel?`style="fill:${css('--accent2')};font-weight:700"`:''}>${esc(p.short || p.label).slice(0, 10)}</text>`);
    out.push(`<line x1="${x(p.fStart)}" x2="${x(p.rEnd+1)}" y1="${yy+5}" y2="${yy+5}" stroke="${css('--ink3')}" stroke-width="1.2" opacity=".6"/>`);
    const fx = x(p.fStart), fw = Math.max(3, x(p.fEnd+1) - fx), rx = x(p.rStart), rw = Math.max(3, x(p.rEnd+1) - rx);
    out.push(`<path d="M${fx} ${yy+1}h${Math.max(0,fw-4)}l4 4l-4 4h-${Math.max(0,fw-4)}z" fill="${css('--fwd')}"/><path d="M${rx+rw} ${yy+1}h-${Math.max(0,rw-4)}l-4 4l4 4h${Math.max(0,rw-4)}z" fill="${css('--rev')}"/>`);
    if (!compact && iw / span > 1.2) out.push(`<text x="${x(p.rEnd+1)+6}" y="${yy+9}" class="axm" style="font-size:9.5px">${p.prod} bp</text>`);
    out.push(`<rect x="0" y="${yy-2}" width="${width}" height="${rowH}" fill="transparent" data-pair="${esc(p.key)}" style="cursor:pointer" data-tip="${esc(`<b>${esc(p.label)}</b><br>${p.prod} bp · F ${gposAny(p,p.fStart)} → R ${gposAny(p,p.rEnd)}<br>Tm ${nf(p.tmF,2)} / ${nf(p.tmR,2)} °C`)}"/>`);
  });
  y += pairs.length * rowH + 6;
  return `<svg width="${width}" height="${y}" viewBox="0 0 ${width} ${y}" role="img" aria-label="Şablon üzerindeki primer çiftleri">${out.join('')}</svg>`;
}
const gposAny = (p, i) => p.tplOff ? (p.tplOff + i + 1).toLocaleString('tr-TR') : (i + 1).toLocaleString('tr-TR');
function varTip(v) {
  const pops = Object.entries(v.pop || {}).sort((a, b) => b[1] - a[1]).slice(0, 3).filter(x => x[1] > 0).map(([k, a]) => `${k} ${afFmt(a)}`).join(' · ');
  return `<b>${esc(v.rsid || 'varyant')}</b> <span class="m">${esc(v.ref)}&gt;${esc(v.alt)}</span><br>chr22:${v.pos.toLocaleString('tr-TR')} · ${vtype(v)}<br>AF ${afFmt(v.af)} · AN ${ni(v.an)}${pops ? '<br>' + pops : ''}`;
}

/* --- baz düzeyinde dizi haritası (iki iplik, primer okları, varyantlar) --- */
function seqMap(width, pm, opts={}) {
  const pad = opts.pad ?? 24; const tpl = pm.tpl;
  const lo = Math.max(0, (opts.lo ?? pm.fStart) - pad), hi = Math.min(tpl.length - 1, (opts.hi ?? pm.rEnd) + pad);
  const CW = 11.2, LM = 116, RM = 14;
  let N = Math.floor((width - LM - RM) / CW); N = Math.max(20, Math.min(100, N - N % 10));
  const vlist = pm.hasVar ? varsIn(lo, hi) : [];
  const vmap = {}; for (const v of vlist) for (let k = 0; k < Math.max(1, v.ref.length === v.alt.length ? 1 : v.ref.length); k++) (vmap[v.idx + k] = vmap[v.idx + k] || []).push(v);
  const probe = opts.probe; // {start,end,strand}
  const rowsH = 106; const rows = Math.ceil((hi - lo + 1) / N); const H = rows * rowsH + 8;
  const fC = css('--fwd'), rC = css('--rev'), pC = css('--probe'), amp = css('--amp'), i3 = css('--ink3'), i2 = css('--ink2');
  const nuc = {A:css('--nA'), C:css('--nC'), G:css('--nG'), T:css('--nT'), N:i3};
  const o = [];
  for (let r = 0; r < rows; r++) {
    const s0 = lo + r * N, s1 = Math.min(hi, s0 + N - 1), Y = 8 + r * rowsH;
    const X = i => LM + (i - s0) * CW;
    const yF = Y + 4, yV = Y + 22, yP = Y + 52, yM = Y + 70, yR = Y + 80;
    // amplikon zemini
    const a0 = Math.max(s0, pm.fStart), a1 = Math.min(s1, pm.rEnd);
    if (a0 <= a1) o.push(`<rect x="${X(a0)-1}" y="${yP-13}" width="${(a1-a0+1)*CW+2}" height="${yM-yP+18}" rx="3" fill="${amp}"/>`);
    // 10'luk işaretler
    for (let i = s0; i <= s1; i++) if ((i + 1) % 10 === 0) o.push(`<line x1="${X(i)+CW/2}" x2="${X(i)+CW/2}" y1="${yV+18}" y2="${yV+22}" stroke="${i3}"/><text x="${X(i)+CW/2}" y="${yV+16}" class="axm" text-anchor="middle" style="font-size:8.5px">${pm.tplOff ? ((pm.tplOff+i+1)%1000).toString().padStart(3,'0') : i+1}</text>`);
    o.push(`<text x="${LM-8}" y="${yP}" class="axm" text-anchor="end">${pm.tplOff ? gposAny(pm, s0) : s0+1}</text>`);
    o.push(`<text x="4" y="${yP}" class="ax" style="font-weight:600">+5′</text><text x="4" y="${yM}" class="ax" style="font-weight:600">−3′</text>`);
    // primerler
    const bar = (st, en, y, col, dir, lab) => { const b0 = Math.max(st, s0), b1 = Math.min(en, s1); if (b0 > b1) return;
      const x0 = X(b0), x1 = X(b1) + CW; const head = (dir > 0 && b1 === en) || (dir < 0 && b0 === st);
      if (dir > 0) o.push(`<path d="M${x0} ${y}h${x1-x0-(head?6:0)}${head?`l6 5l-6 5`:'v10'}h-${x1-x0-(head?6:0)}z" fill="${col}" opacity=".9"/>`);
      else o.push(`<path d="M${x1} ${y}h-${x1-x0-(head?6:0)}${head?`l-6 5l6 5`:'v10'}h${x1-x0-(head?6:0)}z" fill="${col}" opacity=".9"/>`);
      if ((dir > 0 && b0 === st) || (dir < 0 && b1 === en)) o.push(`<text x="${dir>0?x0+4:x1-4}" y="${y+8.2}" text-anchor="${dir>0?'start':'end'}" style="font:700 8.5px var(--sans);fill:${css('--panel')}">${lab}</text>`); };
    bar(pm.fStart, pm.fEnd, yF, fC, 1, 'F  5′→3′');
    bar(pm.rStart, pm.rEnd, yR, rC, -1, '3′←5′  R');
    if (probe) bar(probe.start, probe.end, probe.strand === '-' ? yR + 14 : yF + 12, pC, probe.strand === '-' ? -1 : 1, 'PROB');
    // bazlar
    for (let i = s0; i <= s1; i++) {
      const b = tpl[i], cb = OS.COMP[b] || 'N'; const inF = i >= pm.fStart && i <= pm.fEnd, inR = i >= pm.rStart && i <= pm.rEnd;
      const vv = vmap[i]; const xc = X(i) + CW/2;
      o.push(`<text x="${xc}" y="${yP}" text-anchor="middle" class="seqt" style="fill:${inF ? nuc[b] : (vv ? nuc[b] : i2)};font-weight:${inF?700:400};opacity:${inF||inR||vv?1:.72}">${b}</text>`);
      o.push(`<text x="${xc}" y="${yM}" text-anchor="middle" class="seqt" style="fill:${inR ? nuc[cb] : i3};font-weight:${inR?700:400};opacity:${inR?1:.6}">${cb}</text>`);
      if (vv) { const top = vv.reduce((a, v) => (v.af||0) > (a.af||0) ? v : a); const col = varColor(top.af);
        const big = top.af >= 0.001; o.push(`<path d="M${xc-4} ${yV+2}h8l-4 ${big?8:6}z" fill="${col}" opacity="${big?1:.7}"/>`);
        if (big) o.push(`<line x1="${xc}" x2="${xc}" y1="${yV+10}" y2="${yP-12}" stroke="${col}" stroke-width="1" stroke-dasharray="2 2"/>`);
        o.push(`<rect x="${X(i)}" y="${yV}" width="${CW}" height="${yP-yV+2}" fill="transparent" data-tip="${esc(vv.map(varTip).join('<hr style=&quot;border:0;border-top:1px solid rgba(255,255,255,.2);margin:5px 0&quot;>'))}"/>`); }
    }
  }
  return `<svg width="${width}" height="${H}" viewBox="0 0 ${width} ${H}" role="img" aria-label="Baz düzeyinde dizi haritası">${o.join('')}</svg>`;
}

/* --- Primer-BLAST tarzı bağlanma yeri hizalaması --- */
function siteAlignHTML(primer, site, anchor=12) {
  const L = primer.length; let top = '', bot = '';
  for (let i = 0; i < L; i++) { const az = i >= L - anchor; const p = primer[i], s = site[i] || ' ';
    top += `<span class="${az?'anch':''}">${p}</span>`;
    bot += p === s ? `<span class="${az?'anch':''}">.</span>` : `<span class="mm" title="${i+1}. konum: ${p}→${s}">${s}</span>`; }
  return {top, bot};
}

/* --- IDT OligoAnalyzer tarzı dimer gösterimi --- */
function dimerHTML(al, opts={}) {
  if (!al) return '<div class="faint small">Kararlı yapı bulunamadı</div>';
  const [t, m, b] = al; const w = Math.max(t.length, m.length, b.length);
  const col = ch => ch === ' ' ? ' ' : ch === '-' ? '<span class="lbl">-</span>' : `<span class="n${ch}">${ch}</span>`;
  let T = '', M = '', B = '';
  for (let i = 0; i < w; i++) { const tc = t[i] || ' ', mc = m[i] || ' ', bc = b[i] || ' ';
    T += col(tc); B += col(bc);
    if (mc === '|') { const gc = (tc === 'G' || tc === 'C'); M += `<span class="${gc?'bd':'bdw'}">${gc?'┃':'│'}</span>`; }
    else if (mc === 'X') M += '<span class="bx">×</span>'; else if (mc === 'B') M += '<span class="bx">^</span>'; else M += ' '; }
  return `<div class="aln"><span class="lbl">5′ </span>${T}<span class="lbl"> 3′</span>\n<span class="lbl">   </span>${M}\n<span class="lbl">3′ </span>${B}<span class="lbl"> 5′</span></div>`;
}

/* --- sap-ilmek (hairpin) çizimi nokta-parantezden --- */
function hairpinSVG(seq, db, opts={}) {
  const n = seq.length; const pairs = []; const st = [];
  for (let i = 0; i < n; i++) { if (db[i] === '(') st.push(i); else if (db[i] === ')') { const j = st.pop(); if (j != null) pairs.push([j, i]); } }
  const nuc = {A:css('--nA'), C:css('--nC'), G:css('--nG'), T:css('--nT')};
  const W = Math.max(260, n * 13 + 40), H = 150;
  if (!pairs.length) { // doğrusal + yay diyagramı boş
    const o = seq.split('').map((c, i) => `<text x="${20+i*13}" y="80" class="seqt" text-anchor="middle" style="fill:${nuc[c]}">${c}</text>`).join('');
    return `<svg width="${W}" height="110" viewBox="0 0 ${W} 110">${o}<text x="20" y="40" class="lab">Eşleşmiş baz yok — doğrusal, yapısız (MFE 0)</text><text x="8" y="80" class="ax">5′</text><text x="${24+n*13}" y="80" class="ax">3′</text></svg>`;
  }
  // yay diyagramı: bazlar bir çizgide, eşleşmeler yay
  const X = i => 20 + i * 13; const base = 118; const o = [];
  for (const [i, j] of pairs) { const r = (X(j) - X(i)) / 2; o.push(`<path d="M${X(i)} ${base-12} A ${r} ${Math.min(r, 90)} 0 0 1 ${X(j)} ${base-12}" fill="none" stroke="${css('--accent')}" stroke-width="1.6" opacity=".85"/>`); }
  seq.split('').forEach((c, i) => { const paired = db[i] !== '.'; o.push(`<text x="${X(i)}" y="${base}" class="seqt" text-anchor="middle" style="fill:${nuc[c]};font-weight:${paired?700:400}">${c}</text><text x="${X(i)}" y="${base+14}" class="axm" text-anchor="middle" style="font-size:9px">${db[i]}</text>`); });
  return `<svg width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">${o.join('')}<text x="6" y="${base}" class="ax">5′</text><text x="${X(n-1)+12}" y="${base}" class="ax">3′</text></svg>`;
}
/* sap-ilmek klasik çizimi (yığın + ilmek yayı); hp: hairpinStem sonucu */
function stemLoopSVG(seq, hp) {
  if (!hp) return '';
  const nuc = {A:css('--nA'), C:css('--nC'), G:css('--nG'), T:css('--nT')};
  const L = hp.helix_len, i = hp.i, j = hp.j, loopN = j - i - 1, cx = 150, gap = 19, o = [];
  const R = Math.max(26, loopN * 5.5), y0 = 28 + R * 2;          // en içteki çiftin satırı
  for (let k = 0; k < L; k++) { const a = i - k, b = j + k, y = y0 + k * gap;
    o.push(`<text x="${cx-18}" y="${y}" class="seqt" text-anchor="middle" style="fill:${nuc[seq[a]]};font-weight:700">${seq[a]}</text><text x="${cx+18}" y="${y}" class="seqt" text-anchor="middle" style="fill:${nuc[seq[b]]};font-weight:700">${seq[b]}</text><line x1="${cx-9}" x2="${cx+9}" y1="${y-4}" y2="${y-4}" stroke="${css('--accent')}" stroke-width="${'GC'.includes(seq[a])?2.6:1.3}"/>`); }
  const yc = y0 - 14;
  for (let k = 0; k < loopN; k++) { const ang = Math.PI - (k + 1) / (loopN + 1) * Math.PI;
    o.push(`<text x="${(cx + (R+4) * Math.cos(ang)).toFixed(1)}" y="${(yc - R * Math.sin(ang) + 4).toFixed(1)}" class="seqt" text-anchor="middle" style="fill:${nuc[seq[i+1+k]]}">${seq[i+1+k]}</text>`); }
  const yb = y0 + (L - 1) * gap + 22, t5 = seq.slice(0, i - L + 1), t3 = seq.slice(j + L);
  o.push(`<text x="${cx-26}" y="${yb}" class="seqt" text-anchor="end" style="fill:${css('--ink3')}">5′ ${t5.slice(-7)}</text><text x="${cx+26}" y="${yb}" class="seqt" style="fill:${css('--ink3')}">${t3.slice(0,7)} 3′</text>`);
  o.push(`<line x1="${cx-18}" x2="${cx-24}" y1="${yb-14}" y2="${yb-6}" stroke="${css('--ink3')}"/><line x1="${cx+18}" x2="${cx+24}" y1="${yb-14}" y2="${yb-6}" stroke="${css('--ink3')}"/>`);
  const H = yb + 10;
  return `<svg width="300" height="${H}" viewBox="0 0 300 ${H}" role="img" aria-label="Sap-ilmek yapısı">${o.join('')}</svg>`;
}
/* --- genel çizgi grafiği --- */
function lineChart(width, height, series, o={}) {
  const m = {l:44, r:14, t:12, b:30}; const iw = width - m.l - m.r, ih = height - m.t - m.b;
  const xs = series.flatMap(s => s.pts.map(p => p[0])), ys = series.flatMap(s => s.pts.map(p => p[1]));
  const x0 = o.x0 ?? Math.min(...xs), x1 = o.x1 ?? Math.max(...xs), y0 = o.y0 ?? Math.min(0, ...ys), y1 = o.y1 ?? Math.max(...ys) * 1.06;
  const X = v => m.l + (v - x0) / (x1 - x0) * iw, Y = v => m.t + ih - (v - y0) / (y1 - y0 || 1) * ih;
  const g = []; const xt = o.xticks || niceTicks(x0, x1, 6), yt = o.yticks || niceTicks(y0, y1, 4);
  for (const t of yt) g.push(`<line x1="${m.l}" x2="${width-m.r}" y1="${Y(t)}" y2="${Y(t)}" stroke="${css('--grid')}"/><text x="${m.l-6}" y="${Y(t)+3}" class="axm" text-anchor="end">${o.yfmt ? o.yfmt(t) : t}</text>`);
  for (const t of xt) g.push(`<text x="${X(t)}" y="${height-m.b+14}" class="axm" text-anchor="middle">${o.xfmt ? o.xfmt(t) : t}</text>`);
  g.push(`<line x1="${m.l}" x2="${width-m.r}" y1="${Y(y0)}" y2="${Y(y0)}" stroke="${css('--line2')}"/>`);
  if (o.xlabel) g.push(`<text x="${width-m.r}" y="${height-4}" class="ax" text-anchor="end">${o.xlabel}</text>`);
  if (o.ylabel) g.push(`<text x="${m.l}" y="${m.t-2}" class="ax">${o.ylabel}</text>`);
  for (const v of o.vlines || []) g.push(`<line x1="${X(v.x)}" x2="${X(v.x)}" y1="${m.t}" y2="${m.t+ih}" stroke="${v.c || css('--ink3')}" stroke-dasharray="${v.dash || '3 3'}" opacity="${v.op ?? .9}"/>${v.t ? `<text x="${X(v.x)+(v.end?-4:4)}" y="${m.t+10+(v.dy||0)}" class="axm" text-anchor="${v.end?'end':'start'}" style="fill:${v.c||css('--ink3')}">${v.t}</text>` : ''}`);
  for (const b of o.bands || []) g.push(`<rect x="${X(b.a)}" y="${m.t}" width="${X(b.b)-X(b.a)}" height="${ih}" fill="${b.c}" opacity="${b.op ?? .08}"/>`);
  for (const s of series) { const d = s.pts.map((p, i) => (i ? 'L' : 'M') + X(p[0]).toFixed(1) + ' ' + Y(p[1]).toFixed(1)).join('');
    if (s.area) g.push(`<path d="${d}L${X(s.pts.at(-1)[0])} ${Y(y0)}L${X(s.pts[0][0])} ${Y(y0)}Z" fill="${s.c}" opacity="${s.areaOp ?? .12}"/>`);
    g.push(`<path d="${d}" fill="none" stroke="${s.c}" stroke-width="${s.w || 2}" opacity="${s.op ?? 1}" ${s.dash ? `stroke-dasharray="${s.dash}"` : ''}/>`); }
  for (const p of o.points || []) g.push(`<circle cx="${X(p.x)}" cy="${Y(p.y)}" r="${p.r || 3.5}" fill="${p.c}" stroke="${css('--panel')}" stroke-width="1.5" ${p.tip ? `data-tip="${esc(p.tip)}"` : ''}/>${p.t ? `<text x="${X(p.x)+6}" y="${Y(p.y)-6}" class="axm" style="fill:${p.c}">${p.t}</text>` : ''}`);
  if (o.hover) { const hx = []; const n = series[0].pts.length; const stp = Math.max(1, Math.floor(n / 120));
    for (let i = 0; i < n; i += stp) { const px = series[0].pts[i][0]; const nx = series[0].pts[Math.min(n-1, i+stp)][0];
      hx.push(`<rect x="${X(px)}" y="${m.t}" width="${Math.max(1, X(nx)-X(px))}" height="${ih}" fill="transparent" data-tip="${esc(o.hover(i))}"/>`); } g.push(hx.join('')); }
  return `<svg width="${width}" height="${height}" viewBox="0 0 ${width} ${height}" role="img">${g.join('')}</svg>`;
}
function niceTicks(a, b, n) { const span = b - a || 1; const raw = span / n; const p = Math.pow(10, Math.floor(Math.log10(raw))); const f = raw / p; const st = (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p; const out = []; for (let t = Math.ceil(a / st) * st; t <= b + 1e-9; t += st) out.push(+t.toFixed(10)); return out; }

/* --- yatay çubuk (değer + eşik bantları) --- */
function gaugeHTML(val, lo, hi, bands, fmt=(v)=>nf(v,2), unit='') {
  const pc = v => Math.max(0, Math.min(100, (v - lo) / (hi - lo) * 100));
  const bs = bands.map(b => `<i style="position:absolute;left:${pc(b.a)}%;width:${pc(b.b)-pc(b.a)}%;top:0;bottom:0;background:${b.c};opacity:${b.op??.35}"></i>`).join('');
  return `<div style="position:relative;height:8px;border-radius:5px;background:var(--sunk);overflow:hidden;margin:8px 0 4px">${bs}</div><div style="position:relative;height:0"><span style="position:absolute;left:calc(${pc(val)}% - 6px);top:-15px;width:12px;height:12px;border-radius:50%;background:var(--ink);border:2px solid var(--panel)"></span></div><div class="row tiny faint" style="justify-content:space-between;margin-top:6px"><span class="mono">${fmt(lo)}${unit}</span><span class="mono">${fmt(hi)}${unit}</span></div>`;
}
