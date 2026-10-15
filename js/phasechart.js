// Phase diagram chart (SVG): temperature across, pressure up on a log scale, the solid / liquid /
// gas / supercritical regions labelled in place, the triple and critical points, and a marker for
// where the jar is right now. Hovering reads off temperature, pressure and phase.
(function (CE) {
  'use strict';

  const NS = 'http://www.w3.org/2000/svg';
  const W = 640, H = 420, ML = 62, MR = 18, MT = 14, MB = 46;
  const PW = W - ML - MR, PH = H - MT - MB;
  const PHASE_LABEL = { solid: 'Solid', liquid: 'Liquid', gas: 'Gas', supercritical: 'Supercritical fluid' };

  function svg(tag, attrs, parent) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }

  function fmtP(p) {
    if (p >= 100) return Math.round(p).toLocaleString() + ' atm';
    if (p >= 0.01) return (+p.toPrecision(3)) + ' atm';
    const e = Math.floor(Math.log10(p)), m = p / Math.pow(10, e);
    return m.toFixed(1) + ' × 10' + sup(e) + ' atm';
  }
  function fmtT(t) { return (t >= 100 ? Math.round(t) : +t.toPrecision(3)) + ' K'; }
  function sup(n) { return String(n).replace(/-/g, '⁻').replace(/\d/g, function (d) { return '⁰¹²³⁴⁵⁶⁷⁸⁹'[+d]; }); }
  function decadeLabel(e) { return e === 0 ? '1' : e === 1 ? '10' : '10' + sup(e); }

  function niceStep(span) {
    const raw = span / 6, mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const f = raw / mag;
    return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
  }

  // d: CE.phase entry. Returns { el, setJar(T, P) }.
  function render(container, d) {
    const ph = CE.phase, g = ph.diagram(d);
    const lp0 = Math.log10(g.Pmin), lp1 = Math.log10(g.Pmax);
    const sx = function (T) { return ML + (T - g.Tmin) / (g.Tmax - g.Tmin) * PW; };
    const sy = function (P) { return MT + PH - (Math.log10(Math.max(P, 1e-300)) - lp0) / (lp1 - lp0) * PH; };
    const tOf = function (x) { return g.Tmin + (x - ML) / PW * (g.Tmax - g.Tmin); };
    const pOf = function (y) { return Math.pow(10, lp0 + (MT + PH - y) / PH * (lp1 - lp0)); };
    const clampP = function (P) { return Math.min(g.Pmax, Math.max(g.Pmin, P)); };

    const root = svg('svg', { viewBox: '0 0 ' + W + ' ' + H, class: 'phase-svg', role: 'img',
      'aria-label': 'Phase diagram of ' + d.name + ': temperature in kelvin across, pressure in atmospheres up on a logarithmic scale' });
    const clip = svg('clipPath', { id: 'phase-clip' }, svg('defs', {}, root));
    svg('rect', { x: ML, y: MT, width: PW, height: PH }, clip);

    // gridlines and axes
    const grid = svg('g', { class: 'ph-grid' }, root);
    const axes = svg('g', { class: 'ph-axis' }, root);
    const step = niceStep(g.Tmax - g.Tmin);
    for (let T = 0; T <= g.Tmax + 1e-9; T += step) {
      const x = sx(T);
      svg('line', { x1: x, x2: x, y1: MT, y2: MT + PH }, grid);
      const t = svg('text', { x: x, y: MT + PH + 18, 'text-anchor': 'middle' }, axes);
      t.textContent = Math.round(T * 100) / 100;
    }
    const dec = Math.max(1, Math.ceil((lp1 - lp0) / 8));
    for (let e = Math.ceil(lp0); e <= lp1; e += dec) {
      const y = sy(Math.pow(10, e));
      svg('line', { x1: ML, x2: ML + PW, y1: y, y2: y }, grid);
      const t = svg('text', { x: ML - 8, y: y + 4, 'text-anchor': 'end' }, axes);
      t.textContent = decadeLabel(e);
    }
    svg('line', { x1: ML, x2: ML + PW, y1: MT + PH, y2: MT + PH, class: 'ph-baseline' }, root);
    svg('line', { x1: ML, x2: ML, y1: MT, y2: MT + PH, class: 'ph-baseline' }, root);
    const xl = svg('text', { x: ML + PW / 2, y: H - 6, 'text-anchor': 'middle', class: 'ph-axis-title' }, root);
    xl.textContent = 'Temperature (K)';
    const yl = svg('text', { x: 14, y: MT + PH / 2, 'text-anchor': 'middle', class: 'ph-axis-title',
      transform: 'rotate(-90 14 ' + (MT + PH / 2) + ')' }, root);
    yl.textContent = 'Pressure (atm, log scale)';

    // regions: gas is the background; solid, liquid and supercritical are drawn over it
    const regions = svg('g', { 'clip-path': 'url(#phase-clip)' }, root);
    svg('rect', { x: ML, y: MT, width: PW, height: PH, class: 'ph-gas' }, regions);
    const lines = {};
    g.lines.forEach(function (l) { lines[l.kind] = l.pts; });
    const pt = function (p) { return sx(p[0]).toFixed(1) + ',' + sy(clampP(p[1])).toFixed(1); };
    const mel = lines.melting, vap = lines.vaporisation, sub = lines.sublimation || [];
    const top = sy(g.Pmax), bottom = sy(g.Pmin), left = sx(0), right = sx(g.Tmax);
    const melEnd = mel[mel.length - 1];
    // solid: left of the melting line, above the sublimation line
    let solid = [left + ',' + bottom];
    for (let k = sub.length - 1; k >= 0; k--) solid.push(pt(sub[k]));
    mel.forEach(function (p) { solid.push(pt(p)); });
    solid.push(sx(melEnd[0]).toFixed(1) + ',' + top, left + ',' + top);
    svg('polygon', { points: solid.join(' '), class: 'ph-solid' }, regions);
    // liquid: right of the melting line, above the vaporisation curve, below/left of the critical point
    const tc = sx(d.Tc), pc = sy(d.Pc);
    let liquid = mel.map(pt);
    liquid.push(sx(melEnd[0]).toFixed(1) + ',' + top, tc + ',' + top, tc + ',' + pc);
    for (let k = vap.length - 1; k >= 0; k--) liquid.push(pt(vap[k]));
    svg('polygon', { points: liquid.join(' '), class: 'ph-liquid' }, regions);
    svg('rect', { x: tc, y: top, width: Math.max(0, right - tc), height: Math.max(0, pc - top), class: 'ph-sc' }, regions);

    // reference line at 1 atm
    if (g.Pmin < 1 && g.Pmax > 1) {
      svg('line', { x1: ML, x2: ML + PW, y1: sy(1), y2: sy(1), class: 'ph-ref' }, regions);
      const t = svg('text', { x: ML + PW - 4, y: sy(1) - 5, 'text-anchor': 'end', class: 'ph-ref-label' }, regions);
      t.textContent = '1 atm';
    }
    // the supercritical edges are not phase transitions: thin, dotted
    svg('path', { d: 'M' + tc + ',' + top + 'V' + pc + 'H' + right, class: 'ph-sc-edge' }, regions);
    // phase boundaries
    g.lines.forEach(function (l) {
      svg('polyline', { points: l.pts.map(pt).join(' '), class: 'ph-line' }, regions);
    });

    // region labels, placed inside each region
    function label(text, T, P, cls) {
      const t = svg('text', { x: sx(T), y: sy(P), 'text-anchor': 'middle', class: 'ph-region ' + (cls || '') }, root);
      t.textContent = text;
    }
    if (!d.helium) {
      const Ps = Math.pow(10, Math.log10(d.Pt) + 0.6 * (lp1 - Math.log10(d.Pt)));
      label('SOLID', Math.max(g.Tmax * 0.05, ph.meltT(d, Ps) * 0.5), Ps);
      const Pl = Math.min(g.Pmax / 3, d.Pc * 4);
      label('LIQUID', (Math.max(ph.meltT(d, Pl), d.Tt) + d.Tc) / 2, Pl);
      const Tg = d.Tt + (d.Tc - d.Tt) * 0.7;
      label('GAS', Tg, Math.max(g.Pmin * 30, ph.vaporP(d, Tg) / 300));
    } else {
      label('SOLID', d.Tc * 0.5, g.Pmax / 3);
      label('LIQUID', d.Tc * 0.45, 3);
      label('GAS', d.Tc * 0.75, g.Pmin * 30);
    }
    label('SUPERCRITICAL', (d.Tc + g.Tmax) / 2, Math.min(g.Pmax / 3, Math.sqrt(d.Pc * g.Pmax)), 'ph-region-sm');

    // triple and critical points
    function point(T, P, text, anchor, dy) {
      svg('circle', { cx: sx(T), cy: sy(P), r: 4.5, class: 'ph-point' }, root);
      const t = svg('text', { x: sx(T) + (anchor === 'end' ? -9 : 9), y: sy(P) + (dy || 4), 'text-anchor': anchor || 'start', class: 'ph-point-label' }, root);
      t.textContent = text;
    }
    if (!d.helium) point(d.Tt, d.Pt, 'Triple point', 'start', 14);
    point(d.Tc, d.Pc, 'Critical point', 'start', 18);

    // the jar right now
    const jar = svg('g', { class: 'ph-jar' }, root);
    const jarRing = svg('circle', { r: 7 }, jar);
    const jarDot = svg('circle', { r: 2.5, class: 'ph-jar-dot' }, jar);
    const jarText = svg('text', { class: 'ph-jar-label' }, jar);

    // hover read-out
    const hover = svg('g', { class: 'ph-hover', visibility: 'hidden' }, root);
    const hx = svg('line', { y1: MT, y2: MT + PH }, hover);
    const hy = svg('line', { x1: ML, x2: ML + PW }, hover);
    const tipBg = svg('rect', { rx: 6, class: 'ph-tip' }, hover);
    const tip = svg('text', { class: 'ph-tip-text' }, hover);
    const hit = svg('rect', { x: ML, y: MT, width: PW, height: PH, fill: 'transparent' }, root);
    function showAt(x, y) {
      const T = tOf(x), P = pOf(y);
      hx.setAttribute('x1', x); hx.setAttribute('x2', x);
      hy.setAttribute('y1', y); hy.setAttribute('y2', y);
      tip.textContent = fmtT(T) + ' · ' + fmtP(P) + ' → ' + PHASE_LABEL[ph.phaseAt(d, T, P)];
      const tw = tip.getComputedTextLength ? tip.getComputedTextLength() : 180;
      let tx = x + 12, ty = y - 12;
      if (tx + tw + 16 > ML + PW) tx = x - tw - 20;
      if (ty < MT + 20) ty = y + 26;
      tip.setAttribute('x', tx + 6); tip.setAttribute('y', ty);
      tipBg.setAttribute('x', tx); tipBg.setAttribute('y', ty - 15);
      tipBg.setAttribute('width', tw + 12); tipBg.setAttribute('height', 22);
      hover.setAttribute('visibility', 'visible');
    }
    hit.addEventListener('pointermove', function (e) {
      const r = root.getBoundingClientRect();
      showAt((e.clientX - r.left) / r.width * W, (e.clientY - r.top) / r.height * H);
    });
    hit.addEventListener('pointerleave', function () { hover.setAttribute('visibility', 'hidden'); });

    container.appendChild(root);

    function setJar(T, P) {
      const vacuum = !(P > g.Pmin);
      const x = Math.min(ML + PW, Math.max(ML, sx(T))), y = sy(clampP(vacuum ? g.Pmin : P));
      jarRing.setAttribute('cx', x); jarRing.setAttribute('cy', y);
      jarDot.setAttribute('cx', x); jarDot.setAttribute('cy', y);
      const phase = PHASE_LABEL[ph.phaseAt(d, T, vacuum ? g.Pmin : P)];
      jarText.textContent = 'Jar now: ' + fmtT(T) + ', ' + (vacuum ? 'vacuum' : fmtP(P)) + ' → ' + phase.toLowerCase();
      const right = x > ML + PW * 0.6;
      jarText.setAttribute('text-anchor', right ? 'end' : 'start');
      jarText.setAttribute('x', x + (right ? -12 : 12));
      jarText.setAttribute('y', y < MT + 24 ? y + 18 : y - 10);
      return phase;
    }
    return { el: root, setJar: setJar };
  }

  // Key values as a table (the chart's readable twin).
  function table(d) {
    const est = function (k) { return d.est.indexOf(k) >= 0 ? ' (est.)' : ''; };
    const rows = [];
    if (!d.helium) rows.push(['Triple point', fmtT(d.Tt) + est('Tt') + ', ' + fmtP(d.Pt) + est('Pt')]);
    rows.push(['Critical point', fmtT(d.Tc) + est('Tc') + ', ' + fmtP(d.Pc) + est('Pc')]);
    if (d.Ts) rows.push(['Sublimes at 1 atm', fmtT(d.Ts)]);
    else if (d.Tm) rows.push(['Melts at 1 atm', fmtT(d.Tm)]);
    if (d.Tb) rows.push(['Boils at 1 atm', fmtT(d.Tb)]);
    const t = document.createElement('table');
    t.className = 'ph-table';
    rows.forEach(function (r) {
      const tr = t.insertRow();
      const th = document.createElement('th'); th.scope = 'row'; th.textContent = r[0];
      tr.appendChild(th);
      tr.insertCell().textContent = r[1];
    });
    return t;
  }

  CE.PhaseChart = { render: render, table: table, PHASE_LABEL: PHASE_LABEL };
})(window.CE = window.CE || {});
