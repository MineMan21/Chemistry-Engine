// DOM helpers: periodic table, element card, contents list, inspector.
(function (CE) {
  'use strict';

  const E = CE.elements;
  const UI = {};

  function el(tag, cls, text) {
    const e = document.createElement(tag);
    if (cls) e.className = cls;
    if (text !== undefined) e.textContent = text;
    return e;
  }

  // Standard 18-column layout with the f-block underneath.
  UI.buildPeriodicTable = function (container, onSelect) {
    container.innerHTML = '';
    const cells = {};
    for (let z = 1; z <= 118; z++) {
      const e = E[z];
      const b = el('button', 'cell cat-' + e.category);
      b.innerHTML = '<span class="z">' + z + '</span><span>' + e.symbol + '</span>';
      b.title = e.name + ' (' + e.symbol + '), ' + e.categoryLabel;
      if (e.fIndex >= 0) {
        b.style.gridRow = String(e.period + 3); // rows 9, 10
        b.style.gridColumn = String(e.fIndex + 3);
      } else {
        b.style.gridRow = String(e.period);
        b.style.gridColumn = String(e.group);
      }
      b.addEventListener('click', function () { onSelect(z); });
      container.appendChild(b);
      cells[z] = b;
    }
    // placeholders pointing at the f-block, and a spacer row
    [[6, '57–71'], [7, '89–103']].forEach(function (p) {
      const ph = el('div', 'placeholder', p[1]);
      ph.style.gridRow = String(p[0]); ph.style.gridColumn = '3';
      container.appendChild(ph);
    });
    const gap = el('div', 'gap'); gap.style.gridRow = '8';
    container.appendChild(gap);
    return cells;
  };

  function fmtTemp(k) { return k === null || k === undefined ? '?' : Math.round(k) + ' K'; }

  UI.renderElementInfo = function (container, z, natural) {
    const e = E[z];
    container.innerHTML = '';
    const big = el('div', 'big');
    big.style.background = e.color;
    big.style.color = luminance(e.color) > 140 ? '#111' : '#fff';
    big.innerHTML = '<small>' + z + '</small>' + e.symbol;
    container.appendChild(big);
    container.appendChild(el('div', 'name', e.name));
    const facts1 = el('div', 'facts',
      e.categoryLabel + ' · ' + e.mass + ' u · EN ' + (e.en ? e.en.toFixed(2) : '—') + (e.enEstimated ? ' (est.)' : ''));
    container.appendChild(facts1);
    const val = e.valences.join(', ');
    const form = natural && CE.DIATOMIC.has(z) ? CE.chem.pretty(e.symbol + '2') : e.symbol + (e.isMetal ? ' (metal)' : '');
    container.appendChild(el('div', 'facts',
      'Valence ' + val + ' · mp ' + fmtTemp(e.mp) + ' · bp ' + fmtTemp(e.bp) + (e.radioactive ? ' · radioactive' : '')));
    const f3 = el('div', 'facts', 'Adds as ' + form + ' · covalent radius ' + e.covRadius + ' pm');
    f3.style.gridColumn = '1 / -1';
    container.appendChild(f3);
  };

  function luminance(hex) {
    const r = parseInt(hex.substr(1, 2), 16), g = parseInt(hex.substr(3, 2), 16), b = parseInt(hex.substr(5, 2), 16);
    return 0.299 * r + 0.587 * g + 0.114 * b;
  }

  function packedToCss(c) {
    return 'rgb(' + (c & 255) + ',' + ((c >> 8) & 255) + ',' + ((c >> 16) & 255) + ')';
  }

  UI.renderContents = function (container, species, onPick) {
    container.innerHTML = '';
    if (!species || !species.length) {
      container.appendChild(el('div', 'empty', 'The jar is empty. Pick an element and drop some in.'));
      return;
    }
    species.slice(0, 18).forEach(function (s) {
      const row = el('div', 'sp');
      const sw = el('span', 'sw'); sw.style.background = packedToCss(s.color);
      const nm = el('span', 'nm');
      nm.innerHTML = '<b></b> <small></small>';
      nm.firstChild.textContent = s.formula;
      nm.lastChild.textContent = s.name || '';
      const mol = el('span', 'mol', formatMol(s.mol));
      const ph = el('span', 'ph', s.phase === 'gas' ? 'gas' : 'liq/sol');
      row.title = (s.name ? s.name + ' · ' : '') + s.molecules + ' particles';
      row.append(sw, nm, mol, ph);
      row.addEventListener('click', function () { onPick(s); });
      container.appendChild(row);
    });
    if (species.length > 18) container.appendChild(el('div', 'empty', '+ ' + (species.length - 18) + ' more species'));
  };

  function formatMol(m) {
    if (m >= 10) return m.toFixed(1) + ' mol';
    if (m >= 0.1) return m.toFixed(2) + ' mol';
    return (m * 1000).toFixed(0) + ' mmol';
  }
  UI.formatMol = formatMol;

  UI.renderInspector = function (box, info, pinned) {
    box.innerHTML = '';
    const t = el('div', 'title', info.formula);
    box.appendChild(t);
    box.appendChild(el('div', 'sub', (info.name || (info.atomCount > 1 ? 'Unnamed molecule' : 'Free atom')) +
      ' · ' + info.atomCount + ' atom' + (info.atomCount > 1 ? 's' : '')));
    if (info.bonds.length) {
      const ul = el('ul');
      info.bonds.slice(0, 6).forEach(function (b) {
        ul.appendChild(el('li', null, b.bond + ' ×' + b.count + ' — ' + b.kind + ', ' + b.energy + ' kJ/mol, ' + b.length + ' pm'));
      });
      if (info.bonds.length > 6) ul.appendChild(el('li', null, '… ' + (info.bonds.length - 6) + ' more bond types'));
      box.appendChild(ul);
    }
    const a = info.atom;
    let line = 'Hovered atom: ' + a.name + ' (' + a.symbol + ')';
    if (info.geometry) line += ' · ' + info.geometry.name + ' ' + info.geometry.angle + '°' +
      (info.geometry.lonePairs ? ', ' + info.geometry.lonePairs + ' lone pair' + (info.geometry.lonePairs > 1 ? 's' : '') : '');
    box.appendChild(el('div', null, line));
    const q = info.charge;
    const qs = Math.abs(q) < 0.05 ? 'neutral' : (Math.abs(q) >= 0.95 ? (q > 0 ? '+' : '−') + Math.abs(q).toFixed(0) + ' (ion)' : 'δ' + (q > 0 ? '+' : '−') + Math.abs(q).toFixed(2));
    box.appendChild(el('div', null, 'Charge ' + qs + ' · free valence ' + info.freeValence + ' · local ' + Math.round(info.localT) + ' K'));
    if (info.geometry && info.geometry.angle !== 180) {
      box.appendChild(el('div', 'sub', 'Drawn flattened to 2D; the angle above is the real 3D VSEPR value.'));
    }
    if (pinned) box.appendChild(el('div', 'pinned', 'Pinned · click empty space to release'));
  };

  CE.UI = UI;
})(window.CE = window.CE || {});
