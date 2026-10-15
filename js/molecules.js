// Molecule library: built-in presets plus molecules the user saved from the jar (kept in this
// browser's local storage). Templates are { atoms: [[z, dx, dy], ...], bonds: [[a, b, order]] }
// in world units (1 unit = 10 pm); bond lengths and angles relax to their rest values once added.
(function (CE) {
  'use strict';

  const PRESETS = [
    { id: 'h2o', name: 'Water', formula: 'H₂O',
      atoms: [[8, 0, 0], [1, 9.7, 0], [1, -2.43, 9.39]], bonds: [[0, 1, 1], [0, 2, 1]] },
    { id: 'co2', name: 'Carbon dioxide', formula: 'CO₂', gas: true,
      atoms: [[6, 0, 0], [8, 12.4, 0], [8, -12.4, 0]], bonds: [[0, 1, 2], [0, 2, 2]] },
    { id: 'ch4', name: 'Methane', formula: 'CH₄', gas: true,
      atoms: [[6, 0, 0], [1, 10.7, 0], [1, -10.7, 0], [1, 0, 10.7], [1, 0, -10.7]],
      bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]] },
    { id: 'nh3', name: 'Ammonia', formula: 'NH₃', gas: true,
      atoms: [[7, 0, 0], [1, 10.2, 0], [1, -5.1, 8.83], [1, -5.1, -8.83]],
      bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1]] },
    { id: 'co', name: 'Carbon monoxide', formula: 'CO', gas: true,
      atoms: [[6, 0, 0], [8, 11.2, 0]], bonds: [[0, 1, 3]] },
    { id: 'hcl', name: 'Hydrogen chloride', formula: 'HCl', gas: true,
      atoms: [[1, 0, 0], [17, 13.3, 0]], bonds: [[0, 1, 1]] },
    { id: 'nacl', name: 'Sodium chloride (salt)', formula: 'NaCl',
      atoms: [[11, 0, 0], [17, 29.6, 0]], bonds: [[0, 1, 1]] },
    { id: 'naoh', name: 'Sodium hydroxide', formula: 'NaOH',
      atoms: [[11, 0, 0], [8, 25.5, 0], [1, 35.2, 0]], bonds: [[0, 1, 1], [1, 2, 1]] },
    { id: 'h2o2', name: 'Hydrogen peroxide', formula: 'H₂O₂',
      atoms: [[1, -3, -9.2], [8, 0, 0], [8, 14.2, 0], [1, 17.2, 9.2]],
      bonds: [[0, 1, 1], [1, 2, 1], [2, 3, 1]] },
    { id: 'etoh', name: 'Ethanol', formula: 'C₂H₅OH',
      atoms: [[6, 0, 0], [6, 15.2, 0], [8, 15.2, -14.2], [1, 24.9, -14.2],
        [1, 0, -10.7], [1, 0, 10.7], [1, -10.7, 0], [1, 15.2, 10.7], [1, 25.9, 0]],
      bonds: [[0, 1, 1], [1, 2, 1], [2, 3, 1], [0, 4, 1], [0, 5, 1], [0, 6, 1], [1, 7, 1], [1, 8, 1]] },
    { id: 'c6h6', name: 'Benzene', formula: 'C₆H₆', benzene: true },
    { id: 'h2so4', name: 'Sulfuric acid', formula: 'H₂SO₄',
      atoms: [[16, 0, 0], [8, 0, -14.9], [8, 0, 14.9], [8, -17.1, 0], [1, -26.8, 0], [8, 17.1, 0], [1, 26.8, 0]],
      bonds: [[0, 1, 2], [0, 2, 2], [0, 3, 1], [3, 4, 1], [0, 5, 1], [5, 6, 1]] }
  ];

  // benzene: a hexagon of carbons with alternating double bonds, one H on each
  (function () {
    const b = PRESETS.find(function (p) { return p.benzene; });
    b.atoms = []; b.bonds = [];
    for (let k = 0; k < 6; k++) {
      const a = k * Math.PI / 3;
      b.atoms.push([6, 14 * Math.cos(a), 14 * Math.sin(a)]);
    }
    for (let k = 0; k < 6; k++) {
      const a = k * Math.PI / 3;
      b.atoms.push([1, 24.7 * Math.cos(a), 24.7 * Math.sin(a)]);
      b.bonds.push([k, (k + 1) % 6, k % 2 ? 1 : 2]);
      b.bonds.push([k, k + 6, 1]);
    }
  })();

  const KEY = 'chemistry-engine.savedMolecules';

  function loadSaved() {
    try {
      const raw = window.localStorage.getItem(KEY);
      const list = raw ? JSON.parse(raw) : [];
      return Array.isArray(list) ? list.filter(function (m) { return m && m.atoms && m.bonds; }) : [];
    } catch (e) { return []; }
  }

  function storeSaved(list) {
    try { window.localStorage.setItem(KEY, JSON.stringify(list)); return true; } catch (e) { return false; }
  }

  CE.molecules = { PRESETS: PRESETS, loadSaved: loadSaved, storeSaved: storeSaved };
})(window.CE = window.CE || {});
