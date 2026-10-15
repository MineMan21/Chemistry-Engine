// How substances look to the eye, for the real-world view: colour, transparency, whether a
// condensed lump is a liquid or a solid, and how strongly the vapour tints the air.
(function (CE) {
  'use strict';

  const E = CE.elements;

  // Bulk colour of the pure elements as people see them (most metals are silvery grey).
  const ELEMENT_LOOK = {
    H: '#dfe9f5', He: '#e8f1f8', C: '#1f1f22', N: '#cfe0f0', O: '#9fc3ff', F: '#e8eb9a',
    Ne: '#e8f1f8', Na: '#c9ccd1', Mg: '#b9bec4', Al: '#c6cacf', Si: '#5f6670', P: '#efe6c8',
    S: '#e6d23a', Cl: '#c8dc4a', Ar: '#e8f1f8', K: '#c4c7cc', Ca: '#d2d4d6', Ti: '#a9adb2',
    V: '#a4a9b0', Cr: '#b8bec6', Mn: '#a6a19c', Fe: '#7f8388', Co: '#8e939b', Ni: '#a8a59d',
    Cu: '#b8733d', Zn: '#a9b1b8', Ga: '#c3c7cc', Ge: '#8a8f94', As: '#6e7276', Se: '#4a4a4f',
    Br: '#7a1a0a', Kr: '#e8f1f8', Rb: '#c0c2c6', Sr: '#d0d0cc', Ag: '#cfd2d6', Cd: '#a9adb3',
    In: '#bfc2c6', Sn: '#b8bcc0', Sb: '#8e9196', Te: '#9a9da1', I: '#3a1e44', Xe: '#e8f1f8',
    Cs: '#d8c58a', Ba: '#bfc0bd', W: '#8d9196', Pt: '#c9cbcf', Au: '#e4b84a', Hg: '#c7cad1',
    Pb: '#6f757e', Bi: '#c9b9c9', U: '#7d8186', Pu: '#8a8580', B: '#4a3f36', Os: '#9aa6bd',
    Ir: '#c0c4ca', Pd: '#bfc2c5', Rh: '#c6c9cc', Ru: '#9fa3a8', Mo: '#9ea3a8', Nb: '#9aa0a8',
    Ta: '#8f949a', Re: '#a6aab0', Zr: '#b5b9bd', Hf: '#a9adb2', Li: '#c7c9cc', Be: '#9fa4a8'
  };

  // Compounds by Hill key: [colour, liquid opacity, gas tint strength]
  const COMPOUND_LOOK = {
    H2O: ['#6fb4ff', 0.42, 0], H2O2: ['#a8d4ff', 0.4, 0], ClNa: ['#f4f5f6', 1, 0],
    '(ClNa)n': ['#f4f5f6', 1, 0], ClK: ['#f2f2f2', 1, 0], CaCl2: ['#f0f0f0', 1, 0],
    O2S: ['#f2f2f2', 0.3, 0], NO2: ['#b5541c', 0.6, 0.55], N2O4: ['#e8d8b0', 0.4, 0.05],
    O3: ['#7aa6ff', 0.6, 0.15], ClO2: ['#e0a020', 0.6, 0.4], Fe2O3: ['#8b3a1a', 1, 0],
    FeO: ['#2b2a28', 1, 0], Fe3O4: ['#2a2a2c', 1, 0], CuO: ['#1d1d1f', 1, 0],
    Cu2O: ['#a8321e', 1, 0], Cl2Cu: ['#2f8f8a', 1, 0], Cl3Fe: ['#5a3a16', 1, 0],
    Cl2Fe: ['#a6c48e', 1, 0], MgO: ['#f3f3f1', 1, 0], CaO: ['#f1f1ee', 1, 0],
    Al2O3: ['#f2f2f2', 1, 0], O2Si: ['#e8eef3', 1, 0], OZn: ['#f6f6f2', 1, 0],
    O2Ti: ['#fbfbfb', 1, 0], Cr2O3: ['#3d6e3a', 1, 0], NiO: ['#5a7a3a', 1, 0],
    HgO: ['#d0571f', 1, 0], OPb: ['#d9a400', 1, 0], AgCl: ['#f0f0ec', 1, 0],
    Ag2O: ['#2a2a2a', 1, 0], FeS: ['#3a3530', 1, 0], CuS: ['#1f1f24', 1, 0],
    SZn: ['#f0f0e8', 1, 0], HNaO: ['#f4f4f4', 1, 0], Na2O: ['#f0f0f0', 1, 0],
    O10P4: ['#f6f6f6', 1, 0], B2O3: ['#eef2f4', 1, 0], CS2: ['#f0f0d8', 0.3, 0],
    Br2: ['#7a1a0a', 0.9, 0.6], I2: ['#3a1e44', 1, 0.75], Cl2: ['#c8dc4a', 0.7, 0.28],
    F2: ['#e8eb9a', 0.6, 0.08], S8: ['#e6d23a', 1, 0]
  };

  // boiling points (K) of common compounds, so their vapour is drawn as gas above them
  const COMPOUND_BP = { H2O: 373, H2O2: 423, Br2: 332, I2: 457, Cl2: 239, F2: 85, O3: 161,
    O2S: 263, NO2: 294, N2O4: 294, ClO2: 284, CS2: 319, CH4: 112, CO2: 195, CO: 82, H3N: 240,
    ClH: 188, FH: 293, H2S: 213, CH4O: 338, C2H6O: 351, C6H6: 353, C8H18: 399 };

  const TRANSITION_TINT = { 24: '#3d6e3a', 25: '#6a4a3a', 26: '#8b3a1a', 27: '#3f5f9a', 28: '#5a7a3a', 29: '#2f6f9a' };

  function hex(c) { return [parseInt(c.substr(1, 2), 16), parseInt(c.substr(3, 2), 16), parseInt(c.substr(5, 2), 16)]; }

  // key: species key from Simulation.analyze; counts: {z: n}; atoms: atoms per molecule.
  // Returns { rgb: [r,g,b], opacity (when condensed), solid (bool at temperature T), gas tint }.
  function lookOf(key, counts, atoms, ionic, T) {
    const zs = Object.keys(counts).map(Number);
    const network = key.slice(-2) === ')n';
    if (COMPOUND_LOOK[key]) {
      const c = COMPOUND_LOOK[key];
      let solid = c[1] >= 1;
      let rgb = hex(c[0]), opacity = c[1];
      if (key === 'H2O' && T < 273) { solid = true; rgb = hex('#dcefff'); opacity = 0.85; } // ice
      if ((key === 'Br2' || key === 'Cl2' || key === 'F2' || key === 'O3') && T < 200) solid = false;
      if (key === 'I2') solid = T < 387;
      return { rgb: rgb, opacity: opacity, solid: solid, gas: c[2], bp: COMPOUND_BP[key] || null };
    }
    if (zs.length === 1) {
      const e = E[zs[0]];
      const rgb = hex(ELEMENT_LOOK[e.symbol] || (e.isMetal ? '#b4b8bd' : e.color));
      const solid = network || e.isMetal ? (e.mp === null || T < e.mp) : (e.mp !== null && T < e.mp);
      const liquidOpacity = e.isMetal ? 1 : e.symbol === 'N' || e.symbol === 'O' || e.symbol === 'H' ? 0.3 : 0.6;
      return { rgb: rgb, opacity: solid ? 1 : liquidOpacity, solid: solid, gas: 0, bp: network ? null : e.bp };
    }
    if (network && !ionic) {
      // covalent networks: the count-weighted colour of their elements
      let r = 0, g = 0, b = 0, tot = 0;
      zs.forEach(function (z) {
        const c = hex(ELEMENT_LOOK[E[z].symbol] || E[z].color), w = counts[z];
        r += c[0] * w; g += c[1] * w; b += c[2] * w; tot += w;
      });
      return { rgb: [r / tot, g / tot, b / tot], opacity: 1, solid: true, gas: 0, bp: null };
    }
    if (ionic || network) {
      // salts and oxides: white, unless a coloured transition metal is in them
      let rgb = hex('#ecedeb');
      for (let k = 0; k < zs.length; k++) if (TRANSITION_TINT[zs[k]]) rgb = hex(TRANSITION_TINT[zs[k]]);
      return { rgb: rgb, opacity: 1, solid: true, gas: 0, bp: null };
    }
    // small covalent molecules: clear liquids
    return { rgb: hex('#d6e6f2'), opacity: 0.32, solid: false, gas: 0, bp: COMPOUND_BP[key] || null };
  }

  CE.look = { lookOf: lookOf };
})(typeof window !== 'undefined' ? (window.CE = window.CE || {}) : (globalThis.CE = globalThis.CE || {}));
