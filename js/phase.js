// Phase data for substances with documented triple and critical points, and the phase diagrams
// built from them. Values are standard reference data (NIST Chemistry WebBook / CRC Handbook);
// pressures in atm, temperatures in K. Entries marked est have an estimated value.
(function (CE) {
  'use strict';

  // key: species key (Hill formula, element symbol for monatomic elements, X2 for diatomics)
  // Tt, Pt: triple point; Tc, Pc: critical point; Tm, Tb: normal melting / boiling point
  // Ts: sublimation point at 1 atm (when the triple point lies above 1 atm)
  // hf: enthalpy of fusion (kJ/mol); slope: melting line dT/dP (K/atm), default +0.015
  // meltTop: pressure where this simple melting line stops being valid (other solid phases)
  // est: names of estimated values
  const RAW = {
    H2: { name: 'Hydrogen', Tt: 13.80, Pt: 0.0695, Tc: 33.15, Pc: 12.96, Tm: 13.99, Tb: 20.27, hf: 0.117 },
    He: { name: 'Helium', Tc: 5.195, Pc: 2.245, Tb: 4.222, helium: true,
      note: 'Helium has no triple point with its solid: it stays liquid down to absolute zero unless squeezed above about 25 atm.' },
    N2: { name: 'Nitrogen', Tt: 63.151, Pt: 0.1237, Tc: 126.19, Pc: 33.50, Tm: 63.15, Tb: 77.355, hf: 0.72 },
    O2: { name: 'Oxygen', Tt: 54.36, Pt: 0.00144, Tc: 154.58, Pc: 49.77, Tm: 54.36, Tb: 90.19, hf: 0.444 },
    F2: { name: 'Fluorine', Tt: 53.48, Pt: 0.00249, Tc: 144.41, Pc: 51.04, Tm: 53.48, Tb: 85.03, hf: 0.51 },
    Ne: { name: 'Neon', Tt: 24.556, Pt: 0.428, Tc: 44.49, Pc: 26.44, Tm: 24.56, Tb: 27.10, hf: 0.335 },
    Cl2: { name: 'Chlorine', Tt: 172.22, Pt: 0.01374, Tc: 416.9, Pc: 78.87, Tm: 171.6, Tb: 239.11, hf: 6.40 },
    Ar: { name: 'Argon', Tt: 83.806, Pt: 0.680, Tc: 150.69, Pc: 47.99, Tm: 83.81, Tb: 87.30, hf: 1.18 },
    Br2: { name: 'Bromine', Tt: 265.90, Pt: 0.0572, Tc: 588, Pc: 102, Tm: 265.8, Tb: 332.0, hf: 10.57 },
    Kr: { name: 'Krypton', Tt: 115.78, Pt: 0.7224, Tc: 209.48, Pc: 54.53, Tm: 115.79, Tb: 119.93, hf: 1.64 },
    I2: { name: 'Iodine', Tt: 386.65, Pt: 0.1191, Tc: 819, Pc: 115.5, Tm: 386.85, Tb: 457.4, hf: 15.52 },
    Xe: { name: 'Xenon', Tt: 161.40, Pt: 0.8055, Tc: 289.73, Pc: 57.65, Tm: 161.4, Tb: 165.05, hf: 2.27 },
    Rn: { name: 'Radon', Tt: 202, Pt: 0.69, Tc: 377, Pc: 62, Tm: 202, Tb: 211.5, hf: 3.25, est: ['Pt'] },
    Hg: { name: 'Mercury', Tt: 234.32, Pt: 1.63e-9, Tc: 1750, Pc: 1698, Tm: 234.32, Tb: 629.88, hf: 2.29 },
    Li: { name: 'Lithium', Tm: 453.65, Tb: 1603, Tc: 3223, Pc: 661, hf: 3.0, est: ['Pt', 'Tc', 'Pc'] },
    Na: { name: 'Sodium', Tm: 370.94, Tb: 1156.1, Tc: 2573, Pc: 345, hf: 2.60, est: ['Pt', 'Tc', 'Pc'] },
    K: { name: 'Potassium', Tm: 336.7, Tb: 1032, Tc: 2223, Pc: 158, hf: 2.33, est: ['Pt', 'Tc', 'Pc'] },
    Rb: { name: 'Rubidium', Tm: 312.45, Tb: 961, Tc: 2093, Pc: 158, hf: 2.19, est: ['Pt', 'Tc', 'Pc'] },
    Cs: { name: 'Caesium', Tm: 301.7, Tb: 944, Tc: 1938, Pc: 92.8, hf: 2.09, est: ['Pt', 'Tc', 'Pc'] },
    S: { name: 'Sulfur', Tm: 388.36, Tb: 717.8, Tc: 1314, Pc: 204, hf: 1.72, est: ['Pt', 'Tc', 'Pc'],
      note: 'Sulfur also has two solid forms (rhombic and monoclinic); only one solid region is shown.' },
    As: { name: 'Arsenic', Tt: 1090, Pt: 35.8, Tc: 1673, Pc: 220, Ts: 887, hf: 24.4, est: ['Tc', 'Pc'],
      note: 'Grey arsenic sublimes at 1 atm: it only melts under pressure.' },
    H2O: { name: 'Water', Tt: 273.16, Pt: 0.006037, Tc: 647.096, Pc: 217.75, Tm: 273.15, Tb: 373.12, hf: 6.01,
      slope: -0.0074, meltTop: 2070,
      note: 'Ice melts at a lower temperature when squeezed (the melting line leans left). Above about 2000 atm other forms of ice take over; they are not drawn.' },
    CO2: { name: 'Carbon dioxide', Tt: 216.58, Pt: 5.117, Tc: 304.13, Pc: 72.8, Ts: 194.69, hf: 9.02, slope: 0.014,
      note: 'At 1 atm dry ice turns straight into gas (sublimes) at 195 K: liquid CO₂ needs more than 5.1 atm.' },
    H3N: { name: 'Ammonia', Tt: 195.40, Pt: 0.0598, Tc: 405.4, Pc: 111.3, Tm: 195.42, Tb: 239.82, hf: 5.66 },
    CH4: { name: 'Methane', Tt: 90.69, Pt: 0.1155, Tc: 190.56, Pc: 45.39, Tm: 90.7, Tb: 111.65, hf: 0.94 },
    CO: { name: 'Carbon monoxide', Tt: 68.15, Pt: 0.1515, Tc: 132.86, Pc: 34.48, Tm: 68.13, Tb: 81.6, hf: 0.83 },
    ClH: { name: 'Hydrogen chloride', Tt: 158.97, Pt: 0.1366, Tc: 324.7, Pc: 82.0, Tm: 158.97, Tb: 188.1, hf: 2.0 },
    FH: { name: 'Hydrogen fluoride', Tm: 189.6, Tb: 292.7, Tc: 461, Pc: 64.0, hf: 4.58, est: ['Pt'] },
    H2S: { name: 'Hydrogen sulfide', Tt: 187.66, Pt: 0.2296, Tc: 373.1, Pc: 88.2, Tm: 187.7, Tb: 212.8, hf: 2.38 },
    O2S: { name: 'Sulfur dioxide', Tt: 197.69, Pt: 0.01648, Tc: 430.64, Pc: 77.8, Tm: 197.7, Tb: 263.1, hf: 7.40 },
    N2O: { name: 'Nitrous oxide', Tt: 182.34, Pt: 0.866, Tc: 309.52, Pc: 71.5, Tm: 182.3, Tb: 184.67, hf: 6.54 },
    NO: { name: 'Nitric oxide', Tt: 109.5, Pt: 0.216, Tc: 180, Pc: 64, Tm: 109.5, Tb: 121.4, hf: 2.3 },
    O3: { name: 'Ozone', Tm: 80.7, Tb: 161.1, Tc: 261.1, Pc: 54.6, hf: 2.1, est: ['Pt'] },
    F6S: { name: 'Sulfur hexafluoride', Tt: 223.5, Pt: 2.23, Tc: 318.7, Pc: 37.1, Ts: 209.4, hf: 5.0 },
    C6H6: { name: 'Benzene', Tt: 278.68, Pt: 0.0472, Tc: 562.0, Pc: 48.3, Tm: 278.68, Tb: 353.2, hf: 9.87 },
    C2H6O: { name: 'Ethanol', Tt: 150, Pt: 4.3e-9, Tc: 514.0, Pc: 60.6, Tm: 159.0, Tb: 351.4, hf: 4.9 },
    CH4O: { name: 'Methanol', Tt: 175.6, Pt: 1.84e-6, Tc: 512.6, Pc: 79.9, Tm: 175.6, Tb: 337.8, hf: 3.22 },
    C2H6: { name: 'Ethane', Tt: 90.35, Pt: 1.1e-5, Tc: 305.32, Pc: 48.2, Tm: 90.35, Tb: 184.6, hf: 2.86 },
    C3H8: { name: 'Propane', Tt: 85.5, Pt: 1.7e-9, Tc: 369.8, Pc: 41.9, Tm: 85.5, Tb: 231.1, hf: 3.52 },
    C2H4: { name: 'Ethylene', Tt: 104.0, Pt: 0.00120, Tc: 282.34, Pc: 49.7, Tm: 104.0, Tb: 169.4, hf: 3.35 },
    C2H2: { name: 'Acetylene', Tt: 192.4, Pt: 1.26, Tc: 308.3, Pc: 60.6, Ts: 189.3, hf: 3.8 },
    CCl4: { name: 'Carbon tetrachloride', Tm: 250.3, Tb: 349.9, Tc: 556.6, Pc: 44.9, hf: 2.5, est: ['Pt'] },
    CS2: { name: 'Carbon disulfide', Tm: 161.6, Tb: 319.0, Tc: 552, Pc: 78, hf: 4.39, est: ['Pt'] },
    H2O2: { name: 'Hydrogen peroxide', Tm: 272.7, Tb: 423.4, Tc: 728, Pc: 214, hf: 12.5, est: ['Pt', 'Tc', 'Pc'] }
  };

  const R = 0.008314; // kJ/(mol K)

  function vapAnchors(d) {
    const a = [[d.Tt, d.Pt]];
    if (d.Tb && d.Tb > d.Tt && d.Tb < d.Tc && d.Pt < 1) a.push([d.Tb, 1]);
    a.push([d.Tc, d.Pc]);
    return a;
  }
  // ln P linear in 1/T between anchors (Clausius-Clapeyron, piecewise)
  function lnInterp(a, T) {
    let k = 0;
    while (k < a.length - 2 && T > a[k + 1][0]) k++;
    const [T0, P0] = a[k], [T1, P1] = a[k + 1];
    const f = (1 / T - 1 / T0) / (1 / T1 - 1 / T0);
    return Math.exp(Math.log(P0) + f * (Math.log(P1) - Math.log(P0)));
  }

  // Completes an entry: estimated triple point, enthalpies for the sublimation line.
  function complete(key, raw) {
    const d = Object.assign({ key: key, slope: 0.015, est: [] }, raw);
    if (d.helium) return d;
    if (d.Tt === undefined) {
      // triple point ~ melting point; its pressure from the vapour curve through (Tb, 1 atm)-(Tc, Pc)
      d.Tt = d.Tm;
      d.Pt = lnInterp([[d.Tb, 1], [d.Tc, d.Pc]], d.Tt);
      if (d.est.indexOf('Tt') < 0) d.est = d.est.concat(['Tt']);
    }
    if (d.Tm === undefined) d.Tm = d.Ts ? null : d.Tt;
    const a = vapAnchors(d);
    // enthalpy of vaporisation near the triple point from the first vapour segment's slope
    const hv = -R * (Math.log(a[1][1]) - Math.log(a[0][1])) / (1 / a[1][0] - 1 / a[0][0]);
    d.hv = hv;
    d.hs = hv + (d.hf || 1.5 * R * d.Tt);
    return d;
  }

  const DATA = {};
  Object.keys(RAW).forEach(function (k) { DATA[k] = complete(k, RAW[k]); });

  // Pressure of the gas in equilibrium with the liquid (T >= Tt) or the solid (T < Tt).
  function vaporP(d, T) { return lnInterp(vapAnchors(d), T); }
  function sublimationP(d, T) {
    if (d.Ts) {
      // straight line in (1/T, ln P) through the sublimation point (Ts, 1 atm) and the triple point
      const b = Math.log(d.Pt) / (1 / d.Tt - 1 / d.Ts);
      return Math.exp(Math.log(d.Pt) + b * (1 / T - 1 / d.Tt));
    }
    return Math.exp(Math.log(d.Pt) - d.hs / R * (1 / T - 1 / d.Tt));
  }
  function meltT(d, P) {
    if (d.helium) return 0.4 * Math.sqrt(Math.max(0, P - 25)); // rough helium melting curve
    return d.Tt + d.slope * (P - d.Pt);
  }

  // 'solid' | 'liquid' | 'gas' | 'supercritical'
  function phaseAt(d, T, P) {
    if (d.helium) {
      if (T > d.Tc) return P > d.Pc ? 'supercritical' : 'gas';
      if (P > 25 && T < meltT(d, P)) return 'solid';
      return P >= heliumVap(d, T) ? 'liquid' : 'gas';
    }
    if (P >= d.Pt) {
      if (T < meltT(d, P)) return 'solid';
      if (T >= d.Tc) return P >= d.Pc ? 'supercritical' : 'gas';
      return P >= vaporP(d, T) ? 'liquid' : 'gas';
    }
    if (T < d.Tt && P >= sublimationP(d, T)) return 'solid';
    return 'gas';
  }
  function heliumVap(d, T) { return lnInterp([[d.Tb, 1], [d.Tc, d.Pc]], Math.max(T, 0.5)); }

  // Plot ranges and boundary polylines [[T, P], ...] for drawing.
  function diagram(d) {
    const lines = [];
    const lo = d.helium ? 1e-3 : Math.min(d.Pt, 1e-3);
    const Pmin = Math.pow(10, Math.floor(Math.log10(lo)) - 1);
    const Pmax = Math.pow(10, Math.ceil(Math.log10(Math.max(d.Pc, 1) * 20)));
    const Tmax = Math.ceil(d.Tc * 1.25 / 10) * 10;
    const N = 80;
    if (d.helium) {
      const vap = [];
      for (let k = 0; k <= N; k++) { const T = 0.8 + (d.Tc - 0.8) * k / N; vap.push([T, heliumVap(d, T)]); }
      lines.push({ kind: 'vaporisation', pts: vap });
      const mel = [];
      for (let k = 0; k <= N; k++) { const P = 25 * Math.pow(Pmax / 25, k / N); mel.push([meltT(d, P), P]); }
      lines.push({ kind: 'melting', pts: mel });
    } else {
      const vap = [];
      for (let k = 0; k <= N; k++) { const T = d.Tt + (d.Tc - d.Tt) * k / N; vap.push([T, vaporP(d, T)]); }
      lines.push({ kind: 'vaporisation', pts: vap });
      const sub = [];
      for (let k = 0; k <= N; k++) {
        const T = d.Tt * (1 - k / N * 0.9);
        const P = sublimationP(d, T);
        sub.push([T, P]);
        if (P < Pmin) break;
      }
      lines.push({ kind: 'sublimation', pts: sub });
      const mel = [], top = Math.min(Pmax, d.meltTop || Pmax);
      for (let k = 0; k <= N; k++) { const P = d.Pt * Math.pow(top / d.Pt, k / N); mel.push([meltT(d, P), P]); }
      lines.push({ kind: 'melting', pts: mel });
    }
    return { Tmin: 0, Tmax: Tmax, Pmin: Pmin, Pmax: Pmax, lines: lines };
  }

  // Data for a species key (or null when no documented diagram exists).
  function forKey(key) { return DATA[key] || null; }
  function forElement(z) {
    const e = CE.elements[z];
    return DATA[e.symbol + '2'] && CE.DIATOMIC.has(z) ? DATA[e.symbol + '2'] : DATA[e.symbol] || null;
  }
  // Normal melting point of a species, for telling frozen molecular solids apart from liquids.
  function meltingPoint(key) {
    const d = DATA[key];
    if (!d) return 0;
    return d.helium ? 0 : d.Tm || d.Ts || d.Tt;
  }

  CE.phase = { DATA: DATA, forKey: forKey, forElement: forElement, phaseAt: phaseAt, diagram: diagram,
    vaporP: vaporP, sublimationP: sublimationP, meltT: meltT, meltingPoint: meltingPoint };
})(typeof window !== 'undefined' ? (window.CE = window.CE || {}) : (globalThis.CE = globalThis.CE || {}));
