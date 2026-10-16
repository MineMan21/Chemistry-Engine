// Rule-based chemistry: bond classification, bond energies, bond lengths, valence rules,
// VSEPR geometry and molecule naming. Everything here is pure data/functions (no DOM).
(function (CE) {
  'use strict';

  const E = CE.elements;
  const NZ = 119; // index by atomic number

  // ---------------------------------------------------------------------------------------
  // Bond enthalpies (kJ/mol). '-' single, '=' double, '#' triple. Average literature values.
  // ---------------------------------------------------------------------------------------
  const KNOWN_BONDS = `
H-H 436,H-C 411,H-N 386,H-O 459,H-F 565,H-Cl 428,H-Br 362,H-I 295,H-S 363,H-P 322,H-Si 318,
H-B 389,H-Li 243,H-Na 186,H-K 174,H-Se 276,H-Ge 288,H-As 247,H-Te 238,H-Sn 251,
C-C 346,C=C 602,C#C 835,C-N 305,C=N 615,C#N 887,C-O 358,C=O 799,C#O 1072,C-F 485,C-Cl 327,
C-Br 285,C-I 213,C-S 272,C=S 573,C-Si 318,C-P 264,C-B 356,
N-N 167,N=N 418,N#N 942,N-O 201,N=O 607,N-F 283,N-Cl 313,N-S 280,
O-O 142,O=O 494,O-F 190,O-Cl 218,O-Br 201,O-I 201,O-Si 452,O-P 335,O=P 544,O-S 265,O=S 522,
O-B 536,O-Se 343,O=Se 425,O-Xe 84,O=Xe 170,
F-F 155,Cl-Cl 240,Br-Br 190,I-I 148,F-Cl 253,F-Br 250,F-I 271,Cl-Br 218,Cl-I 208,Br-I 175,
S-S 226,S=S 425,P-P 201,P#P 489,Si-Si 222,B-B 293,Se-Se 172,Ge-Ge 188,As-As 146,Sn-Sn 146,
Te-Te 126,Sb-Sb 121,Si-F 565,Si-Cl 381,S-F 284,S-Cl 255,P-F 490,P-Cl 326,P-Br 264,P-S 230,
P=S 335,B-F 613,B-Cl 456,Xe-F 130,Kr-F 50,Rn-F 150,
Li-Li 105,Na-Na 74,K-K 57,Rb-Rb 49,Cs-Cs 44,Al-Cl 420,Al-O 512
`;

  // Homonuclear single-bond estimates for elements not covered above, by category.
  const SELF_BOND_BY_CATEGORY = {
    alkali: 60, alkaline: 60, transition: 120, 'post-transition': 100, lanthanide: 100,
    actinide: 100, metalloid: 180, nonmetal: 150, halogen: 150, noble: 10
  };

  const known = {};
  KNOWN_BONDS.replace(/\s+/g, '').split(',').forEach(function (entry) {
    if (!entry) return;
    const m = entry.match(/^([A-Z][a-z]?)([-=#])([A-Z][a-z]?)(\d+)$/);
    if (!m) throw new Error('bad bond entry ' + entry);
    const a = CE.bySymbol[m[1]].z, b = CE.bySymbol[m[3]].z;
    const order = m[2] === '-' ? 1 : m[2] === '=' ? 2 : 3;
    known[key(a, b, order)] = +m[4];
  });
  function key(a, b, order) { return (a < b ? a + ':' + b : b + ':' + a) + ':' + order; }

  function selfBond(z) {
    const k = known[key(z, z, 1)];
    if (k) return k;
    return SELF_BOND_BY_CATEGORY[E[z].category];
  }

  // Bond class: 0 = none, 1 = covalent, 2 = ionic.
  function bondClass(a, b) {
    const ea = E[a], eb = E[b];
    if (ea.maxValence === 0 || eb.maxValence === 0) return 0; // noble gases (He, Ne, Ar...)
    if (ea.isMetal && eb.isMetal) return 0; // metallic bonding handled as cohesion
    const d = Math.abs(ea.en - eb.en);
    // a metal with a much more electronegative partner forms an ionic compound: textbook oxides,
    // halides and sulfides of the metals (MgO, Fe2O3, CuO, NaCl). Pauling's 50%-ionic gap is 1.7;
    // 1.5 also catches the transition-metal oxides (Fe-O 1.61, Cu-O 1.54), which are ionic solids
    if (d >= 1.5 && (ea.isMetal || eb.isMetal)) return 2;
    return 1;
  }

  // Elements able to form pi bonds (double/triple) with each other.
  const PI_STRONG = new Set([6, 7, 8]); // C N O
  const PI_WEAK = new Set([5, 15, 16, 34, 54]); // B P S Se Xe (with O)

  function maxOrder(a, b) {
    const cls = bondClass(a, b);
    if (cls === 0) return 0;
    if (cls === 2) {
      // ionic: one bond unit per ion-ion contact, so a doubly or triply charged ion bonds to
      // several counter-ions and the compound grows as an extended lattice with the textbook
      // ratio ((MgO)n, (Al2O3)n, (Fe2O3)n) instead of closing into small Mg=O molecules
      return 1;
    }
    if (a === 1 || b === 1) return 1;
    if (E[a].isMetal || E[b].isMetal) return Math.min(2, Math.min(E[a].maxValence, E[b].maxValence));
    const sa = PI_STRONG.has(a), sb = PI_STRONG.has(b);
    if (sa && sb) {
      if (a === 8 && b === 8) return 2;
      if (a === 8 || b === 8) return (a === 6 || b === 6) ? 3 : 2;
      return 3;
    }
    if ((sa || PI_WEAK.has(a)) && (sb || PI_WEAK.has(b))) {
      // heavier elements form pi bonds with C/N/O (S=O, P=O, C=S) but not with each other:
      // sulfur makes S8 rings and phosphorus P4 / networks of single bonds
      return (sa || sb) ? 2 : 1;
    }
    return 1;
  }

  // Bond enthalpy for a given order (total energy of the bond, not the increment).
  function bondEnergyRaw(a, b, order) {
    const k = known[key(a, b, order)];
    if (k) return k;
    const cls = bondClass(a, b);
    if (cls === 0) return 0;
    let e1 = known[key(a, b, 1)];
    if (!e1) {
      const d = Math.abs(E[a].en - E[b].en);
      // Pauling: D(AB) = average of D(AA), D(BB) + 96.5 * dEN^2
      e1 = (selfBond(a) + selfBond(b)) / 2 + 96.5 * d * d;
      e1 = Math.min(e1, cls === 2 ? 500 : 600);
    }
    if (order === 1) return e1;
    if (cls === 2) return e1 * (1 + 0.75 * (order - 1));
    return e1 * (order === 2 ? 1.75 : 2.45);
  }

  // Bond length (world units; 1 unit = 10 pm) from covalent radii, shortened for multiple bonds.
  function bondLength(a, b, order, cls) {
    let l = (E[a].covRadius + E[b].covRadius) / 10;
    if (cls === 2) l *= 1.1; // ionic distances slightly longer than covalent sums
    if (order === 2) l *= 0.87;
    else if (order === 3) l *= 0.79;
    return Math.max(l, 5);
  }

  // ---------------------------------------------------------------------------------------
  // Pair tables (flattened [a * NZ + b]) used by the simulation's inner loops.
  // ---------------------------------------------------------------------------------------
  const pairClass = new Uint8Array(NZ * NZ);
  const pairMaxOrder = new Uint8Array(NZ * NZ);
  const pairEnergy = new Float32Array(NZ * NZ * 4); // [pair*4 + order]
  const pairLength = new Float32Array(NZ * NZ * 4);
  for (let a = 1; a < NZ; a++) {
    for (let b = 1; b < NZ; b++) {
      const p = a * NZ + b;
      const cls = bondClass(a, b);
      pairClass[p] = cls;
      const mo = maxOrder(a, b);
      pairMaxOrder[p] = mo;
      for (let o = 1; o <= 3; o++) {
        pairEnergy[p * 4 + o] = o <= Math.max(mo, 1) ? bondEnergyRaw(a, b, o) : 0;
        pairLength[p * 4 + o] = bondLength(a, b, o, cls);
      }
    }
  }

  function bondEnergy(a, b, order) { return pairEnergy[(a * NZ + b) * 4 + order]; }
  // Energy needed to remove the top unit of a bond of the given order.
  function stepEnergy(a, b, order) {
    const p = (a * NZ + b) * 4;
    return pairEnergy[p + order] - (order > 1 ? pairEnergy[p + order - 1] : 0);
  }

  // ---------------------------------------------------------------------------------------
  // VSEPR
  // ---------------------------------------------------------------------------------------
  const VSEPR = {
    '1,0': ['linear', 180], '1,1': ['linear', 180], '1,2': ['linear', 180], '1,3': ['linear', 180],
    '2,0': ['linear', 180], '2,1': ['bent', 117], '2,2': ['bent', 104.5], '2,3': ['linear', 180],
    '3,0': ['trigonal planar', 120], '3,1': ['trigonal pyramidal', 107], '3,2': ['T-shaped', 90],
    '4,0': ['tetrahedral', 109.5], '4,1': ['seesaw', 102], '4,2': ['square planar', 90],
    '5,0': ['trigonal bipyramidal', 90], '5,1': ['square pyramidal', 90],
    '6,0': ['octahedral', 90], '6,1': ['pentagonal pyramidal', 72],
    '7,0': ['pentagonal bipyramidal', 72], '8,0': ['square antiprismatic', 70]
  };

  // Lone pairs on a main-group atom given its total covalent bond order.
  function lonePairs(z, bondOrderSum) {
    const ve = E[z].valenceElectrons;
    if (ve === null || E[z].isMetal) return 0;
    return Math.max(0, Math.floor((ve - bondOrderSum) / 2));
  }

  function geometry(z, neighbours, bondOrderSum) {
    const lp = lonePairs(z, bondOrderSum);
    const g = VSEPR[neighbours + ',' + Math.min(lp, 3)];
    if (!g) return { name: neighbours + '-coordinate', angle: 360 / Math.max(neighbours, 1), lonePairs: lp };
    let angle = g[1];
    // Heavier hydrides/bent molecules have angles close to 90 degrees (H2S 92, PH3 93.5).
    if (neighbours === 2 && lp === 2 && E[z].period > 2) angle = 92;
    if (neighbours === 3 && lp === 1 && E[z].period > 2) angle = 94;
    return { name: g[0], angle: angle, lonePairs: lp };
  }

  // Angle (radians) the simulation uses between the two bonds of a 2-coordinate centre.
  function twoBondAngle(z, bondOrderSum) {
    return geometry(z, 2, bondOrderSum).angle * Math.PI / 180;
  }

  // ---------------------------------------------------------------------------------------
  // Formulas & names
  // ---------------------------------------------------------------------------------------
  const SUB = '₀₁₂₃₄₅₆₇₈₉';
  function subscript(n) { return String(n).replace(/\d/g, function (d) { return SUB[+d]; }); }

  // counts: { z: n } -> Hill formula string (plain digits), used as a stable key.
  function hillFormula(counts) {
    const zs = Object.keys(counts).map(Number);
    let order;
    if (counts[6]) {
      order = [6].concat(counts[1] ? [1] : []).concat(
        zs.filter(function (z) { return z !== 6 && z !== 1; })
          .sort(function (a, b) { return E[a].symbol < E[b].symbol ? -1 : 1; }));
    } else {
      order = zs.sort(function (a, b) { return E[a].symbol < E[b].symbol ? -1 : 1; });
    }
    return order.map(function (z) { return E[z].symbol + (counts[z] > 1 ? counts[z] : ''); }).join('');
  }

  // Conventional display order when no name entry exists: organic -> Hill, inorganic -> by EN.
  function displayFormula(counts) {
    const zs = Object.keys(counts).map(Number);
    if (counts[6] && zs.length > 1) return pretty(hillFormula(counts));
    zs.sort(function (a, b) { return E[a].en - E[b].en || a - b; });
    return zs.map(function (z) { return E[z].symbol + (counts[z] > 1 ? subscript(counts[z]) : ''); }).join('');
  }

  function pretty(f) { return f.replace(/\d+/g, function (d) { return subscript(d); }); }

  // Hill key -> [display formula, name]
  const NAMES_RAW = `
H2|H2|Hydrogen gas;O2|O2|Oxygen gas;N2|N2|Nitrogen gas;F2|F2|Fluorine gas;Cl2|Cl2|Chlorine gas;
Br2|Br2|Bromine;I2|I2|Iodine;O3|O3|Ozone;S8|S8|Octasulfur;P4|P4|White phosphorus;S2|S2|Disulfur;
H2O|H2O|Water;H2O2|H2O2|Hydrogen peroxide;HO|OH|Hydroxyl radical;H3N|NH3|Ammonia;H4N2|N2H4|Hydrazine;
CH4|CH4|Methane;C2H6|C2H6|Ethane;C3H8|C3H8|Propane;C4H10|C4H10|Butane;C5H12|C5H12|Pentane;
C6H14|C6H14|Hexane;C8H18|C8H18|Octane;C2H4|C2H4|Ethylene;C2H2|C2H2|Acetylene;C3H6|C3H6|Propene;
C6H6|C6H6|Benzene;CH3|CH3|Methyl radical;CH2|CH2|Methylene;CH|CH|Methylidyne;
CO2|CO2|Carbon dioxide;CO|CO|Carbon monoxide;CH4O|CH3OH|Methanol;C2H6O|C2H5OH|Ethanol / dimethyl ether;
CH2O|CH2O|Formaldehyde;CH2O2|HCOOH|Formic acid;C2H4O2|CH3COOH|Acetic acid;C2H4O|C2H4O|Acetaldehyde;
CH2O3|H2CO3|Carbonic acid;C3H6O|C3H6O|Acetone;CH5N|CH3NH2|Methylamine;CH4N2O|CO(NH2)2|Urea;
CHN|HCN|Hydrogen cyanide;CS2|CS2|Carbon disulfide;CCl4|CCl4|Carbon tetrachloride;CHCl3|CHCl3|Chloroform;
CH3Cl|CH3Cl|Chloromethane;CF4|CF4|Tetrafluoromethane;C2N2|C2N2|Cyanogen;COS|OCS|Carbonyl sulfide;
CCl2O|COCl2|Phosgene;C6H12O6|C6H12O6|Glucose;
NO|NO|Nitric oxide;NO2|NO2|Nitrogen dioxide;N2O|N2O|Nitrous oxide;N2O4|N2O4|Dinitrogen tetroxide;
N2O3|N2O3|Dinitrogen trioxide;N2O5|N2O5|Dinitrogen pentoxide;HNO2|HNO2|Nitrous acid;HNO3|HNO3|Nitric acid;
O2S|SO2|Sulfur dioxide;O3S|SO3|Sulfur trioxide;H2S|H2S|Hydrogen sulfide;H2O4S|H2SO4|Sulfuric acid;
H2O3S|H2SO3|Sulfurous acid;F6S|SF6|Sulfur hexafluoride;F4S|SF4|Sulfur tetrafluoride;
ClH|HCl|Hydrogen chloride;FH|HF|Hydrogen fluoride;BrH|HBr|Hydrogen bromide;HI|HI|Hydrogen iodide;
ClHO|HOCl|Hypochlorous acid;ClF|ClF|Chlorine monofluoride;ClF3|ClF3|Chlorine trifluoride;
ClO2|ClO2|Chlorine dioxide;Cl2O|Cl2O|Dichlorine monoxide;BrCl|BrCl|Bromine monochloride;ClI|ICl|Iodine monochloride;
H3P|PH3|Phosphine;Cl3P|PCl3|Phosphorus trichloride;Cl5P|PCl5|Phosphorus pentachloride;F5P|PF5|Phosphorus pentafluoride;
H3O4P|H3PO4|Phosphoric acid;O10P4|P4O10|Phosphorus pentoxide;O6P4|P4O6|Phosphorus trioxide;
H4Si|SiH4|Silane;O2Si|SiO2|Silicon dioxide;Cl4Si|SiCl4|Silicon tetrachloride;F4Si|SiF4|Silicon tetrafluoride;
BF3|BF3|Boron trifluoride;BCl3|BCl3|Boron trichloride;BH3|BH3|Borane;B2H6|B2H6|Diborane;B2O3|B2O3|Boron trioxide;
F2Xe|XeF2|Xenon difluoride;F4Xe|XeF4|Xenon tetrafluoride;F6Xe|XeF6|Xenon hexafluoride;O3Xe|XeO3|Xenon trioxide;
F2Kr|KrF2|Krypton difluoride;
ClNa|NaCl|Sodium chloride (salt);ClK|KCl|Potassium chloride;ClLi|LiCl|Lithium chloride;FNa|NaF|Sodium fluoride;
BrNa|NaBr|Sodium bromide;INa|NaI|Sodium iodide;BrK|KBr|Potassium bromide;IK|KI|Potassium iodide;
CaCl2|CaCl2|Calcium chloride;Cl2Mg|MgCl2|Magnesium chloride;CaF2|CaF2|Calcium fluoride;BaCl2|BaCl2|Barium chloride;
HNaO|NaOH|Sodium hydroxide;HKO|KOH|Potassium hydroxide;HLiO|LiOH|Lithium hydroxide;CaH2O2|Ca(OH)2|Calcium hydroxide;
H2MgO2|Mg(OH)2|Magnesium hydroxide;HNa|NaH|Sodium hydride;HLi|LiH|Lithium hydride;HK|KH|Potassium hydride;
Na2O|Na2O|Sodium oxide;Na2O2|Na2O2|Sodium peroxide;K2O|K2O|Potassium oxide;Li2O|Li2O|Lithium oxide;
MgO|MgO|Magnesium oxide;CaO|CaO|Calcium oxide;BaO|BaO|Barium oxide;SrO|SrO|Strontium oxide;
Al2O3|Al2O3|Aluminium oxide;AlCl3|AlCl3|Aluminium chloride;Al2Cl6|Al2Cl6|Aluminium chloride dimer;
Fe2O3|Fe2O3|Iron(III) oxide (rust);FeO|FeO|Iron(II) oxide;Fe3O4|Fe3O4|Magnetite;Cl3Fe|FeCl3|Iron(III) chloride;
Cl2Fe|FeCl2|Iron(II) chloride;FeS|FeS|Iron(II) sulfide;FeS2|FeS2|Pyrite;
CuO|CuO|Copper(II) oxide;Cu2O|Cu2O|Copper(I) oxide;Cl2Cu|CuCl2|Copper(II) chloride;CuS|CuS|Copper sulfide;
OZn|ZnO|Zinc oxide;SZn|ZnS|Zinc sulfide;Cl2Zn|ZnCl2|Zinc chloride;O2Ti|TiO2|Titanium dioxide;Cl4Ti|TiCl4|Titanium tetrachloride;
AgCl|AgCl|Silver chloride;Ag2O|Ag2O|Silver oxide;Ag2S|Ag2S|Silver sulfide;OPb|PbO|Lead(II) oxide;
HgO|HgO|Mercury(II) oxide;Cl2Hg|HgCl2|Mercury(II) chloride;O2Sn|SnO2|Tin(IV) oxide;O2U|UO2|Uranium dioxide;
F6U|UF6|Uranium hexafluoride;O3W|WO3|Tungsten trioxide;MnO2|MnO2|Manganese dioxide;Cr2O3|Cr2O3|Chromium(III) oxide;
NiO|NiO|Nickel(II) oxide;CoO|CoO|Cobalt(II) oxide;Na2S|Na2S|Sodium sulfide;CaS|CaS|Calcium sulfide;
Li3N|Li3N|Lithium nitride;Mg3N2|Mg3N2|Magnesium nitride;C2Ca|CaC2|Calcium carbide;
CNa2O3|Na2CO3|Sodium carbonate;CCaO3|CaCO3|Calcium carbonate;CHNaO3|NaHCO3|Sodium bicarbonate;
NNaO3|NaNO3|Sodium nitrate;KNO3|KNO3|Potassium nitrate;Na2O4S|Na2SO4|Sodium sulfate;CaO4S|CaSO4|Calcium sulfate;
ClH4N|NH4Cl|Ammonium chloride;H2Se|H2Se|Hydrogen selenide;AsH3|AsH3|Arsine;GeH4|GeH4|Germane;
H4Sn|SnH4|Stannane;H3Sb|SbH3|Stibine;O2Se|SeO2|Selenium dioxide;Br3P|PBr3|Phosphorus tribromide
`;
  const NAMES = {};
  NAMES_RAW.replace(/\n/g, '').split(';').forEach(function (entry) {
    const p = entry.trim().split('|');
    if (p.length === 3) NAMES[p[0]] = { display: pretty(p[1]), name: p[2] };
  });

  // Describes a molecule/cluster from element counts and whether it contains ionic bonds.
  function describe(counts, atomCount, hasIonic) {
    const zs = Object.keys(counts).map(Number);
    if (atomCount === 1) {
      const el = E[zs[0]];
      return {
        key: el.symbol,
        formula: el.symbol,
        name: el.isMetal ? el.name + ' (metal)' : el.isNoble ? el.name : el.name + ' atom (radical)'
      };
    }
    const hill = hillFormula(counts);
    if (NAMES[hill]) return { key: hill, formula: NAMES[hill].display, name: NAMES[hill].name };
    if (atomCount > 12) {
      // Large networks/crystals/polymers: report the reduced (empirical) formula.
      let g = 0;
      zs.forEach(function (z) { g = gcd(g, counts[z]); });
      const red = {};
      zs.forEach(function (z) { red[z] = counts[z] / g; });
      const redHill = hillFormula(red);
      const base = NAMES[redHill] ? NAMES[redHill] : { display: displayFormula(red), name: null };
      const kind = zs.length === 1 ? 'network' : hasIonic ? 'ionic lattice' : 'network/polymer';
      return {
        key: '(' + redHill + ')n',
        formula: '(' + base.display + ')ₙ',
        name: (base.name ? base.name + ' ' : '') + kind
      };
    }
    return { key: hill, formula: displayFormula(counts), name: null };
  }

  function gcd(a, b) { while (b) { const t = b; b = a % b; a = t; } return a; }

  CE.chem = {
    NZ: NZ,
    pairClass: pairClass,
    pairMaxOrder: pairMaxOrder,
    pairEnergy: pairEnergy,
    pairLength: pairLength,
    bondClass: bondClass,
    bondEnergy: bondEnergy,
    stepEnergy: stepEnergy,
    maxOrder: maxOrder,
    geometry: geometry,
    lonePairs: lonePairs,
    twoBondAngle: twoBondAngle,
    hillFormula: hillFormula,
    describe: describe,
    subscript: subscript,
    pretty: pretty,
    NAMES: NAMES
  };
})(typeof window !== 'undefined' ? (window.CE = window.CE || {}) : (globalThis.CE = globalThis.CE || {}));
