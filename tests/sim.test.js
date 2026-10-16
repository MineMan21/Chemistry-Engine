// Headless checks of the element data, chemistry rules and simulation behaviour.
// Run with: npm test   (or: node tests/sim.test.js)
'use strict';
const assert = require('assert');
const CE = require('./load');

const P = CE.Simulation.params;
const C = CE.chem;
const Z = function (s) { return CE.bySymbol[s].z; };

// deterministic random numbers so runs are repeatable
function rng(seed) {
  return function () {
    seed |= 0; seed = seed + 0x6D2B79F5 | 0;
    let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}
function newSim(seed) { return new CE.Simulation({ rng: rng(seed || 1) }); }
function run(sim, seconds) { const steps = Math.round(seconds / P.DT); for (let s = 0; s < steps; s += 6) sim.step(6); }
function mol(sim, key) { const s = sim.analyze().find(function (x) { return x.key === key; }); return s ? s.mol : 0; }

const tests = [];
function test(name, fn) { tests.push({ name: name, fn: fn }); }

test('all 118 elements have complete data and unique table positions', function () {
  const seen = new Set();
  for (let z = 1; z <= 118; z++) {
    const e = CE.elements[z];
    assert.strictEqual(e.z, z);
    assert.ok(e.symbol && e.name && e.mass > 0 && e.covRadius > 0, e.symbol);
    assert.ok(typeof e.en === 'number', e.symbol + ' electronegativity');
    const pos = e.fIndex >= 0 ? 'f' + e.period + ':' + e.fIndex : e.period + ':' + e.group;
    assert.ok(!seen.has(pos), 'duplicate position ' + pos);
    seen.add(pos);
  }
});

test('bond classes follow electronegativity and metal rules', function () {
  assert.strictEqual(C.bondClass(Z('Na'), Z('Cl')), 2, 'NaCl ionic');
  assert.strictEqual(C.bondClass(Z('H'), Z('O')), 1, 'O-H covalent');
  assert.strictEqual(C.bondClass(Z('Fe'), Z('Fe')), 0, 'metals do not form discrete bonds');
  assert.strictEqual(C.bondClass(Z('He'), Z('F')), 0, 'helium is inert');
  assert.strictEqual(C.maxOrder(Z('N'), Z('N')), 3, 'N triple bond');
  assert.strictEqual(C.maxOrder(Z('O'), Z('O')), 2, 'O double bond');
  assert.strictEqual(C.maxOrder(Z('S'), Z('S')), 1, 'sulfur chains are single bonded');
});

test('bond energies and VSEPR geometry', function () {
  assert.strictEqual(C.bondEnergy(Z('H'), Z('O'), 1), 459);
  assert.strictEqual(C.bondEnergy(Z('N'), Z('N'), 3), 942);
  const water = C.geometry(Z('O'), 2, 2);
  assert.strictEqual(water.name, 'bent');
  assert.strictEqual(water.angle, 104.5);
  assert.strictEqual(C.geometry(Z('C'), 4, 4).name, 'tetrahedral');
  assert.strictEqual(C.geometry(Z('C'), 2, 4).name, 'linear'); // CO2
});

test('formulas and names', function () {
  const w = C.describe({ 1: 2, 8: 1 }, 3, false);
  assert.strictEqual(w.name, 'Water');
  assert.strictEqual(w.formula, 'H₂O');
  assert.strictEqual(C.describe({ 11: 1, 17: 1 }, 2, true).formula, 'NaCl');
  assert.strictEqual(C.describe({ 6: 1, 8: 2 }, 3, false).name, 'Carbon dioxide');
  assert.strictEqual(C.hillFormula({ 6: 2, 1: 6, 8: 1 }), 'C2H6O');
});

test('adding by the mole: 1 mol of nitrogen is 1000 N2 molecules', function () {
  const sim = newSim(2);
  const res = sim.addMoles(Z('N'), 1, P.W / 2, P.H / 2);
  assert.strictEqual(res.units, 1000);
  assert.strictEqual(sim.n, 2000);
  assert.ok(Math.abs(mol(sim, 'N2') - 1) < 1e-9);
});

test('1 mol of N2 at room temperature reads about 1 atm and does not react', function () {
  const sim = newSim(3);
  sim.addMoles(Z('N'), 1, P.W / 2, P.H / 2);
  run(sim, 8);
  assert.strictEqual(sim.reactionCount, 0);
  assert.ok(sim.temperature > 240 && sim.temperature < 340, 'T ' + sim.temperature);
  assert.ok(sim.pressureAtm > 0.5 && sim.pressureAtm < 1.7, 'P ' + sim.pressureAtm);
});

test('hydrogen + oxygen is stable at room temperature', function () {
  const sim = newSim(4);
  sim.addMoles(Z('H'), 0.6, P.W / 2, P.H * 0.35);
  sim.addMoles(Z('O'), 0.3, P.W / 2, P.H * 0.65);
  run(sim, 6);
  assert.strictEqual(mol(sim, 'H2O'), 0);
});

test('hydrogen burns to water on a hot plate', function () {
  const sim = newSim(5);
  sim.addMoles(Z('H'), 0.6, P.W / 2, P.H * 0.35);
  sim.addMoles(Z('O'), 0.3, P.W / 2, P.H * 0.65);
  sim.plateTemp = 2000;
  run(sim, 15);
  assert.ok(mol(sim, 'H2O') > 0.1, 'water ' + mol(sim, 'H2O'));
});

test('a spark burns a 2:1 hydrogen/oxygen mix through in a few seconds', function () {
  const sim = newSim(11);
  sim.addMoles(Z('H'), 1, P.W / 2, P.H * 0.35);
  sim.addMoles(Z('O'), 0.5, P.W / 2, P.H * 0.65);
  run(sim, 2);
  assert.ok(mol(sim, 'H2O') < 0.01, 'no water before the spark');
  for (let f = 0; f < 10; f++) { sim.spark(P.W / 2, P.H * 0.6, 380, 4000); sim.step(9); }
  run(sim, 4);
  // nearly all the oxygen is used up (the rest is radicals still recombining)
  assert.ok(mol(sim, 'O2') < 0.05, 'O2 left ' + mol(sim, 'O2'));
  assert.ok(mol(sim, 'H2O') > 0.6, 'water ' + mol(sim, 'H2O'));
});

test('sodium and chlorine form sodium chloride', function () {
  const sim = newSim(6);
  sim.addMoles(Z('Na'), 0.4, P.W / 2, P.H * 0.85);
  sim.addMoles(Z('Cl'), 0.3, P.W / 2, P.H * 0.4);
  run(sim, 10);
  const salt = sim.analyze().filter(function (s) { return s.key === 'ClNa' || s.key === '(ClNa)n'; });
  assert.ok(salt.length > 0, 'no NaCl formed');
});

test('water is liquid at 300 K and boils away on a hot plate', function () {
  const H2O = { atoms: [[8, 0, 0], [1, 9.7, 0], [1, -2.43, 9.39]], bonds: [[0, 1, 1], [0, 2, 1]] };
  function condensed(T) {
    const sim = newSim(7);
    sim.plateTemp = T;
    sim.addMolecules(H2O, 400, P.W / 2, P.H - 150, { x: 600, y: 140 });
    sim.react = function () {}; sim.bondEvents = function () {}; // phase behaviour only
    run(sim, 6);
    const w = sim.analyze().find(function (s) { return s.key === 'H2O'; });
    return w.condensed / w.molecules;
  }
  const cold = condensed(300), hot = condensed(520);
  assert.ok(cold > 0.8, 'liquid fraction at 300 K ' + cold);
  assert.ok(hot < cold - 0.2, 'liquid fraction at 520 K ' + hot);
});

test('liquid water falls under gravity and wets the glass', function () {
  const H2O = { atoms: [[8, 0, 0], [1, 9.7, 0], [1, -2.43, 9.39]], bonds: [[0, 1, 1], [0, 2, 1]] };
  const sim = newSim(10);
  sim.react = function () {}; sim.bondEvents = function () {};
  sim.addMolecules(H2O, 500, P.W / 2, P.H * 0.3, { x: 250, y: 150 });
  function meanY() { let s = 0, c = 0; for (let i = 0; i < sim.n; i++) if (sim.el[i] === 8) { s += sim.y[i]; c++; } return s / c; }
  const y0 = meanY();
  run(sim, 6);
  assert.ok(meanY() > y0 + 1000, 'water fell ' + (meanY() - y0) + ' units');
  assert.ok(sim.temperature < 380, 'falling did not overheat it (' + sim.temperature + ' K)');
  // a thin film placed against the glass stays there
  const film = newSim(11);
  film.react = function () {}; film.bondEvents = function () {};
  const placed = film.addMolecules(H2O, 200, 30, P.H * 0.35, { x: 18, y: 600 });
  run(film, 3);
  let atWall = 0;
  for (let i = 0; i < film.n; i++) if (film.el[i] === 8 && film.x[i] < 60 && film.y[i] < P.H * 0.6) atWall++;
  assert.ok(atWall >= 0.6 * placed, 'molecules still clinging to the glass: ' + atWall + ' of ' + placed);
});

test('the vacuum pump removes gas and the gas can adds it', function () {
  const sim = newSim(8);
  sim.addMoles(Z('Ar'), 0.5, P.W / 2, P.H / 2);
  sim.vacuumOn = true;
  run(sim, 4);
  assert.ok(sim.n < 500, 'argon left ' + sim.n);
  sim.vacuumOn = false;
  sim.canGas = 'N2'; sim.canRate = 1; sim.canOn = true;
  const before = sim.n;
  run(sim, 1);
  assert.ok(sim.n > before + 1000, 'gas can added ' + (sim.n - before));
});

test('phase data: documented points and phases', function () {
  const ph = CE.phase, w = ph.forKey('H2O');
  assert.strictEqual(ph.phaseAt(w, 293, 1), 'liquid');
  assert.strictEqual(ph.phaseAt(w, 263, 1), 'solid');
  assert.strictEqual(ph.phaseAt(w, 400, 1), 'gas');
  assert.strictEqual(ph.phaseAt(w, 700, 300), 'supercritical');
  assert.strictEqual(ph.phaseAt(w, 100, 1e-9), 'solid', 'ice at 100 K does not evaporate in a vacuum');
  // the vapour curve passes through the normal boiling point; ice's vapour pressure at 250 K is ~7.6e-4 atm
  assert.ok(Math.abs(ph.vaporP(w, 373.12) - 1) < 1e-6);
  const p250 = ph.sublimationP(w, 250);
  assert.ok(p250 > 5e-4 && p250 < 1.1e-3, 'ice vapour pressure at 250 K ' + p250);
  // dry ice sublimes at 1 atm; liquid CO2 needs more than 5.1 atm
  const c = ph.forKey('CO2');
  assert.strictEqual(ph.phaseAt(c, 190, 1), 'solid');
  assert.strictEqual(ph.phaseAt(c, 210, 1), 'gas');
  assert.strictEqual(ph.phaseAt(c, 250, 30), 'liquid');
  assert.ok(ph.forElement(Z('O')) === ph.forKey('O2') && ph.forElement(Z('Ar')) === ph.forKey('Ar'));
  assert.strictEqual(ph.forElement(Z('Fe')), null);
  Object.keys(ph.DATA).forEach(function (k) {
    const d = ph.DATA[k], g = ph.diagram(d);
    g.lines.forEach(function (l) { l.pts.forEach(function (p) { assert.ok(isFinite(p[0]) && isFinite(p[1]) && p[1] > 0, k + ' ' + l.kind); }); });
  });
});

// honeycomb (3-fold) order of water oxygens: 1 when every oxygen has its neighbours at 120 degrees
function waterOrder(sim, onlyCrystal) {
  const all = [];
  for (let i = 0; i < sim.n; i++) if (sim.isWaterO(i)) all.push(i);
  // the crystal's own structure: crystal members and their crystal neighbours
  const O = onlyCrystal ? all.filter(function (i) { return sim.grainOf[i] >= 0; }) : all;
  let order = 0, counted = 0;
  O.forEach(function (i) {
    let re = 0, im = 0, k = 0;
    O.forEach(function (j) {
      if (j === i) return;
      const dx = sim.x[j] - sim.x[i], dy = sim.y[j] - sim.y[i];
      if (dx * dx + dy * dy < 36 * 36) { const a = 3 * Math.atan2(dy, dx); re += Math.cos(a); im += Math.sin(a); k++; }
    });
    if (k >= 2) { order += Math.hypot(re, im) / k; counted++; }
  });
  return { O: all, psi3: order / counted };
}

test('water is a disordered liquid at 300 K and freezes into a honeycomb crystal at 200 K', function () {
  const H2O = { atoms: [[8, 0, 0], [1, 9.7, 0], [1, -2.43, 9.39]], bonds: [[0, 1, 1], [0, 2, 1]] };
  const sim = newSim(12);
  sim.react = function () {}; sim.bondEvents = function () {};
  sim.addMolecules(H2O, 500, P.W / 2, P.H - 200, { x: 350, y: 150 });
  sim.plateTemp = 300;
  run(sim, 4);
  const liquid = waterOrder(sim).psi3;
  assert.ok(liquid < 0.6, 'liquid water is disordered: honeycomb order ' + liquid.toFixed(2));
  sim.plateTemp = 200;
  run(sim, 16);
  const ice = waterOrder(sim, true), O = ice.O;
  let solid = 0, crystal = 0;
  O.forEach(function (i) { if (sim.isSolidAtom(i)) solid++; if (sim.grainOf[i] >= 0) crystal++; });
  assert.ok(sim.temperature < 260, 'cooled to ' + sim.temperature);
  assert.ok(solid > 0.8 * O.length, 'ice counted as solid: ' + solid + ' of ' + O.length);
  assert.ok(crystal > 0.5 * O.length, 'ice is crystalline: ' + crystal + ' of ' + O.length + ' in crystal grains');
  assert.ok(ice.psi3 > 0.75, 'ice crystal is a honeycomb: order ' + ice.psi3.toFixed(2));
});

// order and mobility of the atoms of element z: hexagonal order psi6, fraction in crystal grains,
// fraction with no metal neighbour (vapour) and the 90th-percentile height above the floor
function metalState(sim, z) {
  const A = [];
  for (let i = 0; i < sim.n; i++) if (sim.el[i] === z) A.push(i);
  const sig = 2 * CE.elements[z].radius;
  let order = 0, counted = 0, crystal = 0, vapour = 0;
  const hs = [];
  A.forEach(function (i) {
    let re = 0, im = 0, k = 0, near = 0;
    A.forEach(function (j) {
      if (j === i) return;
      const dx = sim.x[j] - sim.x[i], dy = sim.y[j] - sim.y[i], d2 = dx * dx + dy * dy;
      if (d2 < sig * sig * 2.1) near++;
      if (d2 < sig * sig * 1.56) { const a = 6 * Math.atan2(dy, dx); re += Math.cos(a); im += Math.sin(a); k++; }
    });
    if (k >= 2) { order += Math.hypot(re, im) / k; counted++; }
    if (sim.grainOf[i] >= 0) crystal++;
    if (!near) vapour++;
    hs.push(P.H - sim.y[i]);
  });
  hs.sort(function (a, b) { return a - b; });
  return { psi6: order / counted, crystal: crystal / A.length, vapour: vapour / A.length, h90: hs[Math.floor(hs.length * 0.9)] };
}
function positions(sim) { return { x: Array.from(sim.x.slice(0, sim.n)), y: Array.from(sim.y.slice(0, sim.n)) }; }
function rmsMove(sim, before) {
  let s = 0;
  for (let i = 0; i < before.x.length; i++) s += Math.pow(sim.x[i] - before.x[i], 2) + Math.pow(sim.y[i] - before.y[i], 2);
  return Math.sqrt(s / before.x.length);
}

test('iron melts on a 2100 K plate, stays a liquid, and hardens into crystal grains when cooled', function () {
  const sim = newSim(21);
  const fe = Z('Fe');
  sim.addMoles(fe, 0.3, P.W / 2, P.H - 150);
  run(sim, 2);
  const cold = metalState(sim, fe);
  assert.ok(cold.crystal > 0.8 && cold.psi6 > 0.85, 'solid iron at room temperature is crystalline: ' + JSON.stringify(cold));
  sim.plateTemp = 2100;
  run(sim, 8);
  let before = positions(sim);
  run(sim, 1);
  const hot = metalState(sim, fe), flow = rmsMove(sim, before);
  assert.ok(hot.crystal < 0.1, 'molten iron is not crystalline: ' + JSON.stringify(hot));
  assert.ok(flow > 25, 'molten iron flows: atoms moved ' + flow.toFixed(1) + ' per second');
  assert.ok(hot.vapour < 0.03 && hot.h90 < 600, 'molten iron stays together on the plate: ' + JSON.stringify(hot));
  sim.plateTemp = 293;
  run(sim, 12);
  before = positions(sim);
  run(sim, 1);
  const recast = metalState(sim, fe);
  assert.ok(recast.crystal > 0.6 && recast.psi6 > 0.8, 'cooled iron crystallises: ' + JSON.stringify(recast));
  const still = rmsMove(sim, before);
  assert.ok(still < flow / 2, 'cooled iron is rigid: atoms moved ' + still.toFixed(1) + ' per second (melt: ' + flow.toFixed(1) + ')');
});

test('gallium melts just above its 303 K melting point and freezes below it', function () {
  const sim = newSim(22);
  const ga = Z('Ga');
  sim.addMoles(ga, 0.3, P.W / 2, P.H - 150);
  sim.plateTemp = 340;
  run(sim, 8);
  let before = positions(sim);
  run(sim, 1);
  const melt = metalState(sim, ga);
  assert.ok(melt.crystal < 0.1 && rmsMove(sim, before) > 15, 'liquid gallium at 340 K: ' + JSON.stringify(melt));
  assert.ok(melt.vapour < 0.03, 'liquid gallium does not evaporate: ' + JSON.stringify(melt));
  sim.plateTemp = 240;
  run(sim, 12);
  const frozen = metalState(sim, ga);
  assert.ok(frozen.crystal > 0.6, 'gallium freezes at 240 K: ' + JSON.stringify(frozen));
});

test('mercury is a liquid at room temperature', function () {
  const sim = newSim(23);
  const hg = Z('Hg');
  sim.addMoles(hg, 0.3, P.W / 2, P.H - 150);
  run(sim, 6);
  const before = positions(sim);
  run(sim, 1);
  const s = metalState(sim, hg);
  assert.ok(s.crystal < 0.1 && rmsMove(sim, before) > 15, 'mercury flows at 293 K: ' + JSON.stringify(s));
  assert.ok(s.vapour < 0.05, 'mercury stays mostly liquid: ' + JSON.stringify(s));
});

test('the grab tool lifts a solid iron lump in one piece from a small grip', function () {
  const sim = newSim(24);
  const fe = Z('Fe');
  sim.addMoles(fe, 0.3, P.W / 2, P.H - 150);
  run(sim, 3);
  let top = P.H, cx = 0, c = 0;
  for (let i = 0; i < sim.n; i++) top = Math.min(top, sim.y[i]);
  for (let i = 0; i < sim.n; i++) if (sim.y[i] < top + 60) { cx += sim.x[i]; c++; }
  cx /= c;
  // a small grip on a solid takes the whole connected solid
  assert.strictEqual(sim.grabStart(cx, top + 70, 40), sim.n, 'the whole solid lump is held');
  const steps = Math.round(3 / P.DT);
  for (let s = 0; s < steps; s += 6) { sim.grabMove(cx, top + 70 - 900 * s / steps); sim.step(6); }
  run(sim, 1);
  let lifted = 0;
  for (let i = 0; i < sim.n; i++) if (sim.y[i] < top - 200) lifted++;
  assert.ok(lifted > 0.95 * sim.n, 'the whole lump came up: ' + lifted + ' of ' + sim.n);
  sim.grabEnd();
  assert.strictEqual(sim.grab, null);
});

test('a held solid keeps its shape when shoved into the glass', function () {
  const sim = newSim(25);
  const fe = Z('Fe');
  sim.addMoles(fe, 0.3, P.W / 2, P.H - 150);
  run(sim, 3);
  const A = [];
  for (let i = 0; i < sim.n; i++) A.push(i);
  const pairs = [];
  A.forEach(function (i) {
    A.forEach(function (j) {
      const d = Math.hypot(sim.x[j] - sim.x[i], sim.y[j] - sim.y[i]);
      if (j > i && d < 2 * CE.elements[fe].radius * 1.2) pairs.push([i, j, d]);
    });
  });
  let top = P.H, cx = 0;
  A.forEach(function (i) { top = Math.min(top, sim.y[i]); cx += sim.x[i] / A.length; });
  sim.grabStart(cx, top + 30, 120);
  const steps = Math.round(4 / P.DT);
  for (let s = 0; s < steps; s += 6) {
    const t = s / steps;
    // lift, then drag the hand far past the left wall
    sim.grabMove(t < 0.5 ? cx : cx - 3000 * (t - 0.5), top + 30 - 600 * Math.min(1, 2 * t));
    sim.step(6);
  }
  let broken = 0;
  pairs.forEach(function (p) { if (Math.hypot(sim.x[p[1]] - sim.x[p[0]], sim.y[p[1]] - sim.y[p[0]]) > 1.4 * p[2]) broken++; });
  let left = P.W;
  A.forEach(function (i) { left = Math.min(left, sim.x[i]); });
  assert.ok(left < 60, 'the lump reached the wall (leftmost ' + left.toFixed(0) + ')');
  assert.ok(broken < 0.02 * pairs.length, 'the lump is not crushed: ' + broken + ' of ' + pairs.length + ' contacts broken');
  sim.grabEnd();
});

// ---- textbook reactions ----------------------------------------------------------------
// atoms of metal z bonded to element y, and atoms of y bonded to z: the composition of the
// compound the two have made
function compound(sim, z, y) {
  let nz = 0, ny = 0, zTotal = 0, other = 0;
  for (let i = 0; i < sim.n; i++) {
    const e = sim.el[i];
    if (e !== z && e !== y) continue;
    if (e === z) zTotal++;
    let toPartner = false, toOther = false;
    for (let k = 0; k < sim.nb[i]; k++) {
      const f = sim.el[sim.bp[i * P.MAXB + k]];
      if (f === (e === z ? y : z)) toPartner = true;
      else if (e === z) toOther = true;
    }
    if (toPartner) { if (e === z) nz++; else ny++; }
    if (toOther) other++;
  }
  return { reacted: nz / zTotal, ratio: ny / Math.max(1, nz), metalOtherBonds: other };
}
function speciesWith(sim, test) {
  return sim.analyze().filter(function (s) { return test(s.counts || {}); });
}

test('iron is untouched by oxygen at room temperature and oxidises to Fe2O3/Fe3O4 when heated', function () {
  const fe = Z('Fe'), o = Z('O');
  const cold = newSim(31);
  cold.addMoles(fe, 0.2, P.W / 2, P.H - 150); cold.addMoles(o, 0.2, P.W / 2, P.H * 0.4);
  run(cold, 6);
  assert.ok(compound(cold, fe, o).reacted < 0.02, 'no rusting in dry air at 293 K: ' + JSON.stringify(compound(cold, fe, o)));
  const hot = newSim(32);
  hot.addMoles(fe, 0.2, P.W / 2, P.H - 150); hot.addMoles(o, 0.2, P.W / 2, P.H * 0.4);
  hot.plateTemp = 1100;
  run(hot, 14);
  const c = compound(hot, fe, o);
  assert.ok(c.reacted > 0.1, 'iron oxidises at 1100 K: ' + JSON.stringify(c));
  // hematite Fe2O3 is O:Fe 1.5, magnetite Fe3O4 1.33, wustite FeO 1; a thin scale on the metal
  // reads a little O-rich (its outermost oxygens are bound to one iron only)
  assert.ok(c.ratio > 1.0 && c.ratio < 1.7, 'iron oxide composition O:Fe ' + c.ratio.toFixed(2));
  // no peroxide or superoxide (O-O bonds) on the iron
  let peroxide = 0;
  for (let i = 0; i < hot.n; i++) {
    if (hot.el[i] !== o) continue;
    let toFe = false, toO = false;
    for (let k = 0; k < hot.nb[i]; k++) { const f = hot.el[hot.bp[i * P.MAXB + k]]; if (f === fe) toFe = true; if (f === o) toO = true; }
    if (toFe && toO) peroxide++;
  }
  assert.ok(peroxide === 0, 'no iron peroxides or superoxides: ' + peroxide + ' O-O bonds on iron');
});

test('wet iron rusts to hydroxides and oxides, not iron hydrides', function () {
  const H2O = { atoms: [[8, 0, 0], [1, 9.7, 0], [1, -2.43, 9.39]], bonds: [[0, 1, 1], [0, 2, 1]] };
  const fe = Z('Fe'), o = Z('O');
  const sim = newSim(33);
  sim.addMoles(fe, 0.2, P.W / 2, P.H - 150); sim.addMoles(o, 0.2, P.W / 2, P.H * 0.4);
  sim.addMolecules(H2O, 200, P.W / 2, P.H - 400, { x: 300, y: 100 });
  sim.plateTemp = 350;
  run(sim, 15);
  const hydride = speciesWith(sim, function (k) { return k[fe] && k[1] && !k[o]; });
  assert.ok(hydride.length === 0, 'no iron hydrides: ' + hydride.map(function (s) { return s.key; }).join(' '));
  for (let i = 0; i < sim.n; i++) {
    if (sim.el[i] !== fe) continue;
    for (let k = 0; k < sim.nb[i]; k++) assert.ok(sim.el[sim.bp[i * P.MAXB + k]] !== 1, 'no Fe-H bond');
  }
  assert.ok(compound(sim, fe, o).reacted > 0.02, 'iron corrodes in wet air: ' + JSON.stringify(compound(sim, fe, o)));
});

test('magnesium and aluminium burn to MgO and Al2O3', function () {
  [['Mg', 1.0], ['Al', 1.5]].forEach(function (c, n) {
    const m = Z(c[0]), o = Z('O');
    const sim = newSim(34 + n);
    sim.addMoles(m, 0.2, P.W / 2, P.H - 150); sim.addMoles(o, 0.2, P.W / 2, P.H * 0.4);
    sim.plateTemp = 1300;
    run(sim, 12);
    const r = compound(sim, m, o);
    assert.ok(r.reacted > 0.08, c[0] + ' burns: ' + JSON.stringify(r));
    assert.ok(Math.abs(r.ratio - c[1]) < 0.3, c[0] + ' oxide O:' + c[0] + ' ' + r.ratio.toFixed(2) + ' (expected ' + c[1] + ')');
  });
});

test('sodium and chlorine make solid sodium chloride', function () {
  const na = Z('Na'), cl = Z('Cl');
  const sim = newSim(36);
  sim.addMoles(na, 0.2, P.W / 2, P.H - 150); sim.addMoles(cl, 0.2, P.W / 2, P.H * 0.4);
  run(sim, 25); // the reaction heat first throws the new salt off the metal; it settles as a solid
  const r = compound(sim, na, cl);
  assert.ok(r.reacted > 0.08 && Math.abs(r.ratio - 1) < 0.1, 'NaCl forms 1:1: ' + JSON.stringify(r));
  const salt = sim.analyze().find(function (s) { return s.key === 'ClNa'; });
  assert.ok(salt && salt.condensed / salt.molecules > 0.5, 'salt is a solid: condensed ' + (salt ? (salt.condensed / salt.molecules).toFixed(2) : 'none'));
});

test('hydrogen and chlorine stay apart until sparked, then make HCl; N2 + H2 needs a catalyst', function () {
  const sim = newSim(37);
  sim.addMoles(Z('H'), 0.2, P.W / 2, P.H * 0.4); sim.addMoles(Z('Cl'), 0.2, P.W / 2, P.H * 0.6);
  run(sim, 2);
  assert.ok(mol(sim, 'ClH') < 0.01, 'dark H2/Cl2 mixture is stable');
  sim.spark(P.W / 2, P.H / 2, 550, 4000);
  run(sim, 8);
  assert.ok(mol(sim, 'ClH') > 0.2, 'HCl after the spark: ' + mol(sim, 'ClH').toFixed(3) + ' mol');
  const haber = newSim(38);
  haber.addMoles(Z('N'), 0.2, P.W / 2, P.H * 0.4); haber.addMoles(Z('H'), 0.6, P.W / 2, P.H * 0.6);
  run(haber, 6);
  assert.ok(mol(haber, 'H3N') < 0.005, 'no ammonia without a catalyst');
});

test('removing particles keeps bonds consistent', function () {
  const sim = newSim(9);
  sim.addMoles(Z('O'), 0.1, P.W / 2, P.H / 2);
  for (let k = 0; k < 50; k++) sim.removeMolecule(Math.floor(sim.n / 2));
  for (let i = 0; i < sim.n; i++) {
    for (let k = 0; k < sim.nb[i]; k++) {
      const j = sim.bp[i * P.MAXB + k];
      assert.ok(j >= 0 && j < sim.n, 'dangling bond');
      assert.ok(sim.isBonded(j, i), 'bond not mirrored');
    }
  }
});

(function main() {
  let failed = 0;
  const t0 = Date.now();
  tests.forEach(function (t) {
    const s = Date.now();
    try {
      t.fn();
      console.log('  ok   ' + t.name + ' (' + (Date.now() - s) + ' ms)');
    } catch (e) {
      failed++;
      console.log('  FAIL ' + t.name + '\n       ' + (e && e.message));
    }
  });
  console.log((tests.length - failed) + '/' + tests.length + ' passed in ' + ((Date.now() - t0) / 1000).toFixed(1) + ' s');
  process.exit(failed ? 1 : 0);
})();
