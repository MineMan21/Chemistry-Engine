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
  sim.startBurst(P.W / 2, P.H * 0.6, 550, 4000, 0.16); // the spark tool: heat only
  run(sim, 6);
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
    run(sim, 4);
    let f = 0;
    for (let k = 0; k < 8; k++) { run(sim, 0.5); const w = sim.analyze().find(function (s) { return s.key === 'H2O'; }); f += w.condensed / w.molecules / 8; }
    return f;
  }
  const cold = condensed(300), hot = condensed(520);
  // (real water in this jar would be ~91% liquid; the model's water is somewhat too volatile)
  assert.ok(cold > 0.72, 'liquid fraction at 300 K ' + cold);
  assert.ok(hot < cold - 0.2, 'liquid fraction at 520 K ' + hot);
});

test('hot steam in an insulated jar does not heat itself as it splits', function () {
  // Splitting water takes heat (and re-forming it gives the same back), so steam kicked hot and
  // left alone in an insulated jar cannot run away. (Bonds made between atoms a little apart
  // used to snap shut and fling the atoms off, and the steam ran away to ~19,000 K.)
  const H2O = { atoms: [[8, 0, 0], [1, 9.7, 0], [1, -2.43, 9.39]], bonds: [[0, 1, 1], [0, 2, 1]] };
  const sim = newSim(3);
  sim.addMolecules(H2O, 400, P.W / 2, P.H / 2, { x: 1000, y: 1400 });
  run(sim, 2);
  const saved = [P.PLATE_GAMMA, P.CONDUCTION];
  sim.insulated = true; sim.convection = 0; P.PLATE_GAMMA = 0; P.CONDUCTION = 0;
  try {
    for (let i = 0; i < sim.n; i++) {
      const sd = Math.sqrt(P.KB * 6000 / CE.elements[sim.el[i]].mass) * P.DT;
      sim.px[i] = sim.x[i] - (sim.rng() - 0.5) * 3.46 * sd; sim.py[i] = sim.y[i] - (sim.rng() - 0.5) * 3.46 * sd;
    }
    let peak = 0;
    for (let k = 0; k < 8; k++) { run(sim, 0.5); peak = Math.max(peak, sim.temperature); }
    assert.ok(peak < 10000, 'insulated steam heated itself to ' + Math.round(peak) + ' K');
  } finally { P.PLATE_GAMMA = saved[0]; P.CONDUCTION = saved[1]; }
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
  const placed = film.addMolecules(H2O, 200, 30, P.H * 0.35, { x: 18, y: 600, fixed: true });
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
  // (gallium in this model freezes some way below its real 303 K, see Known_Issues.md; metal
  // grains only record order: a frozen lump is close-packed, much of it in grains, with grain
  // boundaries between)
  sim.plateTemp = 200;
  run(sim, 12);
  const frozen = metalState(sim, ga);
  assert.ok(frozen.psi6 > 0.82 && frozen.crystal > 0.35, 'gallium freezes at 200 K: ' + JSON.stringify(frozen));
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

// a lump of solid iron and the pairs of atoms in contact in it
function ironLump(seed) {
  const sim = newSim(seed);
  const fe = Z('Fe');
  sim.addMoles(fe, 0.3, P.W / 2, P.H - 150);
  run(sim, 4);
  const pairs = [];
  for (let i = 0; i < sim.n; i++) for (let j = i + 1; j < sim.n; j++) {
    const d = Math.hypot(sim.x[j] - sim.x[i], sim.y[j] - sim.y[i]);
    if (d < 2 * CE.elements[fe].radius * 1.2) pairs.push([i, j, d]);
  }
  let top = P.H, cx = 0;
  for (let i = 0; i < sim.n; i++) { top = Math.min(top, sim.y[i]); cx += sim.x[i] / sim.n; }
  return { sim: sim, pairs: pairs, top: top, cx: cx,
    broken: function () { let b = 0; pairs.forEach(function (p) { if (Math.hypot(sim.x[p[1]] - sim.x[p[0]], sim.y[p[1]] - sim.y[p[0]]) > 1.4 * p[2]) b++; }); return b; } };
}

test('solid iron is strong: grabbed by a small patch the whole lump comes along', function () {
  const L = ironLump(24), sim = L.sim;
  const held = sim.grabStart(L.cx, L.top + 20, 120);
  assert.ok(held < 0.2 * sim.n, 'only a patch is under the cursor (' + held + ' of ' + sim.n + ')');
  const y0 = Array.from(sim.y.subarray(0, sim.n));
  const steps = Math.round(2 / P.DT);
  for (let s = 0; s < steps; s += 6) { sim.grabMove(L.cx, L.top + 20 - 700 * Math.min(1, s / (0.75 * steps))); sim.step(6); }
  let lifted = 0; // atoms that rose with the hand (it rose 700)
  for (let i = 0; i < sim.n; i++) if (sim.y[i] < y0[i] - 500) lifted++;
  assert.ok(lifted > 0.85 * sim.n, 'the lump came up with the patch: ' + lifted + ' of ' + sim.n);
  assert.ok(L.broken() < 0.03 * L.pairs.length, 'the lattice held: ' + L.broken() + ' of ' + L.pairs.length + ' contacts broken');
  sim.grabEnd();
});

test('a solid iron lump thrown against the glass stays whole; molten iron does not hold together', function () {
  const L = ironLump(25), sim = L.sim;
  sim.grabStart(L.cx, L.top + 20, 120);
  let steps = Math.round(1.5 / P.DT);
  for (let s = 0; s < steps; s += 6) { sim.grabMove(L.cx, L.top + 20 - 700 * Math.min(1, s / steps)); sim.step(6); }
  // fling it at the left wall and let go
  steps = Math.round(0.4 / P.DT);
  for (let s = 0; s < steps; s += 6) { sim.grabMove(L.cx - (L.cx - 300) * s / steps, L.top - 680); sim.step(6); }
  sim.grabEnd();
  run(sim, 4);
  // it hits the glass at about the speed of a rifle bullet in the real-world jar: it may dent
  // (contacts rearrange where it struck) but it must stay in one piece
  const parent = new Int32Array(sim.n);
  for (let i = 0; i < sim.n; i++) parent[i] = i;
  const find = function (a) { while (parent[a] !== a) { parent[a] = parent[parent[a]]; a = parent[a]; } return a; };
  const reach = 2.8 * CE.elements[Z('Fe')].radius;
  for (let i = 0; i < sim.n; i++) for (let j = i + 1; j < sim.n; j++) if (Math.hypot(sim.x[j] - sim.x[i], sim.y[j] - sim.y[i]) < reach) parent[find(i)] = find(j);
  const size = new Map();
  for (let i = 0; i < sim.n; i++) { const q = find(i); size.set(q, (size.get(q) || 0) + 1); }
  const piece = Math.max.apply(null, Array.from(size.values())) / sim.n;
  assert.ok(piece > 0.98 && L.broken() < 0.15 * L.pairs.length, 'the lump survived the throw: largest piece ' + piece.toFixed(2) + ', ' + L.broken() + ' of ' + L.pairs.length + ' contacts broken');
  // the same lump molten: a grab pulls a piece out instead of lifting it
  const M = newSim(26), fe = Z('Fe');
  M.addMoles(fe, 0.3, P.W / 2, P.H - 150);
  M.plateTemp = 2100;
  run(M, 10);
  let top = P.H, cx = 0;
  for (let i = 0; i < M.n; i++) { top = Math.min(top, M.y[i]); cx += M.x[i] / M.n; }
  M.grabStart(cx, top + 20, 120);
  steps = Math.round(2 / P.DT);
  for (let s = 0; s < steps; s += 6) { M.grabMove(cx, top + 20 - 700 * Math.min(1, s / (0.75 * steps))); M.step(6); }
  let lifted = 0;
  for (let i = 0; i < M.n; i++) if (M.y[i] < top - 300) lifted++;
  assert.ok(lifted < 0.5 * M.n, 'molten iron is not lifted as a lump: ' + lifted + ' of ' + M.n);
  M.grabEnd();
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

// Rusting at room temperature is electrochemical and takes days; in seconds, iron and water react
// only hot, as in the steam-iron process (3 Fe + 4 H2O -> Fe3O4 + 4 H2 above about 840 K)
test('iron in hot steam and air oxidises (releasing hydrogen), not to iron hydrides', function () {
  const H2O = { atoms: [[8, 0, 0], [1, 9.7, 0], [1, -2.43, 9.39]], bonds: [[0, 1, 1], [0, 2, 1]] };
  const fe = Z('Fe'), o = Z('O');
  const sim = newSim(33);
  sim.addMoles(fe, 0.2, P.W / 2, P.H - 150); sim.addMoles(o, 0.2, P.W / 2, P.H * 0.4);
  sim.addMolecules(H2O, 200, P.W / 2, P.H - 400, { x: 300, y: 100 });
  sim.plateTemp = 1000;
  run(sim, 15);
  const hydride = speciesWith(sim, function (k) { return k[fe] && k[1] && !k[o]; });
  assert.ok(hydride.length === 0, 'no iron hydrides: ' + hydride.map(function (s) { return s.key; }).join(' '));
  for (let i = 0; i < sim.n; i++) {
    if (sim.el[i] !== fe) continue;
    for (let k = 0; k < sim.nb[i]; k++) assert.ok(sim.el[sim.bp[i * P.MAXB + k]] !== 1, 'no Fe-H bond');
  }
  assert.ok(compound(sim, fe, o).reacted > 0.02, 'iron oxidises in hot steam: ' + JSON.stringify(compound(sim, fe, o)));
  assert.ok(mol(sim, 'H2') > 0.002, 'hydrogen is released: ' + mol(sim, 'H2').toFixed(3) + ' mol');
});

test('magnesium and aluminium burn to MgO and Al2O3', function () {
  // [metal, textbook O:M, allowed spread]: a thin scale on the metal also holds some interface
  // atoms that are only partly oxidised, which lowers the average a little
  [['Mg', 1.0, 0.3], ['Al', 1.5, 0.45]].forEach(function (c, n) {
    const m = Z(c[0]), o = Z('O');
    const sim = newSim(34 + n);
    sim.addMoles(m, 0.2, P.W / 2, P.H - 150); sim.addMoles(o, 0.2, P.W / 2, P.H * 0.4);
    sim.plateTemp = 1300;
    run(sim, 12);
    const r = compound(sim, m, o);
    assert.ok(r.reacted > 0.08, c[0] + ' burns: ' + JSON.stringify(r));
    assert.ok(r.ratio > c[1] - c[2] && r.ratio < c[1] + 0.3, c[0] + ' oxide O:' + c[0] + ' ' + r.ratio.toFixed(2) + ' (expected ' + c[1] + ')');
  });
});

test('sodium and chlorine make sodium chloride, and salt is a solid', function () {
  const na = Z('Na'), cl = Z('Cl');
  const sim = newSim(36);
  sim.addMoles(na, 0.2, P.W / 2, P.H - 150); sim.addMoles(cl, 0.2, P.W / 2, P.H * 0.4);
  run(sim, 25); // the reaction runs at the surface of the sodium lump
  const r = compound(sim, na, cl);
  assert.ok(r.reacted > 0.05 && Math.abs(r.ratio - 1) < 0.1, 'NaCl forms 1:1: ' + JSON.stringify(r));
  // loose salt at room temperature clumps into a solid
  const NaCl = { atoms: [[na, 0, 0], [cl, 29.5, 0]], bonds: [[0, 1, 1]] };
  const s2 = newSim(37);
  s2.addMolecules(NaCl, 150, P.W / 2, P.H / 2, { x: 600, y: 700 }, 293, 2);
  run(s2, 15);
  const salt = s2.analyze().find(function (x) { return x.key === 'ClNa'; });
  assert.ok(salt && salt.condensed / salt.molecules > 0.5, 'salt clumps into a solid: condensed ' + (salt ? (salt.condensed / salt.molecules).toFixed(2) : 'none'));
});

test('hydrogen and chlorine stay apart until sparked, then make HCl; N2 + H2 needs a catalyst', function () {
  const sim = newSim(37);
  sim.addMoles(Z('H'), 0.2, P.W / 2, P.H * 0.4); sim.addMoles(Z('Cl'), 0.2, P.W / 2, P.H * 0.6);
  run(sim, 2);
  assert.ok(mol(sim, 'ClH') < 0.01, 'dark H2/Cl2 mixture is stable');
  sim.startBurst(P.W / 2, P.H / 2, 550, 4000, 0.16); // the spark tool: heat only
  run(sim, 12);
  assert.ok(mol(sim, 'ClH') > 0.2, 'HCl after the spark: ' + mol(sim, 'ClH').toFixed(3) + ' mol');
  const haber = newSim(38);
  haber.addMoles(Z('N'), 0.2, P.W / 2, P.H * 0.4); haber.addMoles(Z('H'), 0.6, P.W / 2, P.H * 0.6);
  run(haber, 6);
  assert.ok(mol(haber, 'H3N') < 0.005, 'no ammonia without a catalyst');
});

// ---- solids, temperature and bookkeeping ------------------------------------------------
// a lump of `z` (or of the given template) settled on the floor, lifted by a patch at its top and
// thrown at the left wall: returns the share of the lump that rose with the hand and the share in
// its largest piece afterwards
function liftAndThrow(seed, setup) {
  const sim = newSim(seed);
  setup(sim);
  run(sim, 6);
  const lump = [];
  for (let i = 0; i < sim.n; i++) if (sim.condF[i] > 0.5 || sim.isSolidAtom(i)) lump.push(i);
  let top = P.H, cx = 0;
  lump.forEach(function (i) { top = Math.min(top, sim.y[i]); cx += sim.x[i] / lump.length; });
  const y0 = lump.map(function (i) { return sim.y[i]; });
  sim.grabStart(cx, top + 20, 120);
  let steps = Math.round(3 / P.DT);
  for (let s = 0; s < steps; s += 6) { sim.grabMove(cx, top + 20 - 700 * Math.min(1, s / (0.67 * steps))); sim.step(6); }
  let rose = 0;
  lump.forEach(function (i, k) { if (sim.y[i] < y0[k] - 500) rose++; });
  steps = Math.round(0.4 / P.DT);
  for (let s = 0; s < steps; s += 6) { sim.grabMove(cx - (cx - 300) * s / steps, top - 680); sim.step(6); }
  sim.grabEnd();
  run(sim, 3);
  // largest piece: atoms within 1.4 contact distances (or bonded) of each other
  const parent = new Int32Array(sim.n);
  for (let i = 0; i < sim.n; i++) parent[i] = i;
  const find = function (a) { while (parent[a] !== a) { parent[a] = parent[parent[a]]; a = parent[a]; } return a; };
  lump.forEach(function (i) {
    lump.forEach(function (j) {
      if (j <= i) return;
      const d = Math.hypot(sim.x[j] - sim.x[i], sim.y[j] - sim.y[i]);
      if (d < 1.4 * (CE.elements[sim.el[i]].radius + CE.elements[sim.el[j]].radius)) parent[find(i)] = find(j);
    });
  });
  const size = new Map();
  lump.forEach(function (i) { const r = find(i); size.set(r, (size.get(r) || 0) + 1); });
  return { rose: rose / lump.length, piece: Math.max.apply(null, Array.from(size.values())) / lump.length };
}

test('lumps of lead and graphite come up whole when grabbed and survive a throw', function () {
  const lead = liftAndThrow(40, function (sim) { sim.addMoles(Z('Pb'), 0.3, P.W / 2, P.H - 150); });
  assert.ok(lead.rose > 0.85 && lead.piece > 0.85, 'lead: ' + JSON.stringify(lead));
  const graphite = liftAndThrow(41, function (sim) { sim.addMoles(Z('C'), 0.3, P.W / 2, P.H - 150); });
  assert.ok(graphite.rose > 0.85 && graphite.piece > 0.85, 'graphite: ' + JSON.stringify(graphite));
});

test('holding and moving a lump does not warm it, and iron thrown at the glass stays whole', function () {
  // internal temperature of a set of atoms: motion relative to their own centre of mass
  function Tint(sim, set) {
    let vx = 0, vy = 0;
    set.forEach(function (i) { vx += sim.x[i] - sim.px[i]; vy += sim.y[i] - sim.py[i]; });
    vx /= set.length; vy /= set.length;
    let ke = 0;
    set.forEach(function (i) {
      const ux = sim.x[i] - sim.px[i] - vx, uy = sim.y[i] - sim.py[i] - vy;
      ke += 0.5 * CE.elements[sim.el[i]].mass * (ux * ux + uy * uy) / (P.DT * P.DT);
    });
    return ke / (set.length * P.KB);
  }
  function avgT(sim, set, seconds, move) {
    const steps = Math.round(seconds / P.DT);
    let T = 0, c = 0;
    for (let s = 0; s < steps; s += 6) { if (move) move(s / steps); sim.step(6); if (s % 60 === 0) { T += Tint(sim, set); c++; } }
    return T / c;
  }
  ['Fe', 'Au'].forEach(function (sym, k) {
    const sim = newSim(50 + k);
    sim.addMoles(Z(sym), 0.3, P.W / 2, P.H - 150);
    run(sim, 6);
    const saved = [P.CONDUCTION, P.PLATE_GAMMA];
    sim.convection = 0; sim.insulated = true; P.CONDUCTION = 0; P.PLATE_GAMMA = 0;
    run(sim, 2); // (settling into the insulated jar)
    const all = Array.from({ length: sim.n }, function (_, i) { return i; });
    let top = P.H, cx = 0;
    all.forEach(function (i) { top = Math.min(top, sim.y[i]); cx += sim.x[i] / all.length; });
    const T0 = avgT(sim, all, 1);
    sim.grabStart(cx, top + 20, 120);
    // lift, carry sideways, hold still
    avgT(sim, all, 3, function (f) { sim.grabMove(cx - 600 * Math.max(0, f * 2 - 1), top + 20 - 800 * Math.min(1, f * 2)); });
    const T1 = avgT(sim, all, 2);
    P.CONDUCTION = saved[0]; P.PLATE_GAMMA = saved[1];
    let bottom = 0;
    all.forEach(function (i) { bottom = Math.max(bottom, sim.y[i]); });
    assert.ok(bottom < P.H - 400, sym + ' lump did not come up whole (bottom at ' + Math.round(P.H - bottom) + ')');
    assert.ok(T1 - T0 < 30, sym + ' warmed from ' + Math.round(T0) + ' K to ' + Math.round(T1) + ' K while held and carried');
    if (sym !== 'Fe') return;
    // fling it at the left wall and let go
    sim.insulated = false; sim.convection = P.CONVECTION;
    avgT(sim, all, 0.6, function (f) { sim.grabMove(cx - 600 - 3000 * f, top + 20 - 800); });
    sim.grabEnd();
    run(sim, 3);
    const parent = new Int32Array(sim.n);
    for (let i = 0; i < sim.n; i++) parent[i] = i;
    const find = function (a) { while (parent[a] !== a) { parent[a] = parent[parent[a]]; a = parent[a]; } return a; };
    const r = 2.8 * CE.elements[Z(sym)].radius;
    for (let i = 0; i < sim.n; i++) for (let j = i + 1; j < sim.n; j++) if (Math.hypot(sim.x[j] - sim.x[i], sim.y[j] - sim.y[i]) < r) parent[find(i)] = find(j);
    const size = new Map();
    all.forEach(function (i) { const q = find(i); size.set(q, (size.get(q) || 0) + 1); });
    const piece = Math.max.apply(null, Array.from(size.values())) / all.length;
    assert.ok(piece > 0.97, 'thrown iron broke up: largest piece ' + piece.toFixed(2));
  });
});

test('graphite on a hot plate reads the plate temperature', function () {
  const sim = newSim(42);
  sim.insulated = true;
  sim.addMoles(Z('C'), 0.3, P.W / 2, P.H - 200);
  sim.plateTemp = 800;
  run(sim, 12);
  let T = 0;
  for (let k = 0; k < 30; k++) { run(sim, 0.1); T += sim.temperature / 30; }
  assert.ok(T > 500 && T < 1100, 'graphite at ' + Math.round(T) + ' K on an 800 K plate');
});

test('freshly added molecules start near room temperature', function () {
  const CO2 = { atoms: [[6, 0, 0], [8, 12.4, 0], [8, -12.4, 0]], bonds: [[0, 1, 2], [0, 2, 2]] };
  const EtOH = { atoms: [[6, 0, 0], [6, 15.2, 0], [8, 15.2, -14.2], [1, 24.9, -14.2], [1, 0, -10.7], [1, 0, 10.7],
    [1, -10.7, 0], [1, 15.2, 10.7], [1, 25.9, 0]],
    bonds: [[0, 1, 1], [1, 2, 1], [2, 3, 1], [0, 4, 1], [0, 5, 1], [0, 6, 1], [1, 7, 1], [1, 8, 1]] };
  [[CO2, 500], [EtOH, 300]].forEach(function (c, k) {
    const sim = newSim(43 + k);
    sim.react = function () {}; sim.bondEvents = function () {};
    const added = sim.addMolecules(c[0], c[1], P.W / 2, P.H - 250, { x: 400, y: 250 });
    assert.strictEqual(added, c[1], 'all ' + c[1] + ' molecules placed');
    let hottest = 0;
    for (let s = 0; s < 60; s += 6) { sim.step(6); for (let i = 0; i < sim.n; i++) hottest = Math.max(hottest, sim.tLoc[i]); }
    assert.ok(hottest < 1500, 'no overlap flash: hottest spot ' + Math.round(hottest) + ' K');
  });
});

test('carbon monoxide and gold are stable; iron needs heat to burn in chlorine; H2 + F2 explodes in the dark', function () {
  const CO = { atoms: [[6, 0, 0], [8, 11.2, 0]], bonds: [[0, 1, 3]] };
  const co = newSim(45);
  co.addMolecules(CO, 300, P.W / 2, P.H / 2, { x: 900, y: 1200 });
  run(co, 6);
  assert.ok(mol(co, 'CO') > 0.29, 'CO does not polymerise: ' + mol(co, 'CO').toFixed(3) + ' mol left');
  const au = newSim(46);
  au.addMoles(Z('Au'), 0.3, P.W / 2, P.H - 150); au.addMoles(Z('O'), 0.3, P.W / 2, P.H * 0.4);
  au.plateTemp = 900;
  run(au, 8);
  assert.ok(compound(au, Z('Au'), Z('O')).reacted < 0.01, 'gold does not oxidise: ' + JSON.stringify(compound(au, Z('Au'), Z('O'))));
  const fe = newSim(47);
  fe.addMoles(Z('Fe'), 0.3, P.W / 2, P.H - 150); fe.addMoles(Z('Cl'), 0.4, P.W / 2, P.H * 0.4);
  run(fe, 6);
  assert.ok(compound(fe, Z('Fe'), Z('Cl')).reacted < 0.03, 'iron is untouched by chlorine at room temperature: ' + JSON.stringify(compound(fe, Z('Fe'), Z('Cl'))));
  fe.plateTemp = 700;
  run(fe, 10);
  assert.ok(compound(fe, Z('Fe'), Z('Cl')).reacted > 0.05, 'hot iron burns in chlorine: ' + JSON.stringify(compound(fe, Z('Fe'), Z('Cl'))));
  const hf = newSim(48);
  hf.addMoles(Z('H'), 0.3, P.W / 2, P.H * 0.4); hf.addMoles(Z('F'), 0.3, P.W / 2, P.H * 0.6);
  run(hf, 10); // (it goes off once the two layers have mixed: a random 1-6 s)
  assert.ok(mol(hf, 'FH') > 0.2, 'H2 + F2 reacts without a spark: ' + mol(hf, 'FH').toFixed(3) + ' mol HF');
});

// (pure peroxide vapour needs about 200 kJ/mol to split; below its boiling point it lasts, and the
// quick decomposition seen in practice is catalysed by surfaces and traces of metal)
test('hot hydrogen peroxide decomposes to water and oxygen, without O-O-O chains', function () {
  const H2O2 = { atoms: [[1, -3, -9.2], [8, 0, 0], [8, 14.2, 0], [1, 17.2, 9.2]], bonds: [[0, 1, 1], [1, 2, 1], [2, 3, 1]] };
  const sim = newSim(49);
  sim.addMolecules(H2O2, 300, P.W / 2, P.H - 250, { x: 400, y: 250 });
  sim.plateTemp = 750;
  run(sim, 15);
  assert.ok(mol(sim, 'H2O') > 0.1 && mol(sim, 'O2') > 0.03, 'water and oxygen: ' + mol(sim, 'H2O').toFixed(3) + ' / ' + mol(sim, 'O2').toFixed(3));
  // (HO3 is a real, short-lived radical, made from ozone; H2O3 and longer chains are not)
  const chains = speciesWith(sim, function (k) { return k[8] >= 3 && k[1] >= 2; });
  assert.ok(chains.reduce(function (a, s) { return a + s.mol; }, 0) < 0.005, 'no polyoxides: ' + chains.map(function (s) { return s.key; }).join(' '));
});

test('the gauge reads no vapour over a non-volatile solid, and a lump with a stray atom keeps its amount', function () {
  const sim = newSim(50);
  sim.plateTemp = 280;
  sim.addMoles(Z('Ga'), 0.5, P.W / 2, P.H - 150);
  run(sim, 6);
  assert.ok(sim.gasPressure() < 0.01, 'gallium vapour reads ' + sim.gasPressure().toFixed(3) + ' atm');
  const c = newSim(51);
  c.addMoles(Z('C'), 0.3, P.W / 2, P.H - 150);
  run(c, 2);
  // stick one oxygen atom onto an edge carbon of the graphite
  let edge = -1;
  for (let i = 0; i < c.n && edge < 0; i++) if (c.nb[i] >= 2 && c.canAccept(i, Z('O'))) edge = i;
  assert.ok(edge >= 0, 'graphite has an edge carbon with a free valence');
  const o = c.addAtom(Z('O'), c.x[edge], c.y[edge] - 14);
  c.addBond(edge, o, 1, C.bondClass(Z('C'), Z('O')));
  // (the graphite may have grown as more than one sheet: count every network)
  const nets = c.analyze().filter(function (s) { return s.key.slice(-2) === ')n'; });
  const total = nets.reduce(function (a, s) { return a + s.mol; }, 0);
  assert.ok(total > 0.25, 'the lump still counts about 0.3 mol: ' + nets.map(function (s) { return s.key + ' ' + s.mol.toFixed(3); }).join(', '));
});

test('sulfuric acid and benzene do not react with themselves at room temperature', function () {
  // (a concerted four-centre exchange between two whole molecules used to condense sulfuric acid
  // into S-O-S chains and polymerise benzene, which also lacked its aromatic stabilisation)
  const H2SO4 = { atoms: [[16, 0, 0], [8, 0, -14.9], [8, 0, 14.9], [8, -17.1, 0], [1, -26.8, 0], [8, 17.1, 0], [1, 26.8, 0]],
    bonds: [[0, 1, 2], [0, 2, 2], [0, 3, 1], [3, 4, 1], [0, 5, 1], [5, 6, 1]] };
  const C6H6 = { atoms: [[6, 14, 0], [6, 7, 12.1], [6, -7, 12.1], [6, -14, 0], [6, -7, -12.1], [6, 7, -12.1],
    [1, 24.8, 0], [1, 12.4, 21.5], [1, -12.4, 21.5], [1, -24.8, 0], [1, -12.4, -21.5], [1, 12.4, -21.5]],
    bonds: [[0, 1, 2], [1, 2, 1], [2, 3, 2], [3, 4, 1], [4, 5, 2], [5, 0, 1], [0, 6, 1], [1, 7, 1], [2, 8, 1], [3, 9, 1], [4, 10, 1], [5, 11, 1]] };
  [['H2O4S', H2SO4], ['C6H6', C6H6]].forEach(function (c, k) {
    const sim = newSim(60 + k);
    sim.addMolecules(c[1], 200, P.W / 2, P.H - 250, { x: 450, y: 220 });
    run(sim, 10);
    const left = sim.analyze().find(function (s) { return s.key === c[0]; });
    assert.ok(left && left.molecules > 190, c[0] + ' left: ' + (left ? left.molecules : 0) + ' of 200 (T ' + Math.round(sim.temperature) + ' K)');
  });
});

test('a hydrogen flame gives off light, violet-blue', function () {
  const sim = newSim(62);
  sim.addMoles(Z('H'), 1, P.W / 2, P.H / 2); sim.addMoles(Z('O'), 0.5, P.W / 2, P.H / 2);
  run(sim, 1);
  assert.ok(!sim.light || sim.light.n === 0, 'no light before the spark');
  sim.startBurst(P.W / 2, P.H / 2, 550, 4000, 0.16);
  run(sim, 1);
  assert.ok(sim.light && sim.light.n > 50, 'flashes from the flame: ' + (sim.light ? sim.light.n : 0));
  let r = 0, g = 0, b = 0;
  for (let k = 0; k < sim.light.n; k++) { const e = sim.light.buf[6 * k + 2]; r += e * sim.light.buf[6 * k + 3]; g += e * sim.light.buf[6 * k + 4]; b += e * sim.light.buf[6 * k + 5]; }
  assert.ok(b > r && b > g, 'violet-blue: ' + [r, g, b].map(function (v) { return Math.round(v / (r + g + b) * 100); }).join('/'));
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
  // TEST_FILTER=text runs only the tests whose name contains it
  const only = process.env.TEST_FILTER;
  const run = only ? tests.filter(function (t) { return t.name.indexOf(only) >= 0; }) : tests;
  run.forEach(function (t) {
    const s = Date.now();
    try {
      t.fn();
      console.log('  ok   ' + t.name + ' (' + (Date.now() - s) + ' ms)');
    } catch (e) {
      failed++;
      console.log('  FAIL ' + t.name + '\n       ' + (e && e.message));
    }
  });
  console.log((run.length - failed) + '/' + run.length + ' passed in ' + ((Date.now() - t0) / 1000).toFixed(1) + ' s');
  process.exit(failed ? 1 : 0);
})();
