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

test('water freezes into an ordered lattice and is recognised as a solid', function () {
  const H2O = { atoms: [[8, 0, 0], [1, 9.7, 0], [1, -2.43, 9.39]], bonds: [[0, 1, 1], [0, 2, 1]] };
  const sim = newSim(12);
  sim.react = function () {}; sim.bondEvents = function () {};
  sim.addMolecules(H2O, 500, P.W / 2, P.H - 200, { x: 350, y: 150 });
  sim.plateTemp = 200;
  run(sim, 10);
  const O = [];
  for (let i = 0; i < sim.n; i++) if (sim.el[i] === 8) O.push(i);
  let solid = 0, order = 0, counted = 0;
  O.forEach(function (i) {
    if (sim.isSolidAtom(i)) solid++;
    // local hexagonal bond-orientational order of the oxygen lattice
    let re = 0, im = 0, k = 0;
    O.forEach(function (j) {
      if (j === i) return;
      const dx = sim.x[j] - sim.x[i], dy = sim.y[j] - sim.y[i];
      if (dx * dx + dy * dy < 36 * 36) { const a = 6 * Math.atan2(dy, dx); re += Math.cos(a); im += Math.sin(a); k++; }
    });
    if (k >= 2) { order += Math.hypot(re, im) / k; counted++; }
  });
  assert.ok(sim.temperature < 260, 'cooled to ' + sim.temperature);
  assert.ok(solid > 0.8 * O.length, 'ice counted as solid: ' + solid + ' of ' + O.length);
  assert.ok(order / counted > 0.7, 'lattice order psi6 ' + (order / counted).toFixed(2));
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
