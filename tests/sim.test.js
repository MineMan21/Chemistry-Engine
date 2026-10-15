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
