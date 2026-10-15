// Particle simulation of the jar. No DOM access, so it also runs headless under Node.
//
// Units: 1 world unit = 10 pm, mass in u (amu), time in simulated seconds, temperature in K.
// Every particle is one atom on screen but statistically stands for MOL_PER_PARTICLE moles of
// that atom, which is how the jar can hold whole moles while staying real-time.
(function (CE) {
  'use strict';

  const E = CE.elements;
  const C = CE.chem;
  const NZ = C.NZ;

  const P = {
    MAX: 60000,
    MAXB: 8, // max bonds per atom
    KB: 1500, // simulation Boltzmann constant: m*v^2/2 = KB*T per particle (2D)
    R_KJ: 0.008314, // kJ/(mol K)
    DT: 0.0015,
    W: 3200, // jar interior width (world units)
    H: 4000, // jar interior height
    MOL_PER_PARTICLE: 0.001,
    JAR_VOLUME_L: 24, // real volume the jar represents (1 mol of gas at 293 K ~ 1 atm)
    AMBIENT: 293,
    K_REP: 160000, // repulsion stiffness (per unit reduced mass)
    RC_FACTOR: 1.45, // attraction cutoff as a multiple of contact distance
    KC: 2.0e8, // Coulomb constant (sim energy * units / e^2)
    PARTIAL_CHARGE: 0.3, // partial charge per unit electronegativity difference
    VDW_SCALE: 2.25, // calibration factor for van der Waals well depths (2D, short-ranged well)
    MOL_POLARIZABILITY: 0, // (disabled) per-bond dispersion boost; it made bonding exothermic
    HB_EPS: 550, // hydrogen bond well depth (K) for O-H...O; F stronger, N weaker
    HB_CONTACT: 0.78, // hydrogen-bonded H...X pairs may approach to this fraction of contact
    METAL_EPS_PER_BP: 1.05, // metallic cohesion well depth (K) per K of boiling point
    G_GAS: 1, // gravity felt by isolated (gas) particles
    G_CONDENSED: 150, // extra gravity for particles touching others (liquids/solids settle)
    BULK_VMAX: 50, // max bulk (flow) speed of condensed matter, units/s -> gentle settling
    BULK_VMAX_FLOOR: 15, // ... near the floor, so material lands softly
    LN_COLL: 18, // ln of real collisions represented by one simulated contact-step
    LN_UNI: 16, // ln of the unimolecular attempt frequency (bond dissociation)
    ALPHA: 0.15, // fraction of broken-bond energy needed as activation (Evans-Polanyi like)
    ALPHA_4C: 0.22, // same for concerted four-centre exchanges (mostly symmetry forbidden)
    SPIN_BARRIER: 60, // extra activation (kJ/mol) when closed-shell molecules or metal surfaces
    // attack triplet O=O (why O2 is slow with fuels and metals but fast with radicals)
    HEAT_FRACTION: 0.6, // share of reaction enthalpy deposited as translational kinetic energy
    K_WALL: 110000, // soft wall stiffness (per unit mass)
    K_ANGLE: 8000, // bending spring stiffness (per unit reduced mass); bending is soft
    ACCOUNT_ANGLES: 0, // include bending energy in per-event energy accounting
    PLATE_GAMMA: 80, // Langevin coupling (1/s) of atoms near the hot plate
    PLATE_ZONE: 60, // thickness of the layer above the plate that exchanges heat with it
    WALL_GAMMA: 3, // Langevin coupling to the room through the glass
    WALL_ZONE: 30,
    CONVECTION: 0.4, // convective heat exchange (1/s) of the whole jar with the plate, max at bottom
    VALVE_HALF: 140, // half-width of the valve opening in the lid
    MAX_STEP: 6, // max displacement per substep (keeps the integrator stable)
    TCELL: 160, // size of the coarse cells used for the local temperature field
    T_PRIOR: 2,
    SELF_T_RATE: 0.02, // EMA rate of an atom's own kinetic temperature (per substep)
    SELF_T_WEIGHT: 0.5, // weight of that own temperature in the local temperature // pseudo-count pulling sparse cells toward the jar temperature
    ANNEAL_STEPS: 1500, // newly added atoms are held at room temperature while they settle in
    RADICAL_CAPTURE: 0.6, // fraction of the attraction range within which radicals can react
    FRESH: 30, // substeps of softened repulsion after an atom is created or changes bonds
    BOND_RAMP: 0.5, // units per substep a new/changed bond moves toward its rest length
    SHAKE_ITERS: 3, // constraint solver iterations per substep
    MIN_DYN_MASS: 3 // minimum dynamical mass (u), see below
  };
  const ENERGY_PER_KJ = P.KB / P.R_KJ; // sim energy units per kJ/mol

  const GS = 6; // ghost slots per atom

  // Hydrogen-bond acceptor strength (N, O, F lone pairs)
  const HB_ACCEPTOR = new Float32Array(119);
  HB_ACCEPTOR[7] = 0.7; HB_ACCEPTOR[8] = 1.0; HB_ACCEPTOR[9] = 1.3;

  // Per-element constants
  const mass = new Float32Array(NZ), invMass = new Float32Array(NZ), radius = new Float32Array(NZ);
  const en = new Float32Array(NZ), maxVal = new Uint8Array(NZ), normVal = new Uint8Array(NZ);
  const isMetal = new Uint8Array(NZ), vdwEps = new Float32Array(NZ), metalEps = new Float32Array(NZ);
  const elColor = new Uint32Array(NZ);
  for (let z = 1; z < NZ; z++) {
    const e = E[z];
    // Dynamical mass. Hydrogen is made heavier (a standard molecular-dynamics trick) so its
    // fast vibrations stay stable at our time step; chemistry and displayed masses are real.
    mass[z] = Math.max(e.mass, P.MIN_DYN_MASS);
    invMass[z] = 1 / mass[z];
    radius[z] = e.radius;
    en[z] = e.en;
    maxVal[z] = e.maxValence;
    normVal[z] = e.normalValence;
    isMetal[z] = e.isMetal ? 1 : 0;
    vdwEps[z] = e.vdwEps;
    metalEps[z] = e.isMetal ? P.METAL_EPS_PER_BP * (e.bp || (e.mp ? e.mp * 2.5 : 3000)) : 0;
    elColor[z] = packColor(e.color);
  }
  // Pair tables
  const pairSigma = new Float32Array(NZ * NZ), pairRc = new Float32Array(NZ * NZ);
  const pairEps = new Float32Array(NZ * NZ), pairMetalEps = new Float32Array(NZ * NZ);
  for (let a = 1; a < NZ; a++) {
    for (let b = 1; b < NZ; b++) {
      const p = a * NZ + b;
      pairSigma[p] = radius[a] + radius[b];
      pairRc[p] = pairSigma[p] * P.RC_FACTOR;
      pairEps[p] = Math.sqrt(vdwEps[a] * vdwEps[b]) * P.KB;
      pairMetalEps[p] = Math.sqrt(metalEps[a] * metalEps[b]) * P.KB;
    }
  }

  function packColor(hex) {
    const r = parseInt(hex.substr(1, 2), 16), g = parseInt(hex.substr(3, 2), 16), b = parseInt(hex.substr(5, 2), 16);
    return (255 << 24) | (b << 16) | (g << 8) | r; // little-endian ABGR for ImageData
  }

  function gaussian(rng) {
    let u = 0, v = 0;
    while (u === 0) u = rng();
    v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  function Simulation(opts) {
    opts = opts || {};
    this.rng = opts.rng || Math.random;
    const M = P.MAX, MB = P.MAXB;
    this.n = 0;
    this.x = new Float64Array(M); this.y = new Float64Array(M);
    this.px = new Float64Array(M); this.py = new Float64Array(M);
    this.fx = new Float32Array(M); this.fy = new Float32Array(M);
    this.el = new Uint8Array(M);
    this.nb = new Uint8Array(M); // bond count
    this.bp = new Int32Array(M * MB).fill(-1); // bond partner
    this.bo = new Uint8Array(M * MB); // bond order
    this.bc = new Uint8Array(M * MB); // bond class (1 covalent, 2 ionic)
    this.bl = new Float32Array(M * MB); // current target length (ramps to the rest length)
    this.fresh = new Uint8Array(M); // grace counter after creation / bond change
    this.eRef = new Float32Array(M); // kinetic energy allowance while fresh (caps speed)
    // "ghost" exclusions: atoms that were just bonded (or 1-3 neighbours) ignore each other's
    // repulsion until they drift apart, so breaking a bond does not fire the fragments apart
    this.ghost = new Int32Array(M * GS).fill(-1);
    this.ghostT = new Uint8Array(M); // number of active ghost slots
    this.used = new Uint8Array(M); // total bond order (valence used)
    this.q = new Float32Array(M); // charge (ionic) / partial charge (covalent)
    this.tLoc = new Float32Array(M); // local temperature used for kinetics
    this.tSelf = new Float32Array(M); // each atom's own recent kinetic temperature (EMA)
    this.contacts = new Uint8Array(M); // non-bonded neighbours in contact
    this.molId = new Int32Array(M);
    this.hbd = new Uint8Array(M); // hydrogen-bond donor (H bonded to N, O or F)
    // batch tag of atoms added together as an element: bonds between same-element atoms of one
    // batch just assemble the element's standard state (S8, graphite...), which releases no heat
    this.tag = new Int32Array(M);
    this.anneal = new Uint16Array(M); // substeps of strong thermostatting left after being added
    this.batch = 0;
    this.molColor = new Uint32Array(M);
    // spatial grid
    this.cellSize = 40;
    this.gridW = 0; this.gridH = 0;
    this.cellStart = new Int32Array(1);
    this.cellItems = new Int32Array(M);
    this.cellOf = new Int32Array(M);
    this.cellOfT = new Int32Array(M);
    this.present = new Uint8Array(NZ);
    // coarse temperature field (local temperature used for reaction kinetics)
    this.tgw = Math.ceil(P.W / P.TCELL); this.tgh = Math.ceil(P.H / P.TCELL);
    this.tCell = new Float32Array(this.tgw * this.tgh).fill(P.AMBIENT);
    this.tSum = new Float64Array(this.tgw * this.tgh);
    this.tCnt = new Float32Array(this.tgw * this.tgh);
    this.cellKE = new Float64Array(this.tgw * this.tgh); // kinetic energy per coarse cell
    this.cellMom = new Float64Array(this.tgw * this.tgh * 3); // condensed matter momentum, mass
    this.cellDebt = new Float64Array(this.tgw * this.tgh); // energy owed by endothermic events
    // reaction candidates
    this.candI = new Int32Array(M * 2);
    this.candJ = new Int32Array(M * 2);
    this.nCand = 0;
    // controls
    this.plateTemp = P.AMBIENT;
    this.convection = P.CONVECTION;
    this.insulated = false;
    this.vacuumOn = false;
    this.canOn = false;
    this.canGas = 'N2';
    this.canRate = 0.5; // mol per simulated second
    this.canAccumulator = 0;
    this.time = 0;
    this.substep = 0;
    // statistics
    this.impulse = 0; this.impulseTime = 0;
    this.pressureAtm = 0;
    this.temperature = P.AMBIENT;
    this.reactionCount = 0;
    this.reactionsPerSec = 0;
    this.reactionWindow = 0; this.reactionWindowTime = 0;
    this.removedMol = 0;
    this.lastReactions = [];
    this.updateGrid();
  }

  const S = Simulation.prototype;
  Simulation.params = P;
  Simulation.ENERGY_PER_KJ = ENERGY_PER_KJ;
  Simulation.elColor = elColor;
  Simulation.packColor = packColor;

  // ---------------------------------------------------------------------------------------
  // Particle management
  // ---------------------------------------------------------------------------------------
  S.addAtom = function (z, x, y, T, tag) {
    if (this.n >= P.MAX) return -1;
    const i = this.n++;
    const sd = Math.sqrt(P.KB * (T || P.AMBIENT) / mass[z]) * P.DT;
    this.x[i] = x; this.y[i] = y;
    this.px[i] = x - gaussian(this.rng) * sd;
    this.py[i] = y - gaussian(this.rng) * sd;
    this.el[i] = z;
    this.nb[i] = 0; this.used[i] = 0; this.q[i] = 0; this.hbd[i] = 0;
    this.tag[i] = tag || 0;
    this.anneal[i] = P.ANNEAL_STEPS;
    for (let k = 0; k < P.MAXB; k++) this.bp[i * P.MAXB + k] = -1;
    this.tLoc[i] = T || P.AMBIENT;
    this.tSelf[i] = T || P.AMBIENT;
    this.contacts[i] = 0;
    this.fresh[i] = 0;
    this.eRef[i] = P.KB * (T || P.AMBIENT);
    this.ghostT[i] = 0;
    for (let k = 0; k < GS; k++) this.ghost[i * GS + k] = -1;
    this.molColor[i] = elColor[z];
    this.molId[i] = i;
    if (!this.present[z]) { this.present[z] = 1; this.updateGrid(); }
    return i;
  };

  // Adds `count` formula units of an element in its natural form (diatomic or single atoms)
  // spread around (cx, cy). Returns number of atoms added.
  S.addSubstance = function (z, count, cx, cy, opts) {
    opts = opts || {};
    const diatomic = opts.natural !== false && CE.DIATOMIC.has(z);
    const T = opts.T || P.AMBIENT;
    const r = radius[z];
    // gases are released spread out; liquids/solids as a compact blob
    const gas = E[z].bp !== null && E[z].bp < T;
    const spacing = opts.spacing || r * (gas ? 9 : diatomic ? 3.2 : 2.05);
    const added0 = this.n;
    const tag = ++this.batch;
    let placed = 0;
    // occupancy hash of existing + new atoms to avoid overlapping spawns
    const cell = 2 * r + 30, occ = new Map();
    const keyOf = function (cx, cy) { return cx * 100003 + cy; };
    const occupy = function (x, y) {
      const k = keyOf(Math.floor(x / cell), Math.floor(y / cell));
      let a = occ.get(k); if (!a) { a = []; occ.set(k, a); } a.push(x, y);
    };
    for (let i = 0; i < this.n; i++) occupy(this.x[i], this.y[i]);
    const minD = diatomic ? 2 * r + 8 : 2 * r * 0.95;
    const free = function (x, y) {
      const cx = Math.floor(x / cell), cy = Math.floor(y / cell);
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) {
        const arr = occ.get(keyOf(cx + a, cy + b));
        if (!arr) continue;
        for (let m = 0; m < arr.length; m += 2) {
          const dx = arr[m] - x, dy = arr[m + 1] - y;
          if (dx * dx + dy * dy < minD * minD) return false;
        }
      }
      return true;
    };
    // Gases: jittered spiral (sunflower) around the drop point. Solids/liquids: a close-packed
    // triangular lattice, so the material starts near its relaxed density (no heat of compaction).
    // Whatever does not fit is scattered randomly through the jar.
    const golden = Math.PI * (3 - Math.sqrt(5));
    let spiral = 0;
    let lattice = null;
    if (!gas && !diatomic) {
      lattice = [];
      const a = spacing, rows = Math.ceil(Math.sqrt(count * 1.3)) + 2;
      for (let j = -rows; j <= rows; j++) {
        for (let i2 = -rows; i2 <= rows; i2++) {
          const lx = (i2 + (j & 1) * 0.5) * a, ly = j * a * 0.866;
          lattice.push([lx, ly, lx * lx + ly * ly]);
        }
      }
      lattice.sort(function (p, q) { return p[2] - q[2]; });
    }
    for (let k = 0; k < count; k++) {
      let x = 0, y = 0, ok = false;
      for (let attempt = 0; attempt < 40 && !ok; attempt++) {
        if (attempt < 25) {
          if (lattice) {
            if (spiral >= lattice.length) { attempt = 24; continue; }
            x = cx + lattice[spiral][0]; y = cy + lattice[spiral][1];
          } else {
            const rr = spacing * Math.sqrt(spiral + 0.5) * 0.56, th = spiral * golden;
            x = cx + rr * Math.cos(th); y = cy + rr * Math.sin(th);
          }
          spiral++;
          if (x < r * 2 || x > P.W - r * 2 || y < r * 2 || y > P.H - r * 2) continue;
        } else {
          x = r * 2 + this.rng() * (P.W - r * 4);
          y = r * 2 + this.rng() * (P.H - r * 4);
        }
        ok = free(x, y);
      }
      if (!ok) continue;
      occupy(x, y);
      if (diatomic) {
        const L = C.pairLength[(z * NZ + z) * 4 + C.pairMaxOrder[z * NZ + z]] || 8;
        const a = this.rng() * Math.PI * 2;
        const i = this.addAtom(z, x - Math.cos(a) * L / 2, y - Math.sin(a) * L / 2, T, tag);
        if (i < 0) break;
        const j = this.addAtom(z, x + Math.cos(a) * L / 2, y + Math.sin(a) * L / 2, T, tag);
        if (j < 0) { this.removeParticle(i); break; }
        this.addBond(i, j, C.pairMaxOrder[z * NZ + z], C.pairClass[z * NZ + z]);
      } else {
        if (this.addAtom(z, x, y, T, tag) < 0) break;
      }
      placed++;
    }
    return { atoms: this.n - added0, units: placed, requested: count };
  };

  // Adds `count` copies of a small molecule template { atoms: [[z, dx, dy], ...],
  // bonds: [[a, b, order], ...] } scattered around (cx, cy) inside a box of half-size `spread`.
  // Used by the experiment presets (the user otherwise only adds elements).
  S.addMolecules = function (tpl, count, cx, cy, spread, T) {
    T = T || P.AMBIENT;
    let size = 0;
    tpl.atoms.forEach(function (a) { size = Math.max(size, Math.hypot(a[1], a[2]) + radius[a[0]] * 0.5); });
    // triangular lattice filling the box (liquid-like packing), visited in random order
    const a0 = size * 2.05, sites = [];
    for (let j = 0; ; j++) {
      const ly = cy - spread.y + j * a0 * 0.866;
      if (ly > cy + spread.y) break;
      for (let i = 0; ; i++) {
        const lx = cx - spread.x + (i + (j & 1) * 0.5) * a0;
        if (lx > cx + spread.x) break;
        if (lx > size && lx < P.W - size && ly > size && ly < P.H - size) sites.push(lx, ly);
      }
    }
    const nSites = sites.length / 2;
    for (let k = nSites - 1; k > 0; k--) { // shuffle
      const m = Math.floor(this.rng() * (k + 1));
      const tx = sites[2 * k], ty = sites[2 * k + 1];
      sites[2 * k] = sites[2 * m]; sites[2 * k + 1] = sites[2 * m + 1];
      sites[2 * m] = tx; sites[2 * m + 1] = ty;
    }
    let added = 0;
    const tag = ++this.batch;
    for (let k = 0; k < Math.min(count, nSites); k++) {
      const x = sites[2 * k], y = sites[2 * k + 1];
      const ang = this.rng() * Math.PI * 2, ca = Math.cos(ang), sa = Math.sin(ang);
      const ids = [];
      for (let m = 0; m < tpl.atoms.length; m++) {
        const t = tpl.atoms[m];
        const id = this.addAtom(t[0], x + t[1] * ca - t[2] * sa, y + t[1] * sa + t[2] * ca, T, tag);
        if (id < 0) return added;
        ids.push(id);
      }
      for (let m = 0; m < tpl.bonds.length; m++) {
        const b = tpl.bonds[m];
        this.addBond(ids[b[0]], ids[b[1]], b[2], C.pairClass[tpl.atoms[b[0]][0] * NZ + tpl.atoms[b[1]][0]]);
      }
      added++;
    }
    return added;
  };

  // Moles -> particles. Diatomic elements add moles of X2 molecules.
  S.addMoles = function (z, mol, cx, cy, opts) {
    const units = Math.max(1, Math.round(mol / P.MOL_PER_PARTICLE));
    return this.addSubstance(z, units, cx, cy, opts);
  };

  // Adds a gas from a formula spec (used by the pressure can).
  S.addGasMolecule = function (spec, x, y, T) {
    if (spec === 'air') spec = this.rng() < 0.78 ? 'N2' : this.rng() < 0.95 ? 'O2' : 'Ar';
    const z = CE.bySymbol[spec.replace(/2$/, '')].z;
    return this.addSubstance(z, 1, x, y, { T: T, natural: /2$/.test(spec) });
  };

  S.clear = function () {
    this.n = 0;
    this.tCell.fill(this.plateTemp);
    this.cellDebt.fill(0);
    this.present.fill(0);
    this.updateGrid();
  };

  // Removes particle i (breaking its bonds) by swapping the last particle into its slot.
  S.removeParticle = function (i) {
    const MB = P.MAXB;
    this._removing = true;
    while (this.nb[i] > 0) {
      const j = this.bp[i * MB];
      this.setBond(i, j, 0);
    }
    this._removing = false;
    // drop ghost relations of the removed atom
    for (let k = 0; k < GS; k++) {
      const g = this.ghost[i * GS + k];
      if (g >= 0) this.clearGhost(i, g);
    }
    const last = --this.n;
    if (i === last) return;
    this.x[i] = this.x[last]; this.y[i] = this.y[last];
    this.px[i] = this.px[last]; this.py[i] = this.py[last];
    this.el[i] = this.el[last]; this.nb[i] = this.nb[last]; this.used[i] = this.used[last];
    this.q[i] = this.q[last]; this.tLoc[i] = this.tLoc[last]; this.tSelf[i] = this.tSelf[last]; this.contacts[i] = this.contacts[last];
    this.hbd[i] = this.hbd[last]; this.tag[i] = this.tag[last]; this.anneal[i] = this.anneal[last];
    this.fresh[i] = this.fresh[last]; this.eRef[i] = this.eRef[last];
    for (let k = 0; k < GS; k++) {
      const g = this.ghost[last * GS + k];
      this.ghost[i * GS + k] = g;
      this.ghost[last * GS + k] = -1;
      if (g >= 0) for (let m = 0; m < GS; m++) if (this.ghost[g * GS + m] === last) this.ghost[g * GS + m] = i;
    }
    this.ghostT[i] = this.ghostT[last];
    this.ghostT[last] = 0;
    this.molId[i] = this.molId[last]; this.molColor[i] = this.molColor[last];
    for (let k = 0; k < MB; k++) {
      const j = this.bp[last * MB + k];
      this.bp[i * MB + k] = j;
      this.bo[i * MB + k] = this.bo[last * MB + k];
      this.bc[i * MB + k] = this.bc[last * MB + k];
      this.bl[i * MB + k] = this.bl[last * MB + k];
      this.bp[last * MB + k] = -1;
      if (j >= 0) {
        for (let m = 0; m < MB; m++) if (this.bp[j * MB + m] === last) this.bp[j * MB + m] = i;
      }
    }
  };

  // Collects all atoms bonded (directly or indirectly) to i.
  S.moleculeOf = function (i, out) {
    out = out || [];
    out.length = 0;
    const seen = this._seen || (this._seen = new Uint32Array(P.MAX));
    this._seenStamp = (this._seenStamp || 0) + 1;
    const stamp = this._seenStamp;
    out.push(i); seen[i] = stamp;
    for (let h = 0; h < out.length; h++) {
      const a = out[h];
      for (let k = 0; k < this.nb[a]; k++) {
        const b = this.bp[a * P.MAXB + k];
        if (seen[b] !== stamp) { seen[b] = stamp; out.push(b); }
      }
    }
    return out;
  };

  S.removeMolecule = function (i) {
    const atoms = this.moleculeOf(i, []);
    atoms.sort(function (a, b) { return b - a; });
    for (let k = 0; k < atoms.length; k++) this.removeParticle(atoms[k]);
    return atoms.length;
  };

  // ---------------------------------------------------------------------------------------
  // Bonds
  // ---------------------------------------------------------------------------------------
  S.bondSlot = function (i, j) {
    const MB = P.MAXB, base = i * MB;
    for (let k = 0; k < this.nb[i]; k++) if (this.bp[base + k] === j) return k;
    return -1;
  };

  S.bondOrder = function (i, j) {
    const k = this.bondSlot(i, j);
    return k < 0 ? 0 : this.bo[i * P.MAXB + k];
  };

  S.addBond = function (i, j, order, cls) {
    if (this.nb[i] >= P.MAXB || this.nb[j] >= P.MAXB) return false;
    const MB = P.MAXB;
    const dx = this.x[j] - this.x[i], dy = this.y[j] - this.y[i];
    const rest = C.pairLength[(this.el[i] * NZ + this.el[j]) * 4 + order];
    const len = Math.min(Math.max(Math.sqrt(dx * dx + dy * dy), rest), rest * 2.5);
    let k = this.nb[i]++;
    this.bp[i * MB + k] = j; this.bo[i * MB + k] = order; this.bc[i * MB + k] = cls; this.bl[i * MB + k] = len;
    k = this.nb[j]++;
    this.bp[j * MB + k] = i; this.bo[j * MB + k] = order; this.bc[j * MB + k] = cls; this.bl[j * MB + k] = len;
    this.markFresh(i); this.markFresh(j);
    this.used[i] += order; this.used[j] += order;
    this.updateCharge(i); this.updateCharge(j);
    return true;
  };

  // Sets bond order between i and j (0 removes the bond). Bond must exist.
  S.setBond = function (i, j, order) {
    const MB = P.MAXB;
    const ki = this.bondSlot(i, j), kj = this.bondSlot(j, i);
    if (ki < 0 || kj < 0) return;
    const old = this.bo[i * MB + ki];
    this.used[i] += order - old; this.used[j] += order - old;
    if (order > 0) {
      this.bo[i * MB + ki] = order; this.bo[j * MB + kj] = order;
    } else {
      this.removeSlot(i, ki); this.removeSlot(j, kj);
      if (!this._removing) {
        this.addGhost(i, j);
        for (let k = 0; k < this.nb[j]; k++) this.addGhost(i, this.bp[j * MB + k]);
        for (let k = 0; k < this.nb[i]; k++) this.addGhost(j, this.bp[i * MB + k]);
      }
    }
    this.updateCharge(i); this.updateCharge(j);
    this.markFresh(i); this.markFresh(j);
  };

  // Starts the grace period of atom i. Its kinetic energy may not exceed what it had (or the
  // local thermal energy) plus any reaction heat it receives, so geometry snapping after a bond
  // change cannot create energy.
  // Ghost exclusions are symmetric and last until the two atoms have separated to contact distance.
  S.addGhost = function (a, b) {
    if (a === b || this.isGhost(a, b)) return;
    if (!this._ghostSlot(a, b)) return;
    if (!this._ghostSlot(b, a)) { this._unslot(a, b); }
  };
  S._ghostSlot = function (a, b) {
    const g = this.ghost, base = a * GS;
    for (let k = 0; k < GS; k++) if (g[base + k] < 0) { g[base + k] = b; this.ghostT[a]++; return true; }
    return false; // full: the pair simply interacts normally
  };
  S._unslot = function (a, b) {
    const g = this.ghost, base = a * GS;
    for (let k = 0; k < GS; k++) if (g[base + k] === b) { g[base + k] = -1; this.ghostT[a]--; }
  };
  S.isGhost = function (i, j) {
    if (!this.ghostT[i]) return false;
    const g = this.ghost, bi = i * GS;
    for (let k = 0; k < GS; k++) if (g[bi + k] === j) return true;
    return false;
  };
  S.clearGhost = function (i, j) {
    this._unslot(i, j); this._unslot(j, i);
  };

  S.markFresh = function (i) {
    if (!this.fresh[i]) {
      const vx = (this.x[i] - this.px[i]) / P.DT, vy = (this.y[i] - this.py[i]) / P.DT;
      const ke = 0.5 * mass[this.el[i]] * (vx * vx + vy * vy);
      this.eRef[i] = Math.max(ke, 2 * P.KB * Math.max(this.tLoc[i], this.plateTemp));
    }
    this.fresh[i] = P.FRESH;
  };

  S.removeSlot = function (i, k) {
    const MB = P.MAXB, base = i * MB, last = --this.nb[i];
    this.bp[base + k] = this.bp[base + last];
    this.bo[base + k] = this.bo[base + last];
    this.bc[base + k] = this.bc[base + last];
    this.bl[base + k] = this.bl[base + last];
    this.bp[base + last] = -1;
  };

  S.updateCharge = function (i) {
    const MB = P.MAXB, base = i * MB, ei = this.el[i];
    let q = 0;
    for (let k = 0; k < this.nb[i]; k++) {
      const ej = this.el[this.bp[base + k]];
      if (this.bc[base + k] === 2) q += (en[ei] < en[ej] ? 1 : -1) * this.bo[base + k];
      else q += P.PARTIAL_CHARGE * (en[ej] - en[ei]);
    }
    this.q[i] = Math.max(-3, Math.min(3, q));
    let donor = 0;
    if (ei === 1) for (let k = 0; k < this.nb[i]; k++) if (HB_ACCEPTOR[this.el[this.bp[base + k]]]) donor = 1;
    this.hbd[i] = donor;
  };

  // Can atom i (element ei) accept one more bond unit to an atom of element ez?
  S.canAccept = function (i, ez) {
    const ei = this.el[i], u = this.used[i] + 1;
    if (u > maxVal[ei] || this.nb[i] >= P.MAXB) return false;
    if (u > normVal[ei] && !(en[ez] > en[ei] + 0.3)) return false; // hypervalent needs EN partner
    return true;
  };

  S.sameBatch = function (i, j) {
    return this.tag[i] !== 0 && this.tag[i] === this.tag[j] && this.el[i] === this.el[j];
  };

  S.isBonded = function (i, j) {
    const base = i * P.MAXB;
    for (let k = 0; k < this.nb[i]; k++) if (this.bp[base + k] === j) return true;
    return false;
  };

  // True when i and j share a bonded neighbour (1-3 pair in the same molecule).
  S.isOneThree = function (i, j) {
    const MB = P.MAXB, base = i * MB;
    for (let k = 0; k < this.nb[i]; k++) {
      const m = this.bp[base + k];
      const bm = m * MB;
      for (let l = 0; l < this.nb[m]; l++) if (this.bp[bm + l] === j) return true;
    }
    return false;
  };

  // ---------------------------------------------------------------------------------------
  // Spatial grid
  // ---------------------------------------------------------------------------------------
  S.updateGrid = function () {
    let maxRc = 20;
    for (let a = 1; a < NZ; a++) {
      if (!this.present[a]) continue;
      for (let b = a; b < NZ; b++) {
        if (!this.present[b]) continue;
        if (pairRc[a * NZ + b] > maxRc) maxRc = pairRc[a * NZ + b];
      }
    }
    this.cellSize = Math.ceil(maxRc);
    this.gridW = Math.ceil(P.W / this.cellSize) + 1;
    this.gridH = Math.ceil(P.H / this.cellSize) + 1;
    const cells = this.gridW * this.gridH;
    if (this.cellStart.length < cells + 1) this.cellStart = new Int32Array(cells + 1);
  };

  S.buildGrid = function () {
    const n = this.n, cs = this.cellSize, gw = this.gridW, gh = this.gridH;
    const cells = gw * gh;
    const start = this.cellStart, items = this.cellItems, cellOf = this.cellOf;
    start.fill(0, 0, cells + 1);
    for (let i = 0; i < n; i++) {
      let cx = Math.floor(this.x[i] / cs), cy = Math.floor(this.y[i] / cs);
      if (cx < 0) cx = 0; else if (cx >= gw) cx = gw - 1;
      if (cy < 0) cy = 0; else if (cy >= gh) cy = gh - 1;
      const c = cy * gw + cx;
      cellOf[i] = c;
      start[c + 1]++;
    }
    for (let c = 0; c < cells; c++) start[c + 1] += start[c];
    const fill = this._fill && this._fill.length >= cells ? this._fill : (this._fill = new Int32Array(cells));
    fill.set(start.subarray(0, cells));
    for (let i = 0; i < n; i++) items[fill[cellOf[i]]++] = i;
  };

  // Calls fn(j) for particles within radius r of (x, y).
  S.query = function (x, y, r, fn) {
    const cs = this.cellSize, gw = this.gridW, gh = this.gridH;
    const x0 = Math.max(0, Math.floor((x - r) / cs)), x1 = Math.min(gw - 1, Math.floor((x + r) / cs));
    const y0 = Math.max(0, Math.floor((y - r) / cs)), y1 = Math.min(gh - 1, Math.floor((y + r) / cs));
    for (let cy = y0; cy <= y1; cy++) {
      for (let cx = x0; cx <= x1; cx++) {
        const c = cy * gw + cx;
        for (let k = this.cellStart[c]; k < this.cellStart[c + 1]; k++) {
          const j = this.cellItems[k];
          if (j >= this.n) continue;
          const dx = this.x[j] - x, dy = this.y[j] - y;
          if (dx * dx + dy * dy <= r * r) fn(j, Math.sqrt(dx * dx + dy * dy));
        }
      }
    }
  };

  S.nearest = function (x, y, r) {
    let best = -1, bd = Infinity;
    this.query(x, y, r, function (j, d) { if (d < bd) { bd = d; best = j; } });
    return best;
  };

  // ---------------------------------------------------------------------------------------
  // Forces
  // ---------------------------------------------------------------------------------------
  S.computeForces = function () {
    const n = this.n, x = this.x, y = this.y, fx = this.fx, fy = this.fy, el = this.el;
    const q = this.q, nb = this.nb, contacts = this.contacts, ghostT = this.ghostT, used = this.used;
    const gw = this.gridW, gh = this.gridH, start = this.cellStart, items = this.cellItems;
    fx.fill(0, 0, n); fy.fill(0, 0, n);
    contacts.fill(0, 0, n);
    this.nCand = 0;
    const maxCand = this.candI.length;
    const vdwScale = P.VDW_SCALE, metalScale = 1;
    const molPol = P.MOL_POLARIZABILITY, hbEps = P.HB_EPS * P.KB, hbd = this.hbd;
    const self = this;

    function pair(i, j) {
      const dx = x[j] - x[i], dy = y[j] - y[i];
      const d2 = dx * dx + dy * dy;
      const ei = el[i], ej = el[j], p = ei * NZ + ej;
      const rc = pairRc[p];
      if (d2 >= rc * rc) return;
      if ((nb[i] && nb[j]) && (self.isBonded(i, j) || self.isOneThree(i, j))) return;
      const d = Math.sqrt(d2) || 1e-6;
      let sigma = pairSigma[p];
      // hydrogen bond: donor H and acceptor N/O/F may approach closer and attract strongly
      const hbs = hbd[i] ? HB_ACCEPTOR[ej] : hbd[j] ? HB_ACCEPTOR[ei] : 0;
      if (hbs) sigma *= P.HB_CONTACT;
      if ((ghostT[i] || ghostT[j]) && self.isGhost(i, j)) {
        if (d < sigma) return;
        self.clearGhost(i, j);
      }
      let F = 0; // positive = repulsive
      if (d < sigma) {
        const mu = mass[ei] * mass[ej] / (mass[ei] + mass[ej]);
        F = P.K_REP * mu * (sigma - d);
        contacts[i]++; contacts[j]++;
        if (self.nCand < maxCand && C.pairClass[p]) {
          self.candI[self.nCand] = i; self.candJ[self.nCand] = j; self.nCand++;
        }
      } else {
        const w = rc - sigma, s = (d - sigma) / w;
        let eps;
        if (isMetal[ei] && isMetal[ej] && !nb[i] && !nb[j]) eps = pairMetalEps[p] * metalScale;
        else {
          const bi = nb[i] > 2 ? 2 : nb[i], bj = nb[j] > 2 ? 2 : nb[j];
          eps = pairEps[p] * vdwScale * (1 + molPol * bi) * (1 + molPol * bj);
          if (hbs) eps += hbs * hbEps;
        }
        F = -4 * eps * s * (1 - s * s) / w;
        if (s < 0.25) { contacts[i]++; contacts[j]++; }
        // radicals capture partners a little before hard contact (larger reactive cross-section)
        if (s < P.RADICAL_CAPTURE && self.nCand < maxCand && C.pairClass[p] &&
            (used[i] < maxVal[ei] || used[j] < maxVal[ej])) {
          self.candI[self.nCand] = i; self.candJ[self.nCand] = j; self.nCand++;
        }
      }
      const qq = q[i] * q[j];
      if (qq !== 0) {
        // shifted-potential Coulomb, U = KC q q (1/d - 1/rc): keeps the full short-range
        // (hydrogen-bond / ionic) attraction while still going to zero at the cutoff
        const dd = d < sigma * 0.8 ? sigma * 0.8 : d;
        F += P.KC * qq / (dd * dd);
      }
      const fxp = F * dx / d, fyp = F * dy / d;
      fx[i] -= fxp; fy[i] -= fyp;
      fx[j] += fxp; fy[j] += fyp;
    }

    for (let cy = 0; cy < gh; cy++) {
      for (let cx = 0; cx < gw; cx++) {
        const c = cy * gw + cx;
        const s0 = start[c], s1 = start[c + 1];
        if (s0 === s1) continue;
        for (let a = s0; a < s1; a++) {
          const i = items[a];
          for (let b = a + 1; b < s1; b++) pair(i, items[b]);
          // forward half-stencil: right, below-left, below, below-right
          if (cx + 1 < gw) { const c2 = c + 1; for (let b = start[c2]; b < start[c2 + 1]; b++) pair(i, items[b]); }
          if (cy + 1 < gh) {
            const r = c + gw;
            if (cx > 0) for (let b = start[r - 1]; b < start[r]; b++) pair(i, items[b]);
            for (let b = start[r]; b < start[r + 1]; b++) pair(i, items[b]);
            if (cx + 1 < gw) for (let b = start[r + 1]; b < start[r + 2]; b++) pair(i, items[b]);
          }
        }
      }
    }
  };

  // ---------------------------------------------------------------------------------------
  // Integration, constraints, walls
  // ---------------------------------------------------------------------------------------
  S.integrate = function () {
    const n = this.n, x = this.x, y = this.y, px = this.px, py = this.py, el = this.el;
    const dt = P.DT, dt2 = dt * dt, maxS = P.MAX_STEP;
    const contacts = this.contacts;
    const tgw = this.tgw, tgh = this.tgh, tSum = this.tSum, tCnt = this.tCnt, tCell = this.tCell;
    tSum.fill(0); tCnt.fill(0);
    const cellKE = this.cellKE, cellDebt = this.cellDebt;
    cellKE.fill(0);
    const W = P.W, H = P.H, rng = this.rng, wallG = this.insulated ? 0 : P.WALL_GAMMA;
    let ke = 0, bonds = 0, impulse = 0;
    // bulk (flow) velocity of condensed matter per coarse cell
    const mom = this.cellMom, TC = P.TCELL;
    mom.fill(0);
    for (let i = 0; i < n; i++) {
      if (contacts[i] < 2) continue;
      let tcx = Math.floor(x[i] / TC), tcy = Math.floor(y[i] / TC);
      tcx = tcx < 0 ? 0 : tcx >= tgw ? tgw - 1 : tcx;
      tcy = tcy < 0 ? 0 : tcy >= tgh ? tgh - 1 : tcy;
      const c3 = (tcy * tgw + tcx) * 3, m = mass[el[i]];
      mom[c3] += m * (x[i] - px[i]); mom[c3 + 1] += m * (y[i] - py[i]); mom[c3 + 2] += m;
    }
    for (let i = 0; i < n; i++) {
      let vx = x[i] - px[i], vy = y[i] - py[i];
      const v2 = vx * vx + vy * vy;
      if (v2 > maxS * maxS) { const s = maxS / Math.sqrt(v2); vx *= s; vy *= s; }
      const im = invMass[el[i]];
      const c = contacts[i] > 4 ? 4 : contacts[i];
      const g = P.G_GAS + (c >= 2 ? P.G_CONDENSED * c * 0.25 : 0);
      // cap the local bulk flow speed (not thermal motion) of condensed matter, lower near the
      // floor, so settling liquids/solids land softly instead of heating up on impact
      let bvx = 0, bvy = 0;
      if (c >= 2) {
        let tcx = Math.floor(x[i] / TC), tcy = Math.floor(y[i] / TC);
        tcx = tcx < 0 ? 0 : tcx >= tgw ? tgw - 1 : tcx;
        tcy = tcy < 0 ? 0 : tcy >= tgh ? tgh - 1 : tcy;
        const c3 = (tcy * tgw + tcx) * 3;
        if (mom[c3 + 2] > 0) {
          bvx = mom[c3] / mom[c3 + 2]; bvy = mom[c3 + 1] / mom[c3 + 2];
          const floorDist = H - y[i];
          // never act on the cell's random thermal drift (~ sqrt(kT / M_cell))
          const noise = 4 * Math.sqrt(P.KB * this.tLoc[i] / mom[c3 + 2]);
          const vmaxB = Math.max(noise, floorDist < 250 ? P.BULK_VMAX_FLOOR + (P.BULK_VMAX - P.BULK_VMAX_FLOOR) * floorDist / 250 : P.BULK_VMAX) * dt;
          const vb = Math.sqrt(bvx * bvx + bvy * bvy);
          if (vb > vmaxB) { const f = 1 - vmaxB / vb; vx -= bvx * f; vy -= bvy * f; }
        }
      }
      // soft walls (energy conserving with bond constraints, unlike reflections)
      const rw = radius[el[i]] * 0.6;
      let wfx = 0, wfy = 0, Tw = -1, gamma = 0;
      const WZ = P.WALL_ZONE;
      if (x[i] < WZ || x[i] > W - WZ || y[i] < WZ) { Tw = P.AMBIENT; gamma = wallG; }
      if (x[i] < rw) wfx = P.K_WALL * (rw - x[i]);
      else if (x[i] > W - rw) wfx = -P.K_WALL * (x[i] - W + rw);
      if (y[i] < rw) {
        wfy = P.K_WALL * (rw - y[i]);
        impulse += wfy * mass[el[i]] * dt;
      } else if (y[i] > H - P.PLATE_ZONE) {
        if (y[i] > H - rw) wfy = -P.K_WALL * (y[i] - H + rw);
        Tw = this.plateTemp; gamma = P.PLATE_GAMMA;
      }
      if (wfx !== 0 && y[i] < H / 2) impulse += Math.abs(wfx) * mass[el[i]] * dt;
      if (this.anneal[i]) {
        // freshly added material: remove the heat of settling / self-assembly
        this.anneal[i]--;
        Tw = P.AMBIENT; gamma = 40;
      } else if (gamma === 0 && this.convection > 0) {
        // convection: everything exchanges heat with the plate, more strongly near the bottom
        const h = y[i] / H;
        Tw = this.plateTemp; gamma = this.convection * (0.15 + 0.85 * h * h);
      }
      if (gamma > 0) {
        // Langevin (Ornstein-Uhlenbeck) thermostat towards the wall temperature
        const c1 = Math.exp(-gamma * dt);
        const sd = Math.sqrt((1 - c1 * c1) * P.KB * Tw * im) * dt;
        // cheap unit-variance noise (sum of 3 uniforms); exact Gaussians are not needed here
        vx = vx * c1 + sd * 2 * (rng() + rng() + rng() - 1.5);
        vy = vy * c1 + sd * 2 * (rng() + rng() + rng() - 1.5);
      }
      px[i] = x[i]; py[i] = y[i];
      x[i] += vx + (this.fx[i] * im + wfx) * dt2;
      y[i] += vy + (this.fy[i] * im + g + wfy) * dt2;
      // local temperature field: accumulate kinetic energy per coarse cell. Each bond removes
      // one degree of freedom shared by its two atoms, so scale up accordingly.
      const rvx = vx - bvx, rvy = vy - bvy;
      const tInst = mass[el[i]] * (rvx * rvx + rvy * rvy) / (dt2 * 2 * P.KB);
      const dofScale = 2 / Math.max(1, 2 - 0.5 * this.nb[i]);
      let tcx = Math.floor(x[i] / P.TCELL), tcy = Math.floor(y[i] / P.TCELL);
      tcx = tcx < 0 ? 0 : tcx >= tgw ? tgw - 1 : tcx;
      tcy = tcy < 0 ? 0 : tcy >= tgh ? tgh - 1 : tcy;
      const tc = tcy * tgw + tcx;
      this.cellOfT[i] = tc;
      tSum[tc] += tInst * dofScale; tCnt[tc] += 1;
      this.tSelf[i] += P.SELF_T_RATE * (tInst * dofScale - this.tSelf[i]);
      cellKE[tc] += tInst * P.KB;
      ke += tInst;
      bonds += this.nb[i];
    }
    // T = sum(KE) / (KB * DOF / 2) with DOF = 2 per atom minus one per bond constraint
    const globalT = n ? ke / Math.max(1, n - bonds / 4) : this.plateTemp;
    // local temperature: 3x3 blurred cell average, shrunk toward the global temperature when
    // few particles are present (one fast gas atom is not a hot spot), smoothed in time
    const N0 = P.T_PRIOR;
    const tNew = this._tNew || (this._tNew = new Float32Array(tCell.length));
    for (let cy = 0; cy < tgh; cy++) {
      for (let cx = 0; cx < tgw; cx++) {
        let sum = 0, cnt = 0;
        for (let oy = -1; oy <= 1; oy++) {
          const yy = cy + oy;
          if (yy < 0 || yy >= tgh) continue;
          for (let ox = -1; ox <= 1; ox++) {
            const xx = cx + ox;
            if (xx < 0 || xx >= tgw) continue;
            const w = (ox === 0 && oy === 0) ? 1 : 0.35;
            sum += tSum[yy * tgw + xx] * w; cnt += tCnt[yy * tgw + xx] * w;
          }
        }
        tNew[cy * tgw + cx] = (sum + N0 * globalT) / (cnt + N0);
      }
    }
    for (let c = 0; c < tCell.length; c++) tCell[c] += 0.3 * (tNew[c] - tCell[c]);
    // kinetics see a blend of the regional field and the atom's own recent temperature, so a
    // small hot spot (a flame front in a dilute gas) is not averaged away
    const ws = P.SELF_T_WEIGHT;
    for (let i = 0; i < n; i++) {
      // only in the gas phase: in dense matter the regional field is already accurate
      const w = contacts[i] >= 2 ? 0 : ws;
      this.tLoc[i] = (1 - w) * tCell[this.cellOfT[i]] + w * this.tSelf[i];
    }
    // repay local energy debt by gently cooling the particles of indebted cells
    let anyDebt = false;
    for (let c = 0; c < cellDebt.length; c++) if (cellDebt[c] > 0) { anyDebt = true; break; }
    if (anyDebt) {
      const fac = this._fac || (this._fac = new Float32Array(cellDebt.length));
      for (let c = 0; c < cellDebt.length; c++) {
        fac[c] = 1;
        if (cellDebt[c] <= 0 || cellKE[c] <= 0) continue;
        const pay = Math.min(cellDebt[c], cellKE[c] * 0.2);
        cellDebt[c] -= pay;
        fac[c] = Math.sqrt(1 - pay / cellKE[c]);
      }
      for (let i = 0; i < n; i++) {
        const f = fac[this.cellOfT[i]];
        if (f !== 1) { px[i] = x[i] - (x[i] - px[i]) * f; py[i] = y[i] - (y[i] - py[i]) * f; }
      }
    }
    this.impulse += impulse;
    this._keSum = (this._keSum || 0) + ke;
    this._keCount = (this._keCount || 0) + Math.max(1, n - bonds / 4);
  };

  // SHAKE-style distance constraint: the correction acts along the bond direction of the
  // previous step, which (unlike projecting along the new direction) conserves rotational energy.
  S.shake = function (a, b, target, stiffness) {
    const x = this.x, y = this.y, px = this.px, py = this.py;
    const wa = invMass[this.el[a]], wb = invMass[this.el[b]];
    const rx = x[b] - x[a], ry = y[b] - y[a];
    let ox = px[b] - px[a], oy = py[b] - py[a];
    const r2 = rx * rx + ry * ry;
    let dot = rx * ox + ry * oy;
    if (dot * dot < 0.25 * r2 * (ox * ox + oy * oy) || dot <= 0) { ox = rx; oy = ry; dot = r2; }
    if (dot < 1e-9) return;
    const g = stiffness * (target * target - r2) / (2 * dot * (wa + wb));
    x[a] -= g * wa * ox; y[a] -= g * wa * oy;
    x[b] += g * wb * ox; y[b] += g * wb * oy;
  };

  S.solveBonds = function (iterations) {
    const n = this.n, MB = P.MAXB, el = this.el;
    for (let it = 0; it < iterations; it++) {
      for (let i = 0; i < n; i++) {
        const base = i * MB, ei = el[i];
        for (let k = 0; k < this.nb[i]; k++) {
          const j = this.bp[base + k];
          if (j < i) continue;
          this.shake(i, j, this.bl[base + k], 1);
        }
      }
    }
  };

  // Moves each bond's target length toward its rest length. The atoms are moved kinematically
  // (position and previous position together) so the change does not create velocity.
  S.rampBonds = function () {
    const n = this.n, MB = P.MAXB, el = this.el, bl = this.bl, ramp = P.BOND_RAMP;
    const x = this.x, y = this.y, px = this.px, py = this.py;
    for (let i = 0; i < n; i++) {
      const base = i * MB;
      for (let k = 0; k < this.nb[i]; k++) {
        const j = this.bp[base + k];
        if (j < i) continue;
        const rest = C.pairLength[(el[i] * NZ + el[j]) * 4 + this.bo[base + k]];
        const cur = bl[base + k];
        if (cur === rest) continue;
        const next = Math.abs(cur - rest) <= ramp ? rest : cur - Math.sign(cur - rest) * ramp;
        bl[base + k] = next;
        const kj = this.bondSlot(j, i);
        if (kj >= 0) bl[j * MB + kj] = next;
        const dx = x[j] - x[i], dy = y[j] - y[i];
        const d = Math.sqrt(dx * dx + dy * dy) || 1e-6;
        const wi = invMass[el[i]], wj = invMass[el[j]];
        const sh = (d - next) / (d * (wi + wj));
        const sxi = dx * sh * wi, syi = dy * sh * wi, sxj = dx * sh * wj, syj = dy * sh * wj;
        x[i] += sxi; y[i] += syi; px[i] += sxi; py[i] += syi;
        x[j] -= sxj; y[j] -= syj; px[j] -= sxj; py[j] -= syj;
      }
    }
  };

  // VSEPR bond angles as conservative springs between the neighbours of each centre (a bending
  // potential). Two-coordinate centres get their VSEPR angle; higher coordination spreads evenly.
  S.angleForces = function () {
    const n = this.n, MB = P.MAXB, x = this.x, y = this.y, el = this.el, fx = this.fx, fy = this.fy;
    const angs = this._angs || (this._angs = new Float64Array(MB));
    const ids = this._ids || (this._ids = new Int32Array(MB));
    const K = P.K_ANGLE;
    function spring(a, b, target) {
      const dx = x[b] - x[a], dy = y[b] - y[a];
      const d = Math.sqrt(dx * dx + dy * dy) || 1e-6;
      const ma = mass[el[a]], mb = mass[el[b]];
      const F = K * (ma * mb / (ma + mb)) * (d - target) / d;
      fx[a] += F * dx; fy[a] += F * dy;
      fx[b] -= F * dx; fy[b] -= F * dy;
    }
    for (let c = 0; c < n; c++) {
      const m = this.nb[c];
      if (m < 2) continue;
      const base = c * MB, ec = el[c];
      if (m === 2) {
        const theta = angleFor(ec, this.used[c]);
        const la = this.bl[base], lb = this.bl[base + 1];
        spring(this.bp[base], this.bp[base + 1], Math.sqrt(la * la + lb * lb - 2 * la * lb * Math.cos(theta)));
        continue;
      }
      // sort neighbours by angle around the centre, spring consecutive pairs to 360/m
      for (let k = 0; k < m; k++) {
        const j = this.bp[base + k];
        ids[k] = k;
        angs[k] = Math.atan2(y[j] - y[c], x[j] - x[c]);
      }
      for (let k = 1; k < m; k++) { // insertion sort (m <= 8)
        const t = ids[k];
        let l = k - 1;
        while (l >= 0 && angs[ids[l]] > angs[t]) { ids[l + 1] = ids[l]; l--; }
        ids[l + 1] = t;
      }
      const cphi = Math.cos(2 * Math.PI / m);
      for (let k = 0; k < m; k++) {
        const ka = ids[k], kb = ids[(k + 1) % m];
        const la = this.bl[base + ka], lb = this.bl[base + kb];
        spring(this.bp[base + ka], this.bp[base + kb], Math.sqrt(la * la + lb * lb - 2 * la * lb * cphi));
      }
    }
  };

  const angleCache = new Float32Array(NZ * 32).fill(-1);
  function angleFor(z, used) {
    const k = z * 32 + Math.min(used, 31);
    if (angleCache[k] < 0) angleCache[k] = C.twoBondAngle(z, used);
    return angleCache[k];
  }

  // Vacuum valve and a hard safety boundary (the soft walls in integrate() do the real work).
  S.walls = function () {
    const n = this.n, x = this.x, y = this.y, px = this.px, py = this.py;
    const W = P.W, H = P.H;
    const valveL = W / 2 - P.VALVE_HALF, valveR = W / 2 + P.VALVE_HALF;
    const toRemove = this._toRemove || (this._toRemove = []);
    toRemove.length = 0;
    for (let i = 0; i < n; i++) {
      if (this.vacuumOn && y[i] < radius[this.el[i]] && x[i] > valveL && x[i] < valveR) { toRemove.push(i); continue; }
      if (x[i] < -10) { px[i] += -10 - x[i]; x[i] = -10; }
      else if (x[i] > W + 10) { px[i] -= x[i] - W - 10; x[i] = W + 10; }
      if (y[i] < -10) { py[i] += -10 - y[i]; y[i] = -10; }
      else if (y[i] > H + 10) { py[i] -= y[i] - H - 10; y[i] = H + 10; }
    }
    if (toRemove.length) {
      // remove whole molecules sucked into the valve
      const atoms = [];
      for (let k = 0; k < toRemove.length; k++) {
        const i = toRemove[k];
        if (i < this.n) atoms.push.apply(atoms, this.moleculeOf(i, []));
      }
      const uniq = Array.from(new Set(atoms)).sort(function (a, b) { return b - a; });
      for (let k = 0; k < uniq.length; k++) if (uniq[k] < this.n) this.removeParticle(uniq[k]);
      this.removedMol += uniq.length * P.MOL_PER_PARTICLE;
    }
  };


  // ---------------------------------------------------------------------------------------
  // Chemistry: collisions -> reactions, thermal dissociation, bond order changes
  // ---------------------------------------------------------------------------------------

  // Deposits (dE > 0) or withdraws (dE < 0) reaction enthalpy (kJ/mol) as kinetic energy.
  S.heat = function (i, j, dE) {
    this.heatRaw(i, j, dE * P.HEAT_FRACTION * ENERGY_PER_KJ);
  };

  // Same as heat() with the energy already in simulation units.
  S.heatRaw = function (i, j, e) {
    const dt = P.DT, x = this.x, y = this.y, px = this.px, py = this.py;
    const mi = mass[this.el[i]], mj = mass[this.el[j]];
    let vix = (x[i] - px[i]) / dt, viy = (y[i] - py[i]) / dt;
    let vjx = (x[j] - px[j]) / dt, vjy = (y[j] - py[j]) / dt;
    const M = mi + mj;
    const cx = (mi * vix + mj * vjx) / M, cy = (mi * viy + mj * vjy) / M;
    let rix = vix - cx, riy = viy - cy, rjx = vjx - cx, rjy = vjy - cy;
    const keRel = 0.5 * mi * (rix * rix + riy * riy) + 0.5 * mj * (rjx * rjx + rjy * rjy);
    if (e >= 0) {
      if (keRel < 1e-9) {
        const a = this.rng() * Math.PI * 2;
        rix = Math.cos(a) / mi; riy = Math.sin(a) / mi; rjx = -Math.cos(a) / mj; rjy = -Math.sin(a) / mj;
        const k0 = 0.5 * mi * (rix * rix + riy * riy) + 0.5 * mj * (rjx * rjx + rjy * rjy);
        const s = Math.sqrt(e / k0);
        rix *= s; riy *= s; rjx *= s; rjy *= s;
      } else {
        const s = Math.sqrt((keRel + e) / keRel);
        rix *= s; riy *= s; rjx *= s; rjy *= s;
      }
    } else {
      // take what the pair has (keeping 10%); the rest becomes a global energy debt that is
      // drained from all particles, so energy is conserved overall
      const take = Math.min(-e, keRel * 0.9);
      this.cellDebt[this.cellOfT[i]] += -e - take;
      const s = Math.sqrt((keRel - take) / Math.max(keRel, 1e-9));
      rix *= s; riy *= s; rjx *= s; rjy *= s;
    }
    px[i] = x[i] - (cx + rix) * dt; py[i] = y[i] - (cy + riy) * dt;
    px[j] = x[j] - (cx + rjx) * dt; py[j] = y[j] - (cy + rjy) * dt;
    if (e > 0) { this.eRef[i] += e * mj / M; this.eRef[j] += e * mi / M; }
  };

  // Non-bonded potential energy between i and j, mirroring the forces in computeForces().
  S.pairPotential = function (i, j) {
    const ei = this.el[i], ej = this.el[j], p = ei * NZ + ej;
    const dx = this.x[j] - this.x[i], dy = this.y[j] - this.y[i];
    const d = Math.sqrt(dx * dx + dy * dy) || 1e-6;
    const rc = pairRc[p];
    if (d >= rc) return 0;
    if (this.isBonded(i, j) || this.isOneThree(i, j)) return 0;
    const nb = this.nb;
    let sigma = pairSigma[p];
    const hbs = this.hbd[i] ? HB_ACCEPTOR[ej] : this.hbd[j] ? HB_ACCEPTOR[ei] : 0;
    if (hbs) sigma *= P.HB_CONTACT;
    if ((this.ghostT[i] || this.ghostT[j]) && this.isGhost(i, j) && d < sigma) return 0;
    let eps;
    if (isMetal[ei] && isMetal[ej] && !nb[i] && !nb[j]) eps = pairMetalEps[p];
    else {
      const bi = nb[i] > 2 ? 2 : nb[i], bj = nb[j] > 2 ? 2 : nb[j];
      eps = pairEps[p] * P.VDW_SCALE * (1 + P.MOL_POLARIZABILITY * bi) * (1 + P.MOL_POLARIZABILITY * bj);
      if (hbs) eps += hbs * P.HB_EPS * P.KB;
    }
    let U;
    if (d < sigma) {
      const mu = mass[ei] * mass[ej] / (mass[ei] + mass[ej]);
      U = -eps + 0.5 * P.K_REP * mu * (sigma - d) * (sigma - d);
    } else {
      const s = (d - sigma) / (rc - sigma);
      U = -eps * (1 - s * s) * (1 - s * s);
    }
    const qq = this.q[i] * this.q[j];
    if (qq !== 0) {
      const d0 = sigma * 0.8;
      U += d >= d0 ? P.KC * qq * (1 / d - 1 / rc) : P.KC * qq * ((1 / d0 - 1 / rc) + (d0 - d) / (d0 * d0));
    }
    return U;
  };

  // Non-bonded energy of a set of atoms with everything around them (pairs inside counted once).
  S.eventEnergy = function (set, count) {
    let U = 0;
    const self = this, r = this.cellSize;
    for (let a = 0; a < count; a++) {
      const i = set[a];
      this.query(this.x[i], this.y[i], r, function (j) {
        if (j === i) return;
        for (let b = 0; b < count; b++) if (set[b] === j && b < a) return; // already counted
        U += self.pairPotential(i, j);
      });
    }
    return U;
  };

  // Bending (angle spring) energy around centre c, mirroring angleForces().
  S.angleEnergy = function (c) {
    const m = this.nb[c];
    if (m < 2) return 0;
    const MB = P.MAXB, base = c * MB, x = this.x, y = this.y, el = this.el, K = P.K_ANGLE;
    const self = this;
    function springU(a, b, target) {
      const d = Math.hypot(x[b] - x[a], y[b] - y[a]);
      const ma = mass[el[a]], mb = mass[el[b]];
      return 0.5 * K * (ma * mb / (ma + mb)) * (d - target) * (d - target);
    }
    if (m === 2) {
      const theta = angleFor(el[c], this.used[c]);
      const la = this.bl[base], lb = this.bl[base + 1];
      return springU(this.bp[base], this.bp[base + 1], Math.sqrt(la * la + lb * lb - 2 * la * lb * Math.cos(theta)));
    }
    const order = [];
    for (let k = 0; k < m; k++) {
      const j = this.bp[base + k];
      order.push({ k: k, a: Math.atan2(y[j] - y[c], x[j] - x[c]) });
    }
    order.sort(function (p, q) { return p.a - q.a; });
    const cphi = Math.cos(2 * Math.PI / m);
    let U = 0;
    for (let k = 0; k < m; k++) {
      const ka = order[k].k, kb = order[(k + 1) % m].k;
      const la = self.bl[base + ka], lb = self.bl[base + kb];
      U += springU(self.bp[base + ka], self.bp[base + kb], Math.sqrt(la * la + lb * lb - 2 * la * lb * cphi));
    }
    return U;
  };

  S.eventTotalEnergy = function (set, count) {
    let U = this.eventEnergy(set, count);
    if (P.ACCOUNT_ANGLES) for (let a = 0; a < count; a++) U += this.angleEnergy(set[a]);
    return U;
  };

  // Records the non-bonded energy of the atoms about to react...
  S.beginEvent = function (a, b, c, d) {
    const set = this._evSet || (this._evSet = new Int32Array(4));
    let k = 0;
    set[k++] = a; set[k++] = b;
    if (c >= 0) set[k++] = c;
    if (d >= 0 && d !== c) set[k++] = d;
    this._evCount = k;
    this._evU = this.eventTotalEnergy(set, k);
  };
  // ...and afterwards turns the change (charges switching, exclusions changing) into heat, so
  // every reaction conserves energy exactly.
  S.endEvent = function (i, j) {
    const dU = this.eventTotalEnergy(this._evSet, this._evCount) - this._evU;
    if (dU !== 0) this.heatRaw(i, j, -dU);
  };

  // Can the neighbourhood of i, j supply dE kJ/mol (endothermic event)? Energy comes first from
  // the pair's own motion and then from the surrounding cell; events nobody can pay for don't happen.
  S.canAfford = function (i, j, dE) {
    if (dE <= 0) return true;
    const need = dE * P.HEAT_FRACTION * ENERGY_PER_KJ;
    const c = this.cellOfT[i];
    return need < (this.cellKE[c] - this.cellDebt[c]) * 0.5;
  };

  // Weakest bond unit of atom a, excluding partner `skip`. Returns partner index or -1.
  S.weakestBond = function (a, skip) {
    const MB = P.MAXB, base = a * MB, ea = this.el[a];
    let best = -1, bestE = Infinity, spin = false;
    for (let k = 0; k < this.nb[a]; k++) {
      const b = this.bp[base + k];
      if (b === skip) continue;
      const s = C.stepEnergy(ea, this.el[b], this.bo[base + k]);
      if (s < bestE) { bestE = s; best = b; spin = ea === 8 && this.el[b] === 8 && this.bo[base + k] === 2; }
    }
    this._weakE = bestE;
    this._weakSpin = spin; // the unit to break is the pi bond of triplet O2
    return best;
  };

  S.logReaction = function (text) {
    this.reactionCount++;
    this.reactionWindow++;
    if (this.lastReactions.length < 8 && text) this.lastReactions.push(text);
  };

  S.react = function () {
    const rng = this.rng, el = this.el;
    for (let c = 0; c < this.nCand; c++) {
      const i = this.candI[c], j = this.candJ[c];
      if (i >= this.n || j >= this.n || this.isBonded(i, j)) continue;
      const ei = el[i], ej = el[j], p = ei * NZ + ej;
      const cls = C.pairClass[p];
      if (!cls) continue;
      const T = 0.5 * (this.tLoc[i] + this.tLoc[j]);
      const RT = P.R_KJ * Math.max(T, 1);
      const E1 = C.pairEnergy[p * 4 + 1];
      const canI = this.canAccept(i, ej), canJ = this.canAccept(j, ei);
      if (canI && canJ) {
        // radical + radical: nearly barrierless recombination
        if (rng() < Math.exp(Math.min(0, P.LN_COLL - 4 / RT))) {
          this.beginEvent(i, j, -1, -1);
          this.addBond(i, j, 1, cls);
          this.endEvent(i, j);
          if (!this.sameBatch(i, j)) this.heat(i, j, E1);
          this.logReaction();
        }
        continue;
      }
      if (canI !== canJ) {
        // radical abstraction: a free valence on one atom steals a bond unit from the other
        const a = canI ? i : j, b = canI ? j : i; // a has free valence, b must give one up
        const k = this.weakestBond(b, a);
        if (k < 0 || el[k] === el[a]) continue; // X. + X-Y -> X-Y + X. changes nothing
        const step = this._weakE;
        // after losing a unit, b must be able to bond with a (hypervalence rule)
        if (this.used[b] > normVal[el[b]] && !(en[el[a]] > en[el[b]] + 0.3)) continue;
        const dH = step - E1;
        const spin = this._weakSpin && isMetal[el[a]] ? P.SPIN_BARRIER : 0;
        const Ea = P.ALPHA * step + Math.max(0, dH) + spin;
        if (rng() < Math.exp(Math.min(0, P.LN_COLL - Ea / RT)) && this.canAfford(a, b, dH)) {
          this.beginEvent(a, b, k, -1);
          this.setBond(b, k, this.bondOrder(b, k) - 1);
          this.addBond(a, b, 1, cls);
          this.endEvent(a, b);
          this.heat(a, b, -dH);
          this.logReaction();
          continue;
        }
        // atom transfer: the radical takes the whole atom off a multiple bond (H. + O=O ->
        // H-O. + .O.). Endothermic, so it only runs in hot zones, where it branches the chain.
        const ord = this.bondOrder(b, k);
        if (ord >= 2) {
          const full = C.bondEnergy(el[b], el[k], ord);
          const dHt = full - E1;
          const Eat = P.ALPHA * full + Math.max(0, dHt) + spin;
          if (Eat / RT < P.LN_COLL + 12 && rng() < Math.exp(Math.min(0, P.LN_COLL - Eat / RT)) && this.canAfford(a, b, dHt)) {
            this.beginEvent(a, b, k, -1);
            this.setBond(b, k, 0);
            this.addBond(a, b, 1, cls);
            this.endEvent(a, b);
            this.heat(a, b, -dHt);
            this.logReaction();
          }
        }
        continue;
      }
      // both saturated: four-centre exchange AB + CD -> AC + BD (or AC + B. + D.)
      const k = this.weakestBond(i, j);
      if (k < 0) continue;
      const sI = this._weakE, spinI = this._weakSpin;
      const l = this.weakestBond(j, i);
      if (l < 0) continue;
      const sJ = this._weakE, spin4 = (spinI || this._weakSpin) ? P.SPIN_BARRIER : 0;
      // AB + CD -> AC + BD with identical bond types (e.g. I2 + I2) changes nothing
      if (el[i] === el[l] && el[j] === el[k]) continue;
      if (el[i] === el[j] && el[k] === el[l]) continue;
      const breakSum = sI + sJ;
      let formSum = E1, formKL = false;
      if (k !== l && !this.isBonded(k, l)) {
        const pkl = el[k] * NZ + el[l];
        const dx = this.x[l] - this.x[k], dy = this.y[l] - this.y[k];
        const near = dx * dx + dy * dy < 2.25 * pairSigma[pkl] * pairSigma[pkl];
        if (C.pairClass[pkl] && near) { formSum += C.pairEnergy[pkl * 4 + 1]; formKL = true; }
      }
      const dH = breakSum - formSum;
      const Ea = P.ALPHA_4C * breakSum + Math.max(0, dH) + spin4;
      if (Ea / RT > P.LN_COLL + 12) continue;
      if (rng() < Math.exp(Math.min(0, P.LN_COLL - Ea / RT)) && this.canAfford(i, j, dH)) {
        this.beginEvent(i, j, k, l);
        this.setBond(i, k, this.bondOrder(i, k) - 1);
        this.setBond(j, l, this.bondOrder(j, l) - 1);
        this.addBond(i, j, 1, cls);
        if (formKL) this.addBond(k, l, 1, C.pairClass[el[k] * NZ + el[l]]);
        this.endEvent(i, j);
        this.heat(i, j, -dH);
        this.logReaction();
      }
    }
  };

  // Thermal dissociation (Arrhenius) and pi-bond formation between radicals, per bond.
  S.bondEvents = function (dtEff) {
    const n = this.n, MB = P.MAXB, el = this.el, rng = this.rng;
    const lnDt = Math.log(dtEff);
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < this.nb[i]; k++) {
        const j = this.bp[i * MB + k];
        if (j < i) continue;
        const ei = el[i], ej = el[j], p = ei * NZ + ej;
        const order = this.bo[i * MB + k];
        const T = 0.5 * (this.tLoc[i] + this.tLoc[j]);
        const RT = P.R_KJ * Math.max(T, 1);
        // upgrade bond order when both atoms still have free valence
        if (order < C.pairMaxOrder[p] && this.canAccept(i, ej) && this.canAccept(j, ei)) {
          if (rng() < Math.exp(Math.min(0, P.LN_COLL - 8 / RT))) {
            const gain = C.stepEnergy(ei, ej, order + 1);
            this.beginEvent(i, j, -1, -1);
            this.setBond(i, j, order + 1);
            this.endEvent(i, j);
            if (!this.sameBatch(i, j)) this.heat(i, j, gain);
            this.logReaction();
            continue;
          }
        }
        // dissociation of the top bond unit
        let step = C.stepEnergy(ei, ej, order);
        if (step <= 0) step = 1;
        const x = P.LN_UNI + lnDt - step / RT;
        if (x < -12) continue;
        if (rng() < Math.exp(Math.min(0, x)) && this.canAfford(i, j, step)) {
          this.beginEvent(i, j, -1, -1);
          this.setBond(i, j, order - 1);
          this.endEvent(i, j);
          this.heat(i, j, -step);
          this.logReaction();
          if (order === 1) k--; // slot was removed; revisit the one swapped into it
        }
      }
    }
  };

  // Mechanical failure: bonds stretched far beyond their length snap.
  S.breakStretched = function () {
    const n = this.n, MB = P.MAXB, x = this.x, y = this.y, el = this.el;
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < this.nb[i]; k++) {
        const j = this.bp[i * MB + k];
        if (j < i) continue;
        const rest = C.pairLength[(el[i] * NZ + el[j]) * 4 + this.bo[i * MB + k]];
        const dx = x[j] - x[i], dy = y[j] - y[i];
        if (dx * dx + dy * dy > 9 * rest * rest) { this.setBond(i, j, 0); k--; }
      }
    }
  };

  // ---------------------------------------------------------------------------------------
  // Pressure can
  // ---------------------------------------------------------------------------------------
  S.feedGas = function (dt) {
    if (!this.canOn) return;
    this.canAccumulator += this.canRate * dt / P.MOL_PER_PARTICLE;
    while (this.canAccumulator >= 1) {
      this.canAccumulator -= 1;
      const x = P.W / 2 + (this.rng() - 0.5) * P.VALVE_HALF * 1.6;
      const y = 30 + this.rng() * 40;
      const before = this.n;
      this.addGasMolecule(this.canGas, x, y, P.AMBIENT);
      for (let i = before; i < this.n; i++) this.py[i] = this.y[i] - 150 * P.DT; // jet downwards
    }
  };

  // ---------------------------------------------------------------------------------------
  // Main step
  // ---------------------------------------------------------------------------------------
  S.step = function (substeps) {
    for (let s = 0; s < substeps; s++) {
      this.buildGrid();
      this.computeForces();
      this.angleForces();
      this.integrate();
      this.rampBonds();
      this.solveBonds(P.SHAKE_ITERS);
      this.walls();
      this.react();
      if (this.substep % 4 === 0) {
        this.bondEvents(P.DT * 4);
        this.breakStretched();
      }
      this.feedGas(P.DT);
      this.time += P.DT;
      this.substep++;
      this.impulseTime += P.DT;
      this.reactionWindowTime += P.DT;
    }
    this.updateGauges();
  };

  S.updateGauges = function () {
    if (this._keCount > 0) {
      const T = this._keSum / this._keCount;
      this.temperature = this.n ? T : this.plateTemp;
      this._keSum = 0; this._keCount = 0;
    }
    if (this.impulseTime >= 0.15) {
      const len = P.W + P.H; // top wall + two half side walls
      const p2d = this.impulse / (this.impulseTime * len);
      const area = P.W * P.H;
      const atm = p2d * area * P.MOL_PER_PARTICLE * 0.082057 / (P.KB * P.JAR_VOLUME_L);
      this.pressureAtm = this.pressureAtm * 0.85 + atm * 0.15;
      this.impulse = 0; this.impulseTime = 0;
    }
    if (this.reactionWindowTime >= 0.5) {
      this.reactionsPerSec = this.reactionWindow / this.reactionWindowTime;
      this.reactionWindow = 0; this.reactionWindowTime = 0;
    }
  };

  // Local heating/cooling tool (spark igniter / cold finger). A hot spark also splits some bonds
  // in its path, creating the radicals that start chain reactions (as a real spark does).
  S.spark = function (x, y, r, T) {
    const self = this, MB = P.MAXB;
    this.buildGrid();
    const hit = [];
    this.query(x, y, r, function (j) { hit.push(j); });
    for (let k = 0; k < hit.length; k++) {
      const j = hit[k];
      const sd = Math.sqrt(P.KB * T / mass[self.el[j]]) * P.DT;
      self.px[j] = self.x[j] - gaussian(self.rng) * sd;
      self.py[j] = self.y[j] - gaussian(self.rng) * sd;
      self.tLoc[j] = T; self.tSelf[j] = T;
    }
    if (T > 1500) {
      for (let k = 0; k < hit.length; k++) {
        const i = hit[k];
        for (let b = 0; b < this.nb[i]; b++) {
          const j = this.bp[i * MB + b];
          if (j < i || this.rng() > 0.5) continue;
          const order = this.bo[i * MB + b];
          if (C.stepEnergy(this.el[i], this.el[j], order) > 700) continue;
          this.setBond(i, j, order - 1);
          this.logReaction();
          b--;
        }
      }
    }
    // heat the local temperature field as well
    const c0 = Math.max(0, Math.floor((x - r) / P.TCELL)), c1 = Math.min(this.tgw - 1, Math.floor((x + r) / P.TCELL));
    const r0 = Math.max(0, Math.floor((y - r) / P.TCELL)), r1 = Math.min(this.tgh - 1, Math.floor((y + r) / P.TCELL));
    for (let cy = r0; cy <= r1; cy++) for (let cx = c0; cx <= c1; cx++) this.tCell[cy * this.tgw + cx] = T;
  };

  // ---------------------------------------------------------------------------------------
  // Analysis: molecules, species, colours
  // ---------------------------------------------------------------------------------------
  S.analyze = function () {
    const n = this.n, MB = P.MAXB, el = this.el;
    const parent = this._parent && this._parent.length >= n ? this._parent : (this._parent = new Int32Array(P.MAX));
    for (let i = 0; i < n; i++) parent[i] = i;
    function find(a) { while (parent[a] !== a) { parent[a] = parent[parent[a]]; a = parent[a]; } return a; }
    for (let i = 0; i < n; i++) {
      for (let k = 0; k < this.nb[i]; k++) {
        const j = this.bp[i * MB + k];
        if (j > i) { const a = find(i), b = find(j); if (a !== b) parent[a] = b; }
      }
    }
    const mols = new Map();
    for (let i = 0; i < n; i++) {
      const r = find(i);
      this.molId[i] = r;
      let m = mols.get(r);
      if (!m) { m = { atoms: [], counts: {}, ionic: false, contacts: 0 }; mols.set(r, m); }
      m.atoms.push(i);
      m.counts[el[i]] = (m.counts[el[i]] || 0) + 1;
      m.contacts += this.contacts[i];
      if (!m.ionic) for (let k = 0; k < this.nb[i]; k++) if (this.bc[i * MB + k] === 2) { m.ionic = true; break; }
    }
    const species = new Map();
    const descCache = this._descCache || (this._descCache = new Map());
    mols.forEach(function (m) {
      const sig = JSON.stringify(m.counts) + (m.ionic ? 'i' : '');
      let d = descCache.get(sig);
      if (!d) { d = C.describe(m.counts, m.atoms.length, m.ionic); descCache.set(sig, d); }
      let s = species.get(d.key);
      if (!s) {
        s = { key: d.key, formula: d.formula, name: d.name, molecules: 0, atoms: 0, condensed: 0, color: speciesColor(d.key, m.counts), first: m.atoms[0] };
        species.set(d.key, s);
      }
      s.molecules++;
      s.atoms += m.atoms.length;
      if (m.contacts / m.atoms.length >= 1) s.condensed++;
      m.species = s;
    });
    const molColor = this.molColor;
    mols.forEach(function (m) {
      const col = m.species.color;
      for (let k = 0; k < m.atoms.length; k++) molColor[m.atoms[k]] = col;
    });
    const list = Array.from(species.values());
    list.forEach(function (s) {
      // a "molecule" of a network/lattice counts its formula units
      s.mol = (s.key.slice(-2) === ')n' ? s.atoms / formulaSize(s.key) : s.molecules) * P.MOL_PER_PARTICLE;
      s.phase = s.condensed / s.molecules > 0.5 ? 'condensed' : 'gas';
    });
    list.sort(function (a, b) { return b.mol - a.mol; });
    this.species = list;
    return list;
  };

  function formulaSize(key) {
    const inner = key.slice(1, -2);
    let total = 0;
    inner.replace(/([A-Z][a-z]?)(\d*)/g, function (_, s, d) { total += d ? +d : 1; return ''; });
    return Math.max(total, 1);
  }

  const SPECIES_COLORS = {
    H2O: '#4aa3ff', H2: '#9fb7d0', O2: '#ff6b6b', N2: '#7f9cff', CO2: '#b0b0b0', CO: '#8a8a8a',
    ClNa: '#f4f4f4', '(ClNa)n': '#f4f4f4', CH4: '#ffd36b', H3N: '#c08cff', Cl2: '#9be564',
    HO: '#ffb347', H2O2: '#7fd8ff', O3: '#a0f0ff', ClH: '#c6ff8c', FH: '#d8ff9a', H2S: '#ffe066',
    O2S: '#e0c060', HNaO: '#e6e6ff', Br2: '#b5453a', I2: '#8a2a9a', F2: '#ccf08a'
  };

  function speciesColor(key, counts) {
    if (SPECIES_COLORS[key]) return packColor(SPECIES_COLORS[key]);
    const zs = Object.keys(counts);
    if (zs.length === 1) return elColor[+zs[0]];
    let h = 0;
    for (let k = 0; k < key.length; k++) h = (h * 31 + key.charCodeAt(k)) >>> 0;
    const hue = h % 360, sat = 55 + (h >> 9) % 30, lig = 55 + (h >> 15) % 20;
    return packColor(hslToHex(hue, sat, lig));
  }

  function hslToHex(h, s, l) {
    s /= 100; l /= 100;
    const k = function (n) { return (n + h / 30) % 12; };
    const a = s * Math.min(l, 1 - l);
    const f = function (n) { return l - a * Math.max(-1, Math.min(k(n) - 3, Math.min(9 - k(n), 1))); };
    const hx = function (v) { return ('0' + Math.round(v * 255).toString(16)).slice(-2); };
    return '#' + hx(f(0)) + hx(f(8)) + hx(f(4));
  }

  // Detailed description of the molecule containing atom i (for the inspector).
  S.inspect = function (i) {
    if (i < 0 || i >= this.n) return null;
    const MB = P.MAXB, el = this.el;
    const atoms = this.moleculeOf(i, []);
    const counts = {};
    let ionic = false;
    const bondTypes = new Map();
    for (let h = 0; h < atoms.length; h++) {
      const a = atoms[h];
      counts[el[a]] = (counts[el[a]] || 0) + 1;
      for (let k = 0; k < this.nb[a]; k++) {
        const b = this.bp[a * MB + k];
        if (b < a) continue;
        const cls = this.bc[a * MB + k], order = this.bo[a * MB + k];
        if (cls === 2) ionic = true;
        const s1 = E[el[a]].symbol, s2 = E[el[b]].symbol;
        const pairName = (s1 < s2 ? s1 : s2) + ['', '–', '=', '≡'][order] + (s1 < s2 ? s2 : s1);
        const dEN = Math.abs(en[el[a]] - en[el[b]]);
        const kind = cls === 2 ? 'ionic' : dEN < 0.4 ? 'nonpolar covalent' : 'polar covalent';
        const t = bondTypes.get(pairName) || {
          bond: pairName, count: 0, kind: kind,
          energy: Math.round(C.bondEnergy(el[a], el[b], order)),
          length: Math.round(C.pairLength[(el[a] * NZ + el[b]) * 4 + order] * 10)
        };
        t.count++;
        bondTypes.set(pairName, t);
      }
    }
    const d = C.describe(counts, atoms.length, ionic);
    const z = el[i];
    const geo = this.nb[i] >= 2 ? C.geometry(z, this.nb[i], this.used[i]) : null;
    return {
      atom: E[z],
      index: i,
      formula: d.formula,
      name: d.name,
      atomCount: atoms.length,
      bonds: Array.from(bondTypes.values()),
      geometry: geo,
      charge: this.q[i],
      freeValence: maxVal[z] - this.used[i],
      localT: this.tLoc[i]
    };
  };

  CE.Simulation = Simulation;
})(typeof window !== 'undefined' ? (window.CE = window.CE || {}) : (globalThis.CE = globalThis.CE || {}));
