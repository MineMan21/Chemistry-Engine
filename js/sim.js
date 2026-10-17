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
    DT: 0.001, // substep (s): short enough to follow the vibrations of a strongly bound metal
    // lattice without the integrator leaking energy into it (a held lump would warm by itself)
    DT_REF: 0.0015, // the rates below given "per substep" are for a substep this long, and are
    // rescaled to DT so that what happens per simulated second does not depend on DT
    W: 2400, // jar interior width (world units)
    H: 3200, // jar interior height
    MOL_PER_PARTICLE: 0.001,
    JAR_VOLUME_L: 24, // real volume the jar represents (1 mol of gas at 293 K ~ 1 atm)
    AMBIENT: 293,
    K_REP: 160000, // repulsion stiffness (per unit reduced mass)
    RC_FACTOR: 1.45, // attraction cutoff as a multiple of contact distance
    ION_RC: 75, // Coulomb cutoff (units) between ions (atoms with an ionic bond), about 2.5 Na-Cl
    // contacts: ionic lattices are held by the long reach of the full charges, which the short
    // molecular cutoff would cancel (like charges sit just past the nearest opposite ones)
    ION_KC: 3, // Coulomb between two ions, relative to KC (which is kept gentle for the small partial
    // charges of molecules): enough that salts, oxides and other ionic compounds are solids
    ION_CATION: 0.6, ION_ANION: 1.8, // ionic radii as multiples of the covalent radius: cations
    // shrink and anions swell (Na+ 102 pm, Cl- 181 pm; Mg2+ 72, O2- 140), so ions in a lattice
    // touch at the ionic bond length instead of at the neutral atoms' much larger contact
    KC: 2.0e8, // Coulomb constant (sim energy * units / e^2)
    PARTIAL_CHARGE: 0.3, // partial charge per unit electronegativity difference
    VDW_SCALE: 2.25, // calibration factor for van der Waals well depths (2D, short-ranged well)
    MOL_POLARIZABILITY: 0, // (disabled) per-bond dispersion boost; it made bonding exothermic
    HB_EPS: 800, // hydrogen bond well depth (K) for O-H...O; F stronger, N weaker (water gets more of
    // its cohesion from the many-body W2 term below, so lone pairs in the vapour stay apart)
    HB_ANGLE_POW: 2, // directionality of hydrogen bonds: strength ~ cos^n of the D-H...A bend
    HB_ACCEPT_POW: 2, // ...and of the acceptor side: along the lone pair of a two-bonded O/N/F
    HB_FLAT_LO: 0.6, HB_FLAT_HI: 1.15, // (see HB_FLAT_* note below)
    HH_CONTACT: 0.65, // H...H contact of two O-H/N-H/F-H hydrogens (they overlap in a flattened view)
    W3_EPS: 250, // water network angle stiffness (K), see waterAngle(); gentle, so the liquid stays
    // disordered (ice gets its honeycomb from the crystal engine)
    // cohesion (K per molecule) of water inside liquid or ice: -W2_EPS smoothstep over the number
    // of neighbouring water oxygens from W2_RHO_A to W2_RHO_B. A pair or a small clump in the
    // vapour gets none of it, so steam stays a gas of single molecules instead of drifting clusters
    W2_EPS: 900,
    W2_RHO_A: 1.2,
    W2_RHO_B: 3.2,
    W3_RANGE: 40, // ...its range (units)...
    W3_INNER: 32, // ...full strength out to here, past the hydrogen-bonded O...O distance (~28)
    W2_RANGE: 40, // reach of a water neighbour in that count...
    W2_INNER: 32, // ...and where it starts to fade
    HB_REP_SCALE: 1, // stiffness of the H...acceptor contact relative to other contacts
    HB_CONTACT: 0.78, // hydrogen-bonded H...X pairs may approach to this fraction of contact
    ANNEAL_TRIALS: 0.2, // solid-state annealing: move trials per atom every 4 substeps...
    ANNEAL_ANGLE: 0.5, // ...turning a molecule by up to this many radians...
    ANNEAL_SHIFT: 3, // ...and shifting it by up to this many units
    ICE_A: 31, // O-O spacing of the ice honeycomb: wide enough that ice takes ~7% more room than the liquid (real: 9%)
    CRYSTAL_K: 800, // stiffness (1/s^2) pulling ice members to their lattice sites: gentle, so
    // they only order the crystal; what holds a solid together is the attraction of its atoms
    CRYSTAL_K_METAL: 0, // ...and metal members: none, a metal orders and holds by its own
    // cohesion (springs to a template are not interatomic forces, and leak energy as the
    // template refits, which heats a lump that is only being held)
    CRYSTAL_JOIN_MAX_METAL: 2500, // join strain allowance (K) for metal crystals
    CRYSTAL_TRIES: 0.6, // join/nucleation attempts per atom every 4 substeps
    CRYSTAL_JOIN_MAX: 600, // a unit joins only if its new spring is strained by at most this (K)
    CRYSTAL_HYST: 20, // members melt this far (K) above the melting point
    CRYSTAL_NUC_GLASS: 0.02, // chance per attempt that cold material on the glass starts a crystal
    CRYSTAL_NUC_BULK: 0.0001, // ...and in the bulk
    CRYSTAL_GLASS_RANGE: 45, // "on the glass": this close to a wall or the floor
    CRYSTAL_MAX_GRAINS: 40,
    CRYSTAL_MAX_MEMBERS: 4000,
    // Metallic cohesion has two parts. A pair well, METAL_EPS_PER_MP K per K of melting point,
    // sets the melting point. A many-body embedding, -C smoothstep((rho - RHO_A) / (RHO_B - RHO_A))
    // with rho counting an atom's metal neighbours, rewards being inside dense metal however it is
    // packed: crystal and melt gain the same, so it leaves melting alone, but a lone atom, a small
    // cluster or a froth gains nothing. That is the cohesion that keeps real liquid metals together
    // far above their melting points (the wide liquid range a pair well cannot give, least of all
    // in two dimensions). C makes up the rest of the metal's real cohesive energy (METAL_COHESION,
    // the energy to pull the solid apart into free atoms) beyond what the pair well gives, so a
    // lump holds its atoms as firmly as the real metal: iron 413 kJ/mol, about 27 times its
    // melting temperature, where the pair well alone gives under 7
    METAL_EPS_PER_MP: 2.2, // (calibrated on iron, lead and mercury with nothing else ordering the
    // crystal: the pair well alone decides where the solid melts)
    METAL_EMBED: 1, // (0 switches the many-body term off)
    METAL_RHO_A: 0.2, // rho is about 3.3 inside a solid, 2.9 inside a melt and 2.2 at a flat
    METAL_RHO_B: 1.6, // surface: all of them saturate the embedding, so it cannot favour the
    // crystal (a window reaching the melt's density froze mercury at room temperature); it acts on
    // atoms being pulled out of the metal (a corner atom, one leaving for the vapour, the last
    // threads of a crack)
    METAL_RHO_REACH: 2.2, // reach of the embedding density, in contact distances: beyond the
    // nearest neighbours (as the electron density of a real metal reaches), so it tells the inside
    // of a lump (crystal or melt alike) from its surface, and opening a crack costs real cohesion
    METAL_RC_FACTOR: 1.45, // reach of metal-metal attraction, in contact distances
    GRAB_OMEGA: 12, // grab tool: stiffness of the pull toward the hand (1/s)...
    GRAB_MAX_SPEED: 300, // ...the fastest it moves what it holds (units/s; a few times the thermal
    // speed of its atoms, so what it carries through a gas is not scorched by the gas)...
    GRAB_MAX_ACCEL: 5000, // ...and how hard it accelerates it (units/s^2): a patch of liquid is
    // pulled out of the pool when its cohesion cannot keep up
    METAL_G_REF: 45, // pair well per mass (K/u) at and above which a metal feels full settling force
    G_GAS: 1, // gravity felt by isolated (gas) particles
    G_CONDENSED: 650,
    G_SOLID: 110, // gentler settling for solids, so a lump keeps its shape instead of being
    // crushed flat by its own weight (molecular-scale gravity this strong would exceed cohesion) // extra gravity for particles touching others (liquids/solids settle)
    G_EASE: 0.01, // fraction of the way the settling force on an atom moves to its new value per substep
    SETTLE_MIN: 6, // packed atoms a coarse cell needs before settling gravity starts...
    SETTLE_RANGE: 14, // ...and how many more until it is at full strength
    BULK_VMAX: 900, // max bulk (flow) speed of condensed matter, units/s (falling drops)
    LN_COLL: 18, // ln of real collisions represented by one simulated contact-step (of DT_REF)
    LN_UNI: 16, // ln of the unimolecular attempt frequency (bond dissociation)
    ALPHA: 0.15, // fraction of broken-bond energy needed as activation (Evans-Polanyi like)
    ALPHA_ET: 0.03, ET_EN: 1.1, // same for metals below this electronegativity (electron transfer)
    ALPHA_DISS: 0.1, // same for a diatomic splitting across a metal surface (2 M + O2 -> 2 M-O)
    ALPHA_4C: 0.22, // same for concerted four-centre exchanges (mostly symmetry forbidden)
    NETWORK_O2_BARRIER: 60, // extra activation (kJ/mol) for O2 on an edge of a covalent network
    PASSIVE_BARRIER: 85, // extra activation (kJ/mol) for a halogen splitting on a metal that does not
    // react by electron transfer (iron, copper...): its own halide layer slows it at room temperature
    SPIN_BARRIER: 60, // extra activation (kJ/mol) when closed-shell molecules or metal surfaces
    // attack triplet O=O (why O2 is slow with fuels and metals but fast with radicals)
    HEAT_RADIUS: 90, // reaction heat is shared with the particles within this distance...
    HEAT_GROUP: 24, // ...up to this many (plus their bonded partners)
    GAS_CONDUCTION: 20, // heat flow (1/s) between the gas of neighbouring coarse cells
    HEAT_FRACTION: 0.8, // share of reaction enthalpy deposited as translational kinetic energy
    K_WALL: 110000, // soft wall stiffness (per unit mass)
    K_ANGLE: 8000, // bending spring stiffness (per unit reduced mass); bending is soft
    K_ANGLE_2: 8000, // ...for two-bond centres (water's H-O-H): real ones are stiff
    ACCOUNT_ANGLES: 0, // include bending energy in per-event energy accounting
    PLATE_GAMMA: 80, // Langevin coupling (1/s) of atoms near the hot plate
    PLATE_ZONE: 60, // thickness of the layer above the plate that exchanges heat with it
    WALL_GAMMA: 3, // Langevin coupling to the room through the glass
    WALL_ZONE: 30,
    WALL_ADH_RANGE: 14, // range of the glass's attraction (units beyond contact)
    WALL_ADH_SCALE: 0.6, // glass attraction relative to an atom's own dispersion strength
    WALL_HB: 650, // extra attraction (K) of glass silanol groups for H-bond donors/acceptors
    WALL_FRICTION: 250, // drag (1/s) on liquid touching the glass, so drops cling instead of sliding
    WALL_ACCOM_GAMMA: 0, // heat exchange (1/s) of a lone molecule touching the glass (off: it skews the gauge)
    CONDUCTION: 8, // heat exchange (1/s) of liquids/solids with the plate, at the floor...
    CONDUCTION_DEPTH: 350, // ...fading over this height (units)
    CONVECTION: 0.4, // convective heat exchange (1/s) of the whole jar with the plate, max at bottom
    VALVE_HALF: 140,
    PUMP_RATE: 0.7, // vacuum pump: fraction of the gas removed per simulated second // half-width of the valve opening in the lid
    MAX_STEP: 6, // max displacement per substep (keeps the integrator stable; a limit per substep
    // whatever DT is: clipping bonded atoms harder pumps energy into stiff networks such as graphite)
    TCELL: 160, // size of the coarse cells used for the local temperature field
    T_PRIOR: 2,
    SELF_T_RATE: 0.02, // EMA rate of an atom's own kinetic temperature (per substep)
    SELF_T_WEIGHT: 0.5, // weight of that own temperature in the local temperature // pseudo-count pulling sparse cells toward the jar temperature
    GAS_T_WEIGHT: 0.1, // weight of gas in the local temperature where liquid or solid is present
    SPAWN_GAP: 0.9, // added molecules start at least this fraction of a contact apart
    MIN_DOF: 0.08, // fewest degrees of freedom counted for an atom of a rigid network (of 2)
    NET_DOF: 0.55, // in a covalent network rings and stiff bond angles lock up about half of what
    // the bond count leaves (measured on graphite and silicon held at a known temperature)
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
  // |u1 + u2| of an acceptor's two bond directions below which its lone pairs count as pointing
  // out of the plane (bonds opened past ~145 degrees), and above which they are fully in plane (~110)
  let HB_FLAT_LO = 0.6, HB_FLAT_HI = 1.15;
  HB_ACCEPTOR[7] = 0.7; HB_ACCEPTOR[8] = 1.0; HB_ACCEPTOR[9] = 1.3;

  // Cohesion factors of documented molecular substances, scaling their van der Waals and hydrogen-
  // bond attraction so that each boils near its real boiling point (calibrated against the
  // documented vapour pressures; see README). Substances not listed use 1.
  // Each factor was bisected so that 400 molecules in the insulated jar are half condensed at the
  // temperature where the real substance would be (from its vapour-pressure curve). Water has its
  // own many-body cohesion (W2) and is left at 1.
  const COHESION = {
    N2: 1.82, O2: 2.30, F2: 2.77, Cl2: 2.00, Br2: 0.59, I2: 0.71,
    CO2: 1.66, CH4: 1.19, CO: 0.68, ClH: 0.90, H3N: 0.43, C2H6O: 2.17, C6H6: 2.24, H2O2: 0.93
  };

  const HALOGENS = new Set([9, 17, 35, 53, 85]);

  // Cohesive energies of metals (kJ/mol: the energy to take the solid apart into free atoms; Kittel,
  // Introduction to Solid State Physics). Other metals use 0.12 kJ/mol per K of boiling point,
  // which fits these to within about 20%.
  const METAL_COHESION = {
    Li: 158, Be: 320, Na: 107, Mg: 145, Al: 327, K: 90.1, Ca: 178, Sc: 376, Ti: 468, V: 512, Cr: 395,
    Mn: 282, Fe: 413, Co: 424, Ni: 428, Cu: 336, Zn: 130, Ga: 271, Rb: 82.2, Sr: 166, Y: 422, Zr: 603,
    Nb: 730, Mo: 658, Tc: 661, Ru: 650, Rh: 554, Pd: 376, Ag: 284, Cd: 112, In: 243, Sn: 303, Cs: 77.6,
    Ba: 183, La: 431, Hf: 621, Ta: 782, W: 859, Re: 775, Os: 788, Ir: 670, Pt: 564, Au: 368, Hg: 65,
    Tl: 182, Pb: 196, Bi: 210, Th: 598, U: 536, Pu: 347
  };

  // Per-element constants
  const mass = new Float32Array(NZ), invMass = new Float32Array(NZ), radius = new Float32Array(NZ);
  const en = new Float32Array(NZ), maxVal = new Uint8Array(NZ), normVal = new Uint8Array(NZ);
  const isMetal = new Uint8Array(NZ), vdwEps = new Float32Array(NZ), metalEps = new Float32Array(NZ);
  const elColor = new Uint32Array(NZ), metalA = new Float32Array(NZ), metalG = new Float32Array(NZ);
  const wallAdh = new Float32Array(NZ); // glass attraction well depth (K) per element
  const meltK = new Float32Array(NZ); // melting point (K)
  const pairSigma = new Float32Array(NZ * NZ), pairRc = new Float32Array(NZ * NZ), pairRcIon = new Float32Array(NZ * NZ);
  const pairRhoRc = new Float32Array(NZ * NZ); // reach of the metallic embedding density
  const ionCap = new Uint8Array(NZ); // can form ionic bonds with some element
  const catR = new Float32Array(NZ), anR = new Float32Array(NZ); // ionic contact radii (units)
  const pairEps = new Float32Array(NZ * NZ), pairMetalEps = new Float32Array(NZ * NZ);
  // (re)builds the per-element and pair tables from P; tests and tuning scripts call it again
  // after changing parameters
  function buildTables() {
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
    metalEps[z] = e.isMetal ? P.METAL_EPS_PER_MP * (e.mp || 1000) : 0;
    // cohesive energy (K per atom) minus the pair well's share: half of 6 close-packed neighbours
    const coh = (METAL_COHESION[e.symbol] || 0.12 * (e.bp || 3000)) / P.R_KJ;
    metalA[z] = e.isMetal ? P.METAL_EMBED * Math.max(0, coh - 3 * metalEps[z]) * P.KB : 0;
    // the settling force on a free metal atom, scaled down for heavy, weakly bound metals (it is
    // proportional to mass), so that it never squeezes liquid mercury or gallium into a crystal
    metalG[z] = e.isMetal ? Math.min(1, metalEps[z] / mass[z] / P.METAL_G_REF) : 1;
    elColor[z] = packColor(e.color);
    meltK[z] = e.mp || 1000;
    // metals (mercury!) barely wet glass; everything else feels a fraction of its dispersion
    wallAdh[z] = e.isMetal ? 20 : P.WALL_ADH_SCALE * e.vdwEps * P.VDW_SCALE;
  }
  // Pair tables
  for (let a = 1; a < NZ; a++) {
    for (let b = 1; b < NZ; b++) {
      const p = a * NZ + b;
      pairSigma[p] = radius[a] + radius[b];
      // metallic attraction reaches further: a longer-ranged attraction gives a wide liquid range
      // (a short one sublimes the solid almost as soon as it melts, especially in two dimensions)
      pairRc[p] = pairSigma[p] * (isMetal[a] && isMetal[b] ? P.METAL_RC_FACTOR : P.RC_FACTOR);
      pairRhoRc[p] = isMetal[a] && isMetal[b] ? pairSigma[p] * P.METAL_RHO_REACH : 0;
      pairEps[p] = Math.sqrt(vdwEps[a] * vdwEps[b]) * P.KB;
      pairMetalEps[p] = Math.sqrt(metalEps[a] * metalEps[b]) * P.KB;
    }
  }
  for (let a = 1; a < NZ; a++) { catR[a] = E[a].covRadius / 10 * P.ION_CATION; anR[a] = E[a].covRadius / 10 * P.ION_ANION; }
  for (let a = 1; a < NZ; a++) { ionCap[a] = 0; for (let b = 1; b < NZ; b++) if (C.pairClass[a * NZ + b] === 2) { ionCap[a] = 1; break; } }
  for (let a = 1; a < NZ; a++) for (let b = 1; b < NZ; b++) {
    const p = a * NZ + b;
    pairRcIon[p] = ionCap[a] && ionCap[b] ? Math.max(pairRc[p], P.ION_RC) : 0;
  }
  }
  buildTables();

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
    this.ion = new Uint8Array(M); // 1 for an ion (an atom with an ionic bond): long-range Coulomb
    this.tLoc = new Float32Array(M); // local temperature used for kinetics
    this.tSelf = new Float32Array(M); // each atom's own recent kinetic temperature (EMA)
    this.contacts = new Uint8Array(M); // non-bonded neighbours in contact
    this.condF = new Float32Array(M); // smoothed 'part of a liquid/solid' flag (0..1), sets gravity
    this.molId = new Int32Array(M);
    this.hbd = new Uint8Array(M); // hydrogen-bond donor (H bonded to N, O or F)
    // batch tag of atoms added together as an element: bonds between same-element atoms of one
    // batch just assemble the element's standard state (S8, graphite...), which releases no heat
    this.tag = new Int32Array(M);
    this.anneal = new Uint16Array(M); // substeps of strong thermostatting left after being added
    this.gAcc = new Float32Array(M); // gravity each atom felt in the last substep (units/s^2)
    this.meltT = new Float32Array(M); // melting point of the molecule an atom belongs to (0: unknown)
    this.netF = new Uint8Array(M); // 1: part of a covalent network (graphite, silicon, sulfur chains)
    this.coh = new Float32Array(M).fill(1); // cohesion factor of the molecule (see COHESION)
    this.grainOf = new Int32Array(M).fill(-1); // crystal grain an atom belongs to (-1: none)
    this.grab = null; // the grab tool's hold (see grabStart)
    this.held = new Uint8Array(M); // 1 for atoms the grab tool holds
    this.carried = new Uint8Array(M); // 1 for the solid body hanging from them (see updateCarried)
    this.rho = new Float64Array(M); // metallic neighbour density of free metal atoms
    this.fsP = new Float64Array(M); // minus the derivative of its embedding energy
    this.fsI = new Int32Array(16 * M); this.fsJ = new Int32Array(16 * M);
    this.fsD = new Float64Array(32 * M); // metal pairs inside the density range: i, j, dphi/dr * (dx, dy)/d
    this.mark = new Int32Array(M); this.markTag = 0;
    this.grains = []; // crystal grains, see the crystal engine below
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
    this.pumpRate = P.PUMP_RATE; // vacuum pump strength (fraction of the gas removed per second)
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
  Simulation.buildTables = buildTables;
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
    this.anneal[i] = Math.round(P.ANNEAL_STEPS * P.DT_REF / P.DT);
    this.meltT[i] = 0; this.netF[i] = 0; this.gAcc[i] = 0; this.grainOf[i] = -1; this.coh[i] = 1;
    for (let k = 0; k < P.MAXB; k++) this.bp[i * P.MAXB + k] = -1;
    this.tLoc[i] = T || P.AMBIENT;
    this.tSelf[i] = T || P.AMBIENT;
    this.contacts[i] = 0; this.condF[i] = 0;
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
    // occupancy hash of existing + new atoms (position and radius) to avoid overlapping spawns:
    // two atoms starting inside each other's repulsion would turn into a hot spot
    const cell = 2 * Math.max(r, 26) + 30, occ = new Map();
    const keyOf = function (cx, cy) { return cx * 100003 + cy; };
    const occupy = function (x, y, ra) {
      const k = keyOf(Math.floor(x / cell), Math.floor(y / cell));
      let a = occ.get(k); if (!a) { a = []; occ.set(k, a); } a.push(x, y, ra);
    };
    for (let i = 0; i < this.n; i++) occupy(this.x[i], this.y[i], radius[this.el[i]]);
    // gases keep clear of everything; condensed lattices may pack at their own spacing
    const gap = gas ? 1.1 : 0.95;
    const free1 = function (x, y) {
      const cx = Math.floor(x / cell), cy = Math.floor(y / cell);
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) {
        const arr = occ.get(keyOf(cx + a, cy + b));
        if (!arr) continue;
        for (let m = 0; m < arr.length; m += 3) {
          const dx = arr[m] - x, dy = arr[m + 1] - y, dmin = (r + arr[m + 2]) * gap;
          if (dx * dx + dy * dy < dmin * dmin) return false;
        }
      }
      return true;
    };
    const rng = this.rng;
    const Ld = diatomic ? (C.pairLength[(z * NZ + z) * 4 + C.pairMaxOrder[z * NZ + z]] || 8) : 0;
    let ang = 0;
    const free = function (x, y) {
      if (!diatomic) return free1(x, y);
      ang = rng() * Math.PI * 2;
      const ox = Math.cos(ang) * Ld / 2, oy = Math.sin(ang) * Ld / 2;
      return free1(x - ox, y - oy) && free1(x + ox, y + oy);
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
      if (diatomic) {
        const L = Ld, a = ang;
        occupy(x - Math.cos(a) * L / 2, y - Math.sin(a) * L / 2, r);
        occupy(x + Math.cos(a) * L / 2, y + Math.sin(a) * L / 2, r);
        const i = this.addAtom(z, x - Math.cos(a) * L / 2, y - Math.sin(a) * L / 2, T, tag);
        if (i < 0) break;
        const j = this.addAtom(z, x + Math.cos(a) * L / 2, y + Math.sin(a) * L / 2, T, tag);
        if (j < 0) { this.removeParticle(i); break; }
        this.addBond(i, j, C.pairMaxOrder[z * NZ + z], C.pairClass[z * NZ + z]);
      } else {
        occupy(x, y, r);
        if (this.addAtom(z, x, y, T, tag) < 0) break;
      }
      placed++;
    }
    return { atoms: this.n - added0, units: placed, requested: count };
  };

  // A template drawn by hand (flat, with square angles) can put atoms inside each other's reach:
  // ethanol's hydrogens on one carbon sit on the oxygen of the next. Released as drawn, that
  // strain turns into a flash of a few thousand kelvin. Each template is relaxed once in a private
  // one-molecule jar (forces only, motion removed every step) and the result is reused.
  const relaxCache = new WeakMap();
  function relaxedTemplate(tpl) {
    let r = relaxCache.get(tpl);
    if (r) return r;
    r = tpl;
    if (tpl.atoms.length > 2) {
      const max = P.MAX;
      P.MAX = 64;
      try {
        const s = new Simulation({ rng: function () { return 0.5; } });
        s.react = s.bondEvents = s.crystalUpdate = s.solidAnneal = function () {};
        s.convection = 0; s.insulated = true;
        for (const a of tpl.atoms) s.addAtom(a[0], P.W / 2 + a[1], P.H / 2 + a[2], 1);
        for (const b of tpl.bonds) s.addBond(b[0], b[1], b[2], C.pairClass[tpl.atoms[b[0]][0] * NZ + tpl.atoms[b[1]][0]]);
        for (let k = 0; k < 400; k++) {
          s.step(1);
          for (let i = 0; i < s.n; i++) { s.px[i] = s.x[i]; s.py[i] = s.y[i]; }
        }
        let mx = 0, my = 0;
        for (let i = 0; i < s.n; i++) { mx += s.x[i] / s.n; my += s.y[i] / s.n; }
        if (s.n === tpl.atoms.length) r = { atoms: tpl.atoms.map(function (a, i) { return [a[0], s.x[i] - mx, s.y[i] - my]; }), bonds: tpl.bonds };
      } finally { P.MAX = max; }
    }
    relaxCache.set(tpl, r);
    return r;
  }

  // Adds `count` copies of a small molecule template { atoms: [[z, dx, dy], ...],
  // bonds: [[a, b, order], ...] } scattered around (cx, cy) inside a box of half-size `spread`.
  // Used by the experiment presets (the user otherwise only adds elements).
  S.addMolecules = function (tpl, count, cx, cy, spread, T, spacing) {
    T = T || P.AMBIENT;
    tpl = relaxedTemplate(tpl);
    let size = 0;
    tpl.atoms.forEach(function (a) { size = Math.max(size, Math.hypot(a[1], a[2]) + radius[a[0]] * 0.5); });
    const a0 = size * 2.05 * (spacing || 1);
    let added = 0;
    const tag = ++this.batch;
    // occupancy hash of existing and placed atoms: a long molecule (CO2, ethanol) turned toward its
    // neighbour on the next site would start inside it, and that overlap flashes into thousands of
    // kelvin the moment it is released. Each molecule tries a few turns and takes one that fits.
    const cell = 40, occ = new Map();
    const occupy = function (x, y, r) {
      const k = Math.floor(x / cell) * 100003 + Math.floor(y / cell);
      let a = occ.get(k); if (!a) { a = []; occ.set(k, a); } a.push(x, y, r);
    };
    const fits = function (x, y, r) {
      const cx = Math.floor(x / cell), cy = Math.floor(y / cell);
      for (let a = -1; a <= 1; a++) for (let b = -1; b <= 1; b++) {
        const arr = occ.get((cx + a) * 100003 + cy + b);
        if (!arr) continue;
        for (let m = 0; m < arr.length; m += 3) {
          const dx = arr[m] - x, dy = arr[m + 1] - y, dmin = P.SPAWN_GAP * (r + arr[m + 2]);
          if (dx * dx + dy * dy < dmin * dmin) return false;
        }
      }
      return true;
    };
    for (let i = 0; i < this.n; i++) occupy(this.x[i], this.y[i], radius[this.el[i]]);
    // triangular lattice filling the box (liquid-like packing unless a wider spacing is asked for,
    // e.g. for gases), visited in random order; when the box is full before `count` molecules are
    // in, it grows (the caller asked for an amount, not a box), unless spread.fixed asks for the
    // box as given (a film against the glass)
    for (let pass = 0; pass < (spread.fixed ? 1 : 8) && added < count; pass++) {
      const grow = Math.pow(1.3, pass), sx = spread.x * grow, sy = spread.y * grow, sites = [];
      for (let j = 0; ; j++) {
        const ly = cy - sy + j * a0 * 0.866;
        if (ly > cy + sy) break;
        for (let i = 0; ; i++) {
          const lx = cx - sx + (i + (j & 1) * 0.5) * a0;
          if (lx > cx + sx) break;
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
      for (let k = 0; k < nSites && added < count; k++) {
        const x = sites[2 * k], y = sites[2 * k + 1];
        let ang = 0, ok = false;
        for (let t = 0; t < 8 && !ok; t++) {
          ang = this.rng() * Math.PI * 2;
          const ca = Math.cos(ang), sa = Math.sin(ang);
          ok = tpl.atoms.every(function (a) { return fits(x + a[1] * ca - a[2] * sa, y + a[1] * sa + a[2] * ca, radius[a[0]]); });
        }
        if (!ok) continue;
        const ca = Math.cos(ang), sa = Math.sin(ang);
        tpl.atoms.forEach(function (a) { occupy(x + a[1] * ca - a[2] * sa, y + a[1] * sa + a[2] * ca, radius[a[0]]); });
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
    }
    return added;
  };

  // Copies the molecule containing atom i as a template for addMolecules (positions relative to
  // its centre, bonds with their orders). Returns null for clusters too big to reuse.
  S.captureMolecule = function (i, maxAtoms) {
    if (i < 0 || i >= this.n) return null;
    const atoms = this.moleculeOf(i, []);
    if (atoms.length > (maxAtoms || 120)) return null;
    let cx = 0, cy = 0;
    for (let k = 0; k < atoms.length; k++) { cx += this.x[atoms[k]]; cy += this.y[atoms[k]]; }
    cx /= atoms.length; cy /= atoms.length;
    const index = new Map();
    const tplAtoms = atoms.map(function (a, k) { index.set(a, k); return [this.el[a], +(this.x[a] - cx).toFixed(2), +(this.y[a] - cy).toFixed(2)]; }, this);
    const bonds = [];
    const MB = P.MAXB;
    for (let k = 0; k < atoms.length; k++) {
      const a = atoms[k];
      for (let m = 0; m < this.nb[a]; m++) {
        const b = this.bp[a * MB + m];
        if (b > a) bonds.push([k, index.get(b), this.bo[a * MB + m]]);
      }
    }
    return { atoms: tplAtoms, bonds: bonds };
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
    this.grains = [];
    this.grabEnd();
    this.tCell.fill(this.plateTemp);
    this.cellDebt.fill(0);
    this.present.fill(0);
    this.updateGrid();
  };

  // Removes particle i (breaking its bonds) by swapping the last particle into its slot.
  S.removeParticle = function (i) {
    const MB = P.MAXB;
    if (this.grainOf[i] >= 0) this.crystalEject(i);
    if (this.grab) this.grabRenumber(i, this.n - 1);
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
    this.condF[i] = this.condF[last];
    this.q[i] = this.q[last]; this.tLoc[i] = this.tLoc[last]; this.tSelf[i] = this.tSelf[last]; this.contacts[i] = this.contacts[last];
    this.hbd[i] = this.hbd[last]; this.tag[i] = this.tag[last]; this.anneal[i] = this.anneal[last]; this.meltT[i] = this.meltT[last]; this.netF[i] = this.netF[last]; this.gAcc[i] = this.gAcc[last];
    this.grainOf[i] = this.grainOf[last];
    if (this.grainOf[last] >= 0) this.crystalRenumber(last, i);
    this.grainOf[last] = -1;
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
    // starts at the atoms' current distance (also when they overlap) and eases to its rest
    // length kinematically; snapping it there in one step would fling the atoms apart
    const len = Math.min(Math.max(Math.sqrt(dx * dx + dy * dy), rest * 0.4), rest * 2.5);
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
    // once a bond has been broken the atoms are reaction products, no longer freshly added
    // material assembling itself: from now on re-bonding releases its heat
    if (order < old) { this.tag[i] = 0; this.tag[j] = 0; }
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
    if (this.grainOf[i] >= 0 && !this._crystalBusy) this.crystalEject(i);
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
    let q = 0, ion = 0;
    for (let k = 0; k < this.nb[i]; k++) {
      const ej = this.el[this.bp[base + k]];
      if (this.bc[base + k] === 2) { q += (en[ei] < en[ej] ? 1 : -1) * this.bo[base + k]; ion = 1; }
      else q += P.PARTIAL_CHARGE * (en[ej] - en[ei]);
    }
    this.q[i] = Math.max(-3, Math.min(3, q));
    this.ion[i] = ion;
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

  // Is atom i part of a solid? Metals below their melting point, ionic crystals and covalent
  // networks (chains/rings/lattices of multiply bonded atoms) are; molecular liquids are not.
  S.isSolidAtom = function (i) {
    const z = this.el[i], nb = this.nb[i];
    if (this.grainOf[i] >= 0) return true; // part of a crystal
    if (this.nb[i] === 1 && this.grainOf[this.bp[i * P.MAXB]] >= 0) return true; // ...its H
    // a molecular solid (ice, dry ice, solid iodine...) below its documented melting point
    if (this.meltT[i] && this.tLoc[i] < this.meltT[i] && this.condF[i] > 0.5) return true;
    if (isMetal[z] && nb === 0) return this.tLoc[i] < meltK[z];
    if (nb === 0) return false;
    const MB = P.MAXB, base = i * MB;
    for (let k = 0; k < nb; k++) {
      if (this.bc[base + k] === 2) return true;
      if (nb >= 2 && this.nb[this.bp[base + k]] >= 2) return true;
    }
    return false;
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
        maxRc = Math.max(maxRc, pairRc[a * NZ + b], pairRhoRc[a * NZ + b]);
        if (pairRcIon[a * NZ + b] > maxRc) maxRc = pairRcIon[a * NZ + b];
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
  // Particles move a little between grid rebuilds, so the cells scanned reach a margin beyond r
  // (a missed neighbour would make energy bookkeeping, and Metropolis moves, quietly wrong).
  S.query = function (x, y, r, fn) {
    const cs = this.cellSize, gw = this.gridW, gh = this.gridH, rm = r + 10;
    const x0 = Math.max(0, Math.floor((x - rm) / cs)), x1 = Math.min(gw - 1, Math.floor((x + rm) / cs));
    const y0 = Math.max(0, Math.floor((y - rm) / cs)), y1 = Math.min(gh - 1, Math.floor((y + rm) / cs));
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
    const q = this.q, nb = this.nb, contacts = this.contacts, ghostT = this.ghostT, used = this.used, ion = this.ion, coh = this.coh;
    const bp = this.bp, MB = P.MAXB;
    const gw = this.gridW, gh = this.gridH, start = this.cellStart, items = this.cellItems;
    fx.fill(0, 0, n); fy.fill(0, 0, n);
    contacts.fill(0, 0, n);
    this.nCand = 0;
    const maxCand = this.candI.length;
    const vdwScale = P.VDW_SCALE;
    const molPol = P.MOL_POLARIZABILITY, hbEps = P.HB_EPS * P.KB, hbd = this.hbd;
    const self = this;
    const rho = this.rho, fsI = this.fsI, fsJ = this.fsJ, fsD = this.fsD, fsMax = fsI.length;
    rho.fill(0, 0, n);
    let nFs = 0;

    // atoms three bonds apart (opposite corners of a ring) touch within one molecule: that is not
    // being packed in a liquid, so it must not count as a contact (a lone benzene molecule would
    // otherwise count as condensed and drop like a liquid)
    function oneFour(i, j) {
      const bi = i * MB;
      for (let k = 0; k < nb[i]; k++) {
        const a = bp[bi + k], ba = a * MB;
        for (let l = 0; l < nb[a]; l++) {
          const b = bp[ba + l], bb = b * MB;
          if (b === i) continue;
          for (let m = 0; m < nb[b]; m++) if (bp[bb + m] === j) return true;
        }
      }
      return false;
    }

    function pair(i, j) {
      const dx = x[j] - x[i], dy = y[j] - y[i];
      const d2 = dx * dx + dy * dy;
      const ei = el[i], ej = el[j], p = ei * NZ + ej;
      const rc = pairRc[p], rr = pairRhoRc[p];
      const ionic = ion[i] && ion[j], rcq = ionic ? pairRcIon[p] : rc;
      if (d2 >= rc * rc && (!ionic || d2 >= rcq * rcq) && d2 >= rr * rr) return;
      if (nb[i] && nb[j]) {
        // excluded pairs: bonded (1-2) or sharing a bonded neighbour (1-3); inlined, it is hot
        const bi = i * MB, nbi = nb[i];
        for (let k = 0; k < nbi; k++) {
          const m = bp[bi + k];
          if (m === j) return;
          const bm = m * MB, nbm = nb[m];
          for (let l = 0; l < nbm; l++) if (bp[bm + l] === j) return;
        }
      }
      const d = Math.sqrt(d2) || 1e-6;
      const fm = isMetal[ei] && isMetal[ej] && !nb[i] && !nb[j];
      if (fm && d2 < rr * rr) {
        // metallic density phi(d) = (1 - (d/R)^2)^3 out to its reach R: smooth (twice
        // differentiable everywhere), so the integrator keeps the energy of these strong forces
        const u = 1 - d2 / (rr * rr), phi = u * u * u;
        rho[i] += phi; rho[j] += phi;
        if (nFs < fsMax) {
          const dphi = -6 * u * u / (rr * rr); // (dphi/dd) / d
          fsI[nFs] = i; fsJ[nFs] = j; fsD[2 * nFs] = dphi * dx; fsD[2 * nFs + 1] = dphi * dy; nFs++;
        }
      }
      if (d2 >= rc * rc) {
        // two ions beyond the contact range: long-range Coulomb only
        if (ionic) {
          const F = P.KC * P.ION_KC * q[i] * q[j] * (1 / d2 - 1 / (rcq * rcq));
          const fxp = F * dx / d, fyp = F * dy / d;
          fx[i] -= fxp; fy[i] -= fyp; fx[j] += fxp; fy[j] += fyp;
        }
        return;
      }
      let sigma = ionic ? (q[i] > 0 ? catR[ei] : anR[ei]) + (q[j] > 0 ? catR[ej] : anR[ej]) : pairSigma[p];
      // hydrogen bond: donor H and acceptor N/O/F may approach closer and attract strongly
      const hbs = hbd[i] ? HB_ACCEPTOR[ej] : hbd[j] ? HB_ACCEPTOR[ei] : 0;
      if (hbs) sigma *= P.HB_CONTACT;
      else if (hbd[i] && hbd[j]) sigma *= P.HH_CONTACT;
      if ((ghostT[i] || ghostT[j]) && self.isGhost(i, j)) {
        if (d < sigma) return;
        self.clearGhost(i, j);
      }
      let F = 0; // positive = repulsive
      if (d < sigma) {
        const mu = mass[ei] * mass[ej] / (mass[ei] + mass[ej]);
        F = (hbs ? P.HB_REP_SCALE : 1) * P.K_REP * mu * (sigma - d);
        if (!(nb[i] && nb[j] && oneFour(i, j))) { contacts[i]++; contacts[j]++; }
        if (self.nCand < maxCand && C.pairClass[p]) {
          self.candI[self.nCand] = i; self.candJ[self.nCand] = j; self.nCand++;
        }
      } else {
        const w = rc - sigma, s = (d - sigma) / w;
        let eps;
        if (fm) eps = pairMetalEps[p];
        else {
          const bi = nb[i] > 2 ? 2 : nb[i], bj = nb[j] > 2 ? 2 : nb[j];
          eps = pairEps[p] * vdwScale * (1 + molPol * bi) * (1 + molPol * bj);
          if (coh[i] !== 1 || coh[j] !== 1) eps *= coh[i] === coh[j] ? coh[i] : Math.sqrt(coh[i] * coh[j]);
        }
        F = -4 * eps * s * (1 - s * s) / w;
        if (s < 0.25 && !(nb[i] && nb[j] && oneFour(i, j))) { contacts[i]++; contacts[j]++; }
        // radicals capture partners a little before hard contact (larger reactive cross-section)
        if (s < P.RADICAL_CAPTURE && self.nCand < maxCand && C.pairClass[p] && nb[i] <= 1 && nb[j] <= 1 &&
            (used[i] < maxVal[ei] || used[j] < maxVal[ej])) {
          self.candI[self.nCand] = i; self.candJ[self.nCand] = j; self.nCand++;
        }
      }
      const qq = q[i] * q[j];
      if (qq !== 0) {
        // shifted-force Coulomb, U = KC q q (1/d - 1/rc + (d - rc)/rc^2): energy and force both
        // reach zero smoothly at the cutoff (an abrupt force cutoff slowly heats liquids)
        const dd = d < sigma * 0.8 ? sigma * 0.8 : d;
        F += (ionic ? P.KC * P.ION_KC : P.KC) * qq * (1 / (dd * dd) - 1 / (rcq * rcq));
      }
      const fxp = F * dx / d, fyp = F * dy / d;
      fx[i] -= fxp; fy[i] -= fyp;
      fx[j] += fxp; fy[j] += fyp;
      if (hbs) hbond(i, j, d, sigma, rc, hbs * hbEps * (coh[i] === coh[j] ? coh[i] : Math.sqrt(coh[i] * coh[j])));
    }

    // Directional hydrogen bond D-H...A: strongest when the H points straight at the acceptor,
    // U = -G w(d) a(cos), with w the usual well shape in the H...A distance and a = cos^m of the
    // angle between D->H and H->A. This is what builds open, ordered ice instead of a dense heap.
    const hbPow = P.HB_ANGLE_POW, accPow = P.HB_ACCEPT_POW;
    HB_FLAT_LO = P.HB_FLAT_LO; HB_FLAT_HI = P.HB_FLAT_HI;
    function hbond(i, j, d, sigma, rc, G) {
      const h = hbd[i] && HB_ACCEPTOR[el[j]] ? i : j, A = h === i ? j : i;
      const D = bp[h * MB];
      let ux = x[h] - x[D], uy = y[h] - y[D];
      const L = Math.sqrt(ux * ux + uy * uy) || 1e-6;
      ux /= L; uy /= L;
      const vx = (x[A] - x[h]) / d, vy = (y[A] - y[h]) / d;
      const c = ux * vx + uy * vy;
      if (c <= 0) return;
      const a = Math.pow(c, hbPow), da = hbPow * a / c;
      // acceptor side: an acceptor with two bonds (the O of water, alcohols, ethers) takes the H
      // along its lone-pair direction l, opposite the bisector of its bonds. A molecule seen
      // flattened (its bonds opened wide, as a tilted molecule looks from above) has its lone
      // pairs pointing out of the plane and accepts from any side: lam blends the two.
      let b = 1, q = 0, N1 = -1, N2 = -1, lx = 0, ly = 0, S = 0, u1x = 0, u1y = 0, u2x = 0, u2y = 0, L1 = 1, L2 = 1;
      let dbq = 0, dbS = 0;
      if (nb[A] === 2 && accPow) {
        N1 = bp[A * MB]; N2 = bp[A * MB + 1];
        u1x = x[N1] - x[A]; u1y = y[N1] - y[A]; L1 = Math.sqrt(u1x * u1x + u1y * u1y) || 1e-6; u1x /= L1; u1y /= L1;
        u2x = x[N2] - x[A]; u2y = y[N2] - y[A]; L2 = Math.sqrt(u2x * u2x + u2y * u2y) || 1e-6; u2x /= L2; u2y /= L2;
        const sx = u1x + u2x, sy = u1y + u2y;
        S = Math.sqrt(sx * sx + sy * sy);
        const t = (S - HB_FLAT_LO) / (HB_FLAT_HI - HB_FLAT_LO);
        if (t <= 0) N1 = -1; // fully flattened: an isotropic acceptor
        else {
          const tc = t >= 1 ? 1 : t, lam = tc * tc * (3 - 2 * tc);
          const dlam = t >= 1 ? 0 : 6 * tc * (1 - tc) / (HB_FLAT_HI - HB_FLAT_LO);
          lx = -sx / S; ly = -sy / S;
          q = -(lx * vx + ly * vy); // cos between l and A->H (= -v)
          const qn = q > 0 ? Math.pow(q, accPow) : 0;
          b = 1 - lam + lam * qn;
          if (b <= 0) return;
          dbq = q > 0 ? lam * accPow * qn / q : 0;
          dbS = dlam * (qn - 1);
        }
      }
      let w = 1, dw = 0;
      if (d > sigma) {
        const W = rc - sigma, s = (d - sigma) / W;
        w = (1 - s * s) * (1 - s * s); dw = -4 * s * (1 - s * s) / W;
      }
      // U = -G w a b; F_X = G (a b dw dd/dX + w b da dc/dX + w a db dq/dX)
      const ka = G * a * b * dw, kc = G * w * b * da;
      const ax = (ux - c * vx) / d, ay = (uy - c * vy) / d; // dc/dA
      const bx = (vx - c * ux) / L, by = (vy - c * uy) / L; // -dc/dD
      fx[A] += ka * vx + kc * ax; fy[A] += ka * vy + kc * ay;
      fx[h] += -ka * vx + kc * (bx - ax); fy[h] += -ka * vy + kc * (by - ay);
      fx[D] -= kc * bx; fy[D] -= kc * by;
      if (N1 >= 0) {
        const kq = G * w * a * dbq, kS = G * w * a * dbS;
        // with e = A->H unit vector (= -v): dq/dH = (l - q e)/d, dq/dA gets the opposite of that
        const hx = kq * (lx + q * vx) / d, hy = kq * (ly + q * vy) / d;
        // through the bond directions s = u1 + u2: dq/ds = -(e - q l)/S and dS/ds = s/S = -l;
        // d/dN = (g - (u.g)u)/L for a gradient g with respect to s
        const gx = kq * (vx + q * lx) / S - kS * lx, gy = kq * (vy + q * ly) / S - kS * ly;
        const g1 = u1x * gx + u1y * gy, g2 = u2x * gx + u2y * gy;
        const n1x = (gx - g1 * u1x) / L1, n1y = (gy - g1 * u1y) / L1;
        const n2x = (gx - g2 * u2x) / L2, n2y = (gy - g2 * u2y) / L2;
        fx[h] += hx; fy[h] += hy;
        fx[A] += -hx - n1x - n2x; fy[A] += -hy - n1y - n2y;
        fx[N1] += n1x; fy[N1] += n1y;
        fx[N2] += n2x; fy[N2] += n2y;
      }
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
    if (nFs) {
      // embedding forces: U = sum F_i(rho_i); on the pair i, j the force is
      // -(F'_i + F'_j) dphi/dr along the pair
      const fsP = this.fsP;
      for (let i = 0; i < n; i++) fsP[i] = rho[i] > 0 ? -embedD(el[i], rho[i]) : 0;
      for (let k = 0; k < nFs; k++) {
        const i = fsI[k], j = fsJ[k], c = fsP[i] + fsP[j];
        const gx = c * fsD[2 * k], gy = c * fsD[2 * k + 1];
        fx[j] += gx; fy[j] += gy;
        fx[i] -= gx; fy[i] -= gy;
      }
    }
    if (P.W3_EPS) for (let i = 0; i < n; i++) if (this.isWaterO(i)) this.waterAngle(i, true);
    if (P.W2_EPS) this.waterDenseForces();
    if (this.grains.length) this.crystalForces();
  };

  // ---------------------------------------------------------------------------------------
  // Integration, constraints, walls
  // ---------------------------------------------------------------------------------------
  S.integrate = function () {
    const n = this.n, x = this.x, y = this.y, px = this.px, py = this.py, el = this.el;
    const dt = P.DT, dt2 = dt * dt, sub = dt / P.DT_REF, maxS = P.MAX_STEP;
    const condRate = perSub(0.04), gEase = perSub(P.G_EASE), selfRate = perSub(P.SELF_T_RATE);
    const contacts = this.contacts;
    const tgw = this.tgw, tgh = this.tgh, tSum = this.tSum, tCnt = this.tCnt, tCell = this.tCell;
    // condensed matter and gas are summed apart: where both share a cell, the few gas molecules
    // (which hold very little heat next to a liquid or solid) barely move the local temperature
    // of the condensed matter; a lump swung through vapour does not melt from the vapour it stirs
    const tSumG = this._tSumG || (this._tSumG = new Float64Array(tSum.length));
    const tCntG = this._tCntG || (this._tCntG = new Float32Array(tSum.length));
    tSum.fill(0); tCnt.fill(0); tSumG.fill(0); tCntG.fill(0);
    const cellKE = this.cellKE, cellDebt = this.cellDebt;
    cellKE.fill(0);
    const W = P.W, H = P.H, rng = this.rng, wallG = this.insulated ? 0 : P.WALL_GAMMA;
    const grab = this.grab, held = grab ? this.held : null;
    let ke = 0, dofHalf = 0, impulse = 0;
    // bulk (flow) velocity of condensed matter per coarse cell
    const mom = this.cellMom, TC = P.TCELL;
    mom.fill(0);
    // number of packed atoms per coarse cell: only real liquids/solids settle, not a pair of gas
    // molecules that briefly stick together (strong gravity switching on and off pumps heat in)
    const cndN = this.cellCnd || (this.cellCnd = new Uint16Array(tgw * tgh));
    cndN.fill(0);
    // a molecule settles as a whole: an atom counts as condensed if it or a bonded neighbour
    // touches something; smoothed in time so gravity does not flicker (which would heat liquids)
    const condF = this.condF, MB = P.MAXB;
    for (let i = 0; i < n; i++) {
      // (an atom of a covalent network is condensed matter even where its only neighbours are the
      // ones it is bonded to: otherwise a lump of graphite or silicon floats as if it were a gas)
      let c = contacts[i] >= 2 || this.netF[i] ? 1 : 0;
      if (!c) for (let k = 0; k < this.nb[i]; k++) if (contacts[this.bp[i * MB + k]] >= 2) { c = 1; break; }
      // changes gradually both ways: a force that flickers with contacts would pump energy into
      // surface molecules
      condF[i] += condRate * (c - condF[i]);
    }
    for (let i = 0; i < n; i++) {
      // a held body moves with the hand: it is not part of the flow of what it passes through
      if (contacts[i] < 2 || (held && held[i])) continue;
      let tcx = Math.floor(x[i] / TC), tcy = Math.floor(y[i] / TC);
      tcx = tcx < 0 ? 0 : tcx >= tgw ? tgw - 1 : tcx;
      tcy = tcy < 0 ? 0 : tcy >= tgh ? tgh - 1 : tcy;
      const c3 = (tcy * tgw + tcx) * 3, m = mass[el[i]];
      mom[c3] += m * (x[i] - px[i]); mom[c3 + 1] += m * (y[i] - py[i]); mom[c3 + 2] += m;
      cndN[tcy * tgw + tcx]++;
    }
    for (let i = 0; i < n; i++) {
      let vx = x[i] - px[i], vy = y[i] - py[i];
      const v2 = vx * vx + vy * vy;
      if (v2 > maxS * maxS) { const s = maxS / Math.sqrt(v2); vx *= s; vy *= s; }
      const im = invMass[el[i]];
      const c = contacts[i] > 4 ? 4 : contacts[i];
      // settling gravity only for material that is clearly packed (passing gas collisions give
      // a small condF that must not pile the vapour up near the floor)
      let cs = (condF[i] - 0.5) * 3.3;
      if (cs > 0) {
        let tcx = Math.floor(x[i] / TC), tcy = Math.floor(y[i] / TC);
        tcx = tcx < 0 ? 0 : tcx >= tgw ? tgw - 1 : tcx;
        tcy = tcy < 0 ? 0 : tcy >= tgh ? tgh - 1 : tcy;
        const nc = (cndN[tcy * tgw + tcx] - P.SETTLE_MIN) / P.SETTLE_RANGE;
        cs *= nc <= 0 ? 0 : nc >= 1 ? 1 : nc;
      }
      let g = (this.isSolidAtom(i) ? P.G_SOLID : P.G_CONDENSED) * (cs <= 0 ? 0 : cs >= 1 ? 1 : cs);
      if (this.nb[i] === 0) g *= metalG[el[i]];
      // the hand carries the weight of what it holds: the patch it grips and the solid body
      // hanging from it (the settling force is far stronger than real gravity on this scale;
      // its weight would tear a held lump apart, which real weight never does to a small one)
      if (held && (held[i] || this.carried[i])) g = 0;
      g += P.G_GAS;
      // eased in time: a settling force that flickers as material moves between crowded and empty
      // regions (a lump carried across the jar) would pump energy into it
      g = this.gAcc[i] + gEase * (g - this.gAcc[i]);
      this.gAcc[i] = g;
      // cap the local bulk flow speed (not thermal motion) of condensed matter
      let bvx = 0, bvy = 0, noise = -1;
      if (held && held[i]) {
        // a held body moves as one with the hand: that is its flow, under the same speed limits
        bvx = grab.vx; bvy = grab.vy; noise = 0;
      } else if (c >= 2 || condF[i] > 0.5) {
        // (packed matter, its surface included: an edge atom of a moving lump moves with it)
        let tcx = Math.floor(x[i] / TC), tcy = Math.floor(y[i] / TC);
        tcx = tcx < 0 ? 0 : tcx >= tgw ? tgw - 1 : tcx;
        tcy = tcy < 0 ? 0 : tcy >= tgh ? tgh - 1 : tcy;
        const c3 = (tcy * tgw + tcx) * 3;
        if (mom[c3 + 2] > 0) {
          bvx = mom[c3] / mom[c3 + 2]; bvy = mom[c3 + 1] / mom[c3 + 2];
          // never act on the cell's random thermal drift (~ sqrt(kT / M_cell))
          noise = 4 * Math.sqrt(P.KB * this.tLoc[i] / mom[c3 + 2]);
        }
      }
      if (noise >= 0) {
        const vmaxB = Math.max(noise, P.BULK_VMAX) * dt;
        const vb = Math.sqrt(bvx * bvx + bvy * bvy);
        if (vb > vmaxB) { const f = 1 - vmaxB / vb; vx -= bvx * f; vy -= bvy * f; bvx *= 1 - f; bvy *= 1 - f; }
      }
      // soft walls (energy conserving with bond constraints, unlike reflections)
      const rw = radius[el[i]] * 0.6;
      let wfx = 0, wfy = 0, Tw = -1, gamma = 0;
      const WZ = P.WALL_ZONE;
      if (x[i] < WZ || x[i] > W - WZ || y[i] < WZ) { Tw = P.AMBIENT; gamma = wallG; }
      if (x[i] < rw) wfx = P.K_WALL * (rw - x[i]);
      else if (x[i] > W - rw) wfx = -P.K_WALL * (x[i] - W + rw);
      // pressure gauge: momentum pushed into the upper half of the side walls (contact only)
      if (wfx !== 0 && y[i] < H / 2) impulse += Math.abs(wfx) * mass[el[i]] * dt;
      // glass attracts atoms near the side walls and floor (polar molecules such as water wet
      // it, which makes the meniscus and lets drops cling); conservative, like the pair wells
      const adhR = P.WALL_ADH_RANGE;
      let adhered = false;
      if (x[i] < rw + adhR || x[i] > W - rw - adhR || y[i] > H - rw - adhR) {
        const z = el[i];
        // glass is a polar, ionic surface: hydrogen bonders and ions (salt deposits) cling to it
        const eA = (wallAdh[z] + ((this.hbd[i] || (HB_ACCEPTOR[z] && this.nb[i]) || this.ion[i]) ? P.WALL_HB : 0)) * P.KB * im;
        let d = x[i] - rw, pull = 0;
        if (d >= 0 && d < adhR) { const s = d / adhR; pull = 4 * eA * s * (1 - s * s) / adhR; wfx -= pull; adhered = true; }
        d = W - rw - x[i];
        if (d >= 0 && d < adhR) { const s = d / adhR; pull = 4 * eA * s * (1 - s * s) / adhR; wfx += pull; adhered = true; }
        // the glass is pulled inward by what it attracts: subtract it from the gauge reading
        if (pull && y[i] < H / 2) impulse -= pull * mass[z] * dt;
        d = H - rw - y[i];
        if (d >= 0 && d < adhR) { const s = d / adhR; wfy += 4 * eA * s * (1 - s * s) / adhR; }
      }
      if (y[i] < rw) {
        wfy = P.K_WALL * (rw - y[i]);
        impulse += wfy * mass[el[i]] * dt;
      } else if (y[i] > H - P.PLATE_ZONE) {
        if (y[i] > H - rw) wfy = -P.K_WALL * (y[i] - H + rw);
        Tw = this.plateTemp; gamma = P.PLATE_GAMMA;
      }

      if (this.anneal[i] && (--this.anneal[i], contacts[i] >= 1 || this.nb[i] >= 2)) {
        // freshly added condensed material: remove the heat of settling / self-assembly by
        // rescaling only the random motion (relative to the local flow), toward room temperature,
        // so a freshly dropped blob still falls freely
        let network = false;
        if (this.nb[i] >= 2) for (let k = 0; k < this.nb[i]; k++) if (this.nb[this.bp[i * MB + k]] >= 2) { network = true; break; }
        if (network) {
          // covalent networks (carbon, sulfur...): a plain strong thermostat while they assemble
          Tw = P.AMBIENT; gamma = 40;
        } else {
          const tc = tCell[this.cellOfT[i]] > 1 ? tCell[this.cellOfT[i]] : 1;
          let f = Math.sqrt(1 + 0.06 * (P.AMBIENT / tc - 1));
          f = f < 0.9 ? 0.9 : f > 1.03 ? 1.03 : f;
          vx = bvx + (vx - bvx) * f; vy = bvy + (vy - bvy) * f;
        }
      } else if (adhered) {
        // touching the glass: heat exchange with the room-temperature glass (so vapour condenses
        // on the walls of a hot jar, and lone adsorbed molecules do not rattle hot in the well),
        // plus a firm grip on liquid so drops cling
        Tw = this.insulated ? this.tLoc[i] : P.AMBIENT;
        gamma = contacts[i] >= 1 ? P.WALL_FRICTION : P.WALL_ACCOM_GAMMA;
      } else if (gamma === 0 && this.convection > 0) {
        // convection: everything exchanges heat with the plate, more strongly near the bottom;
        // liquids and solids also conduct the plate's heat upward through themselves
        const h = y[i] / H;
        Tw = this.plateTemp; gamma = this.convection * (0.4 + 0.6 * h * h);
        if (condF[i] > 0.5) gamma += P.CONDUCTION * Math.exp(-(H - y[i]) / P.CONDUCTION_DEPTH);
      }
      if (gamma > 0) {
        // Langevin (Ornstein-Uhlenbeck) thermostat towards the wall temperature
        const c1 = Math.exp(-gamma * dt);
        const sd = Math.sqrt((1 - c1 * c1) * P.KB * Tw * im) * dt;
        // cheap unit-variance noise (sum of 3 uniforms); exact Gaussians are not needed here
        // a held body is thermostatted about its own motion: heat exchange is not friction, so
        // the plate, glass or air does not anchor what the hand is lifting (dragging the rest
        // of a lump against the patch in the hand would strain and heat it)
        const hvx = held ? (held[i] ? grab.vx : this.carried[i] ? grab.bvx : 0) : 0;
        const hvy = held ? (held[i] ? grab.vy : this.carried[i] ? grab.bvy : 0) : 0;
        vx = hvx + (vx - hvx) * c1 + sd * 2 * (rng() + rng() + rng() - 1.5);
        vy = hvy + (vy - hvy) * c1 + sd * 2 * (rng() + rng() + rng() - 1.5);
      }
      px[i] = x[i]; py[i] = y[i];
      x[i] += vx + (this.fx[i] * im + wfx) * dt2;
      y[i] += vy + (this.fy[i] * im + g + wfy) * dt2;
      // local temperature field: accumulate kinetic energy per coarse cell. Each rigid bond
      // removes one degree of freedom shared by its two atoms, so an atom keeps 2 - nb/2 of its
      // 2 (equipartition); a covalent network keeps very few (graphite 0.5 per atom), and the
      // temperature is its kinetic energy over the degrees of freedom actually left
      // (motion relative to the local flow: a held body's own motion is the hand's, not heat)
      const rvx = vx - bvx, rvy = vy - bvy;
      const tInst = mass[el[i]] * (rvx * rvx + rvy * rvy) / (dt2 * 2 * P.KB);
      // (an atom's own reading shares the bonds' constraints out more evenly: a carbon with four
      // bonds still moves, the count only balances over the whole molecule or region)
      const dof = Math.max(P.MIN_DOF, (2 - 0.5 * this.nb[i]) * (this.netF[i] ? P.NET_DOF : 1));
      const dofScale = 2 / Math.max(1, 2 - 0.5 * this.nb[i]);
      let tcx = Math.floor(x[i] / P.TCELL), tcy = Math.floor(y[i] / P.TCELL);
      tcx = tcx < 0 ? 0 : tcx >= tgw ? tgw - 1 : tcx;
      tcy = tcy < 0 ? 0 : tcy >= tgh ? tgh - 1 : tcy;
      const tc = tcy * tgw + tcx;
      this.cellOfT[i] = tc;
      if (condF[i] <= 0.5 && !(held && held[i])) { tSumG[tc] += tInst; tCntG[tc] += 0.5 * dof; }
      else { tSum[tc] += tInst; tCnt[tc] += 0.5 * dof; }
      this.tSelf[i] += selfRate * (tInst * dofScale - this.tSelf[i]);
      cellKE[tc] += tInst * P.KB;
      ke += tInst;
      dofHalf += 0.5 * dof;
    }
    // T = sum(KE) / (KB * DOF / 2) with DOF = 2 per atom minus one per bond constraint
    const globalT = n ? ke / dofHalf : this.plateTemp;
    // local temperature: 3x3 blurred cell average, shrunk toward the global temperature when
    // few particles are present (one fast gas atom is not a hot spot), smoothed in time
    const N0 = P.T_PRIOR;
    const tNew = this._tNew || (this._tNew = new Float32Array(tCell.length));
    for (let cy = 0; cy < tgh; cy++) {
      for (let cx = 0; cx < tgw; cx++) {
        let sum = 0, cnt = 0, sumG = 0, cntG = 0;
        for (let oy = -1; oy <= 1; oy++) {
          const yy = cy + oy;
          if (yy < 0 || yy >= tgh) continue;
          for (let ox = -1; ox <= 1; ox++) {
            const xx = cx + ox;
            if (xx < 0 || xx >= tgw) continue;
            const w = (ox === 0 && oy === 0) ? 1 : 0.35, c = yy * tgw + xx;
            sum += tSum[c] * w; cnt += tCnt[c] * w; sumG += tSumG[c] * w; cntG += tCntG[c] * w;
          }
        }
        const wg = cnt > 3 ? P.GAS_T_WEIGHT : 1;
        tNew[cy * tgw + cx] = (sum + wg * sumG + N0 * globalT) / (cnt + wg * cntG + N0);
      }
    }
    const tRate = perSub(0.3);
    for (let c = 0; c < tCell.length; c++) tCell[c] += tRate * (tNew[c] - tCell[c]);
    // kinetics see a blend of the regional field and the atom's own recent temperature, so a
    // small hot spot (a flame front in a dilute gas) is not averaged away
    const ws = P.SELF_T_WEIGHT;
    for (let i = 0; i < n; i++) {
      // only in the gas phase: in dense matter the regional field is already accurate
      // (and only lone atoms and molecule ends: a bonded centre's own reading is too rough)
      const w = contacts[i] >= 2 || this.nb[i] >= 2 ? 0 : ws;
      const tc = tCell[this.cellOfT[i]];
      // at a liquid's surface one fast-jiggling molecule is not a fire; in a gas flame it may be
      const own = condF[i] > 0.3 && this.tSelf[i] > 3 * tc ? 3 * tc : this.tSelf[i];
      this.tLoc[i] = (1 - w) * tc + w * own;
    }
    // repay local energy debt by gently cooling the particles of indebted cells
    let anyDebt = false;
    for (let c = 0; c < cellDebt.length; c++) if (cellDebt[c] > 0) { anyDebt = true; break; }
    if (anyDebt) {
      const fac = this._fac || (this._fac = new Float32Array(cellDebt.length));
      for (let c = 0; c < cellDebt.length; c++) {
        fac[c] = 1;
        if (cellDebt[c] <= 0 || cellKE[c] <= 0) continue;
        const pay = Math.min(cellDebt[c], cellKE[c] * 0.2 * sub);
        cellDebt[c] -= pay;
        fac[c] = Math.sqrt(1 - pay / cellKE[c]);
      }
      for (let i = 0; i < n; i++) {
        const f = fac[this.cellOfT[i]];
        if (f === 1 || (held && held[i])) continue;
        // (what the hand carries is cooled about its own motion, not slowed down)
        const ux = held && this.carried[i] ? grab.bvx : 0, uy = held && this.carried[i] ? grab.bvy : 0;
        px[i] = x[i] - ux - (x[i] - px[i] - ux) * f; py[i] = y[i] - uy - (y[i] - py[i] - uy) * f;
      }
    }
    this.impulse += impulse;
    this._keSum = (this._keSum || 0) + ke;
    this._keCount = (this._keCount || 0) + dofHalf;
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

  // Heat conduction through the gas between neighbouring coarse cells. Simulated gas molecules
  // collide only a few times per simulated second, far too rarely to carry heat (a flame front)
  // across the jar the way a real gas does, so energy also flows from hot gas to cooler gas
  // next to it. Each cell's gas atoms are rescaled together: energy is conserved exactly.
  S.gasConduction = function () {
    const kappa = P.GAS_CONDUCTION;
    if (!kappa) return;
    const n = this.n, x = this.x, y = this.y, px = this.px, py = this.py, el = this.el;
    const tgw = this.tgw, tgh = this.tgh, nc = tgw * tgh, condF = this.condF, cellOf = this.cellOfT;
    if (!this._gKE || this._gKE.length !== nc) {
      this._gKE = new Float64Array(nc); this._gN = new Float64Array(nc); this._gdE = new Float64Array(nc);
    }
    const gKE = this._gKE, gN = this._gN, gdE = this._gdE;
    gKE.fill(0); gN.fill(0); gdE.fill(0);
    const idt2 = 1 / (P.DT * P.DT);
    // (a held body is not gas: its motion is the hand's, not heat to pass on)
    const held = this.grab ? this.held : null;
    for (let i = 0; i < n; i++) {
      if (condF[i] > 0.5 || (held && held[i])) continue;
      const vx = x[i] - px[i], vy = y[i] - py[i], c = cellOf[i];
      gKE[c] += 0.5 * mass[el[i]] * (vx * vx + vy * vy) * idt2;
      gN[c]++;
    }
    const f = kappa * P.DT;
    for (let cy = 0; cy < tgh; cy++) {
      for (let cx = 0; cx < tgw; cx++) {
        const a = cy * tgw + cx;
        if (gN[a] < 1) continue;
        const Ta = gKE[a] / gN[a];
        for (let d = 0; d < 2; d++) {
          const b = d === 0 ? (cx + 1 < tgw ? a + 1 : -1) : (cy + 1 < tgh ? a + tgw : -1);
          if (b < 0 || gN[b] < 1) continue;
          // flow ~ temperature difference times the smaller heat capacity
          const flow = f * (Ta - gKE[b] / gN[b]) * Math.min(gN[a], gN[b]);
          gdE[a] -= flow; gdE[b] += flow;
        }
      }
    }
    for (let c = 0; c < nc; c++) {
      if (gN[c] < 1 || gdE[c] === 0) { gdE[c] = 1; continue; }
      const ke = gKE[c], e = ke + gdE[c];
      gdE[c] = ke > 0 && e > 0 ? Math.sqrt(e / ke) : 1; // reuse as the velocity scale
    }
    for (let i = 0; i < n; i++) {
      if (condF[i] > 0.5 || (held && held[i])) continue;
      const s = gdE[cellOf[i]];
      if (s === 1) continue;
      px[i] = x[i] - (x[i] - px[i]) * s;
      py[i] = y[i] - (y[i] - py[i]) * s;
    }
  };

  // Moves each bond's target length toward its rest length, then re-satisfies all bonds with a
  // kinematic constraint pass (positions and previous positions move together), so the
  // geometry change of a new/changed bond never turns into velocity - even in networks where
  // many bonds ramp at once.
  S.rampBonds = function () {
    const n = this.n, MB = P.MAXB, el = this.el, bl = this.bl, ramp = P.BOND_RAMP * P.DT / P.DT_REF;
    let any = false;
    for (let i = 0; i < n; i++) {
      const base = i * MB;
      for (let k = 0; k < this.nb[i]; k++) {
        const rest = C.pairLength[(el[i] * NZ + el[this.bp[base + k]]) * 4 + this.bo[base + k]];
        const cur = bl[base + k];
        if (cur === rest) continue;
        bl[base + k] = Math.abs(cur - rest) <= ramp ? rest : cur - Math.sign(cur - rest) * ramp;
        any = true;
      }
    }
    if (!any) return;
    const x = this.x, y = this.y, px = this.px, py = this.py;
    for (let it = 0; it < P.SHAKE_ITERS; it++) {
      for (let i = 0; i < n; i++) {
        const base = i * MB;
        for (let k = 0; k < this.nb[i]; k++) {
          const j = this.bp[base + k];
          if (j < i) continue;
          const dx = x[j] - x[i], dy = y[j] - y[i];
          const d = Math.sqrt(dx * dx + dy * dy) || 1e-6;
          const diff = d - bl[base + k];
          if (diff > -1e-4 && diff < 1e-4) continue;
          const wi = invMass[el[i]], wj = invMass[el[j]];
          const sh = diff / (d * (wi + wj));
          const sxi = dx * sh * wi, syi = dy * sh * wi, sxj = dx * sh * wj, syj = dy * sh * wj;
          x[i] += sxi; y[i] += syi; px[i] += sxi; py[i] += syi;
          x[j] -= sxj; y[j] -= syj; px[j] -= sxj; py[j] -= syj;
        }
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
    function spring(a, b, target, k) {
      const dx = x[b] - x[a], dy = y[b] - y[a];
      const d = Math.sqrt(dx * dx + dy * dy) || 1e-6;
      const ma = mass[el[a]], mb = mass[el[b]];
      const F = (k || K) * (ma * mb / (ma + mb)) * (d - target) / d;
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
        spring(this.bp[base], this.bp[base + 1], Math.sqrt(la * la + lb * lb - 2 * la * lb * Math.cos(theta)), P.K_ANGLE_2);
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

  // Same as heat() with the energy already in simulation units. The energy is shared with the
  // molecules around the reacting pair (the third bodies that carry off the heat of a real
  // recombination): the two atoms alone cannot hold it, because a new rigid bond leaves them
  // only rotation, and the constraint would throw most of it away.
  S.heatRaw = function (i, j, e) {
    const dt = P.DT, x = this.x, y = this.y, px = this.px, py = this.py, MB = P.MAXB;
    const grp = this._grp || (this._grp = new Int32Array(P.HEAT_GROUP + 2 * MB + 2));
    let g = 0;
    grp[g++] = i; grp[g++] = j;
    const self = this;
    if (this.cellItems) {
      this.query(0.5 * (x[i] + x[j]), 0.5 * (y[i] + y[j]), P.HEAT_RADIUS, function (k) {
        if (g >= P.HEAT_GROUP || k === i || k === j) return;
        for (let a = 0; a < g; a++) if (grp[a] === k) return;
        grp[g++] = k;
      });
    }
    // whole small molecules move together, so scaling their velocities keeps bonds rigid
    for (let a = 0, g0 = g; a < g0; a++) {
      const k = grp[a];
      for (let b = 0; b < this.nb[k] && g < grp.length; b++) {
        const m = this.bp[k * MB + b];
        let seen = false;
        for (let c = 0; c < g; c++) if (grp[c] === m) { seen = true; break; }
        if (!seen) grp[g++] = m;
      }
    }
    // a new bond i-j holds the pair at fixed distance: its closing speed along the bond would be
    // removed by the constraint, so turn it into heat here instead of losing it
    if (this.isBonded(i, j)) {
      const mi = mass[this.el[i]], mj = mass[this.el[j]], mu = mi * mj / (mi + mj);
      const dx = x[j] - x[i], dy = y[j] - y[i];
      const d = Math.sqrt(dx * dx + dy * dy) || 1e-6, ux = dx / d, uy = dy / d;
      const vr = ((x[j] - px[j]) - (x[i] - px[i])) * ux + ((y[j] - py[j]) - (y[i] - py[i])) * uy;
      // remove the radial relative motion (momentum conserving) and bank its energy
      px[i] -= vr * ux * mj / (mi + mj); py[i] -= vr * uy * mj / (mi + mj);
      px[j] += vr * ux * mi / (mi + mj); py[j] += vr * uy * mi / (mi + mj);
      e += 0.5 * mu * vr * vr / (dt * dt);
    }
    let M = 0, cx = 0, cy = 0;
    for (let a = 0; a < g; a++) {
      const k = grp[a], m = mass[this.el[k]];
      M += m; cx += m * (x[k] - px[k]); cy += m * (y[k] - py[k]);
    }
    cx /= M; cy /= M;
    let K = 0;
    for (let a = 0; a < g; a++) {
      const k = grp[a], rx = x[k] - px[k] - cx, ry = y[k] - py[k] - cy;
      K += 0.5 * mass[this.el[k]] * (rx * rx + ry * ry);
    }
    K /= dt * dt;
    let s;
    if (e >= 0) {
      if (K < 1e-9) return this.heatPair(i, j, e);
      s = Math.sqrt((K + e) / K);
    } else {
      // take up to half of what the group has; the rest becomes a local energy debt that is
      // drained from the cell's particles, so energy is conserved overall
      const take = Math.min(-e, K * 0.5);
      this.cellDebt[this.cellOfT[i]] += -e - take;
      s = Math.sqrt((K - take) / Math.max(K, 1e-9));
    }
    for (let a = 0; a < g; a++) {
      const k = grp[a];
      px[k] = x[k] - (cx + (x[k] - px[k] - cx) * s);
      py[k] = y[k] - (cy + (y[k] - py[k] - cy) * s);
      if (e > 0) this.eRef[k] += e / g;
    }
  };

  // Fallback: deposit into the relative motion of the pair alone (a group with no motion).
  S.heatPair = function (i, j, e) {
    const dt = P.DT, x = this.x, y = this.y, px = this.px, py = this.py;
    const mi = mass[this.el[i]], mj = mass[this.el[j]];
    const vix = (x[i] - px[i]) / dt, viy = (y[i] - py[i]) / dt;
    const vjx = (x[j] - px[j]) / dt, vjy = (y[j] - py[j]) / dt;
    const M = mi + mj;
    const cx = (mi * vix + mj * vjx) / M, cy = (mi * viy + mj * vjy) / M;
    const a = this.rng() * Math.PI * 2;
    let rix = Math.cos(a) / mi, riy = Math.sin(a) / mi, rjx = -Math.cos(a) / mj, rjy = -Math.sin(a) / mj;
    const k0 = 0.5 * mi * (rix * rix + riy * riy) + 0.5 * mj * (rjx * rjx + rjy * rjy);
    const s = Math.sqrt(e / k0);
    rix *= s; riy *= s; rjx *= s; rjy *= s;
    px[i] = x[i] - (cx + rix) * dt; py[i] = y[i] - (cy + riy) * dt;
    px[j] = x[j] - (cx + rjx) * dt; py[j] = y[j] - (cy + rjy) * dt;
  };

  // Non-bonded potential energy between i and j, mirroring the forces in computeForces().
  S.pairPotential = function (i, j) {
    const ei = this.el[i], ej = this.el[j], p = ei * NZ + ej;
    const dx = this.x[j] - this.x[i], dy = this.y[j] - this.y[i];
    const d = Math.sqrt(dx * dx + dy * dy) || 1e-6;
    const rc = pairRc[p];
    const ionic = this.ion[i] && this.ion[j], rcq = ionic ? pairRcIon[p] : rc;
    if (d >= rc && (!ionic || d >= rcq)) return 0;
    if (this.isBonded(i, j) || this.isOneThree(i, j)) return 0;
    if (d >= rc) return P.KC * P.ION_KC * this.q[i] * this.q[j] * (1 / d - 1 / rcq + (d - rcq) / (rcq * rcq));
    const nb = this.nb, q = this.q;
    let sigma = ionic ? (q[i] > 0 ? catR[ei] : anR[ei]) + (q[j] > 0 ? catR[ej] : anR[ej]) : pairSigma[p];
    const hbs = this.hbd[i] ? HB_ACCEPTOR[ej] : this.hbd[j] ? HB_ACCEPTOR[ei] : 0;
    if (hbs) sigma *= P.HB_CONTACT;
    else if (this.hbd[i] && this.hbd[j]) sigma *= P.HH_CONTACT;
    if ((this.ghostT[i] || this.ghostT[j]) && this.isGhost(i, j) && d < sigma) return 0;
    let eps;
    if (isMetal[ei] && isMetal[ej] && !nb[i] && !nb[j]) eps = pairMetalEps[p];
    else {
      const bi = nb[i] > 2 ? 2 : nb[i], bj = nb[j] > 2 ? 2 : nb[j];
      eps = pairEps[p] * P.VDW_SCALE * (1 + P.MOL_POLARIZABILITY * bi) * (1 + P.MOL_POLARIZABILITY * bj);
    }
    const cf = this.coh[i] === this.coh[j] ? this.coh[i] : Math.sqrt(this.coh[i] * this.coh[j]);
    if (!(isMetal[ei] && isMetal[ej] && !nb[i] && !nb[j])) eps *= cf;
    let U;
    if (d < sigma) {
      const mu = mass[ei] * mass[ej] / (mass[ei] + mass[ej]);
      U = -eps + 0.5 * (hbs ? P.HB_REP_SCALE : 1) * P.K_REP * mu * (sigma - d) * (sigma - d);
    } else {
      const s = (d - sigma) / (rc - sigma);
      U = -eps * (1 - s * s) * (1 - s * s);
    }
    if (hbs) U += this.hbondEnergy(i, j, d, sigma, rc, hbs * P.HB_EPS * P.KB * cf);
    const qq = this.q[i] * this.q[j];
    if (qq !== 0) {
      const d0 = sigma * 0.8, kc = ionic ? P.KC * P.ION_KC : P.KC;
      const irc2 = 1 / (rcq * rcq);
      U += d >= d0 ? kc * qq * (1 / d - 1 / rcq + (d - rcq) * irc2)
        : kc * qq * ((1 / d0 - 1 / rcq + (d0 - rcq) * irc2) + (d0 - d) * (1 / (d0 * d0) - irc2));
    }
    return U;
  };

  // Energy of the directional hydrogen bond between donor H and acceptor (see hbond in computeForces).
  S.hbondEnergy = function (i, j, d, sigma, rc, G) {
    const h = this.hbd[i] && HB_ACCEPTOR[this.el[j]] ? i : j, A = h === i ? j : i;
    const D = this.bp[h * P.MAXB], x = this.x, y = this.y;
    const ux = x[h] - x[D], uy = y[h] - y[D], L = Math.sqrt(ux * ux + uy * uy) || 1e-6;
    const c = (ux * (x[A] - x[h]) + uy * (y[A] - y[h])) / (L * d);
    if (c <= 0) return 0;
    let b = 1;
    if (this.nb[A] === 2 && P.HB_ACCEPT_POW) {
      const N1 = this.bp[A * P.MAXB], N2 = this.bp[A * P.MAXB + 1];
      const a1x = x[N1] - x[A], a1y = y[N1] - y[A], l1 = Math.sqrt(a1x * a1x + a1y * a1y) || 1e-6;
      const a2x = x[N2] - x[A], a2y = y[N2] - y[A], l2 = Math.sqrt(a2x * a2x + a2y * a2y) || 1e-6;
      const sx = a1x / l1 + a2x / l2, sy = a1y / l1 + a2y / l2, S = Math.sqrt(sx * sx + sy * sy);
      const t = (S - HB_FLAT_LO) / (HB_FLAT_HI - HB_FLAT_LO);
      if (t > 0) {
        const tc = t >= 1 ? 1 : t, lam = tc * tc * (3 - 2 * tc);
        const q = -(sx * (x[h] - x[A]) + sy * (y[h] - y[A])) / (S * d);
        b = 1 - lam + lam * (q > 0 ? Math.pow(q, P.HB_ACCEPT_POW) : 0);
        if (b <= 0) return 0;
      }
    }
    let w = 1;
    if (d > sigma) { const s = (d - sigma) / (rc - sigma); w = (1 - s * s) * (1 - s * s); }
    return -G * w * b * Math.pow(c, P.HB_ANGLE_POW);
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
    if (P.W3_EPS) U += this.waterAngleEnergyNear(set, count);
    if (P.W2_EPS) U += this.waterDenseEnergyNear(set, count);
    if (P.METAL_EMBED) U += this.metalEnergyNear(set, count);
    return U;
  };

  // Metallic density of free metal atom i (0 for anything else); with `tag`, only neighbours
  // stamped (sign > 0) or not stamped (sign < 0) count
  S.metalRho = function (i, tag, sign) {
    const el = this.el, nb = this.nb, x = this.x, y = this.y, mark = this.mark;
    const ei = el[i];
    if (!isMetal[ei] || nb[i]) return 0;
    let r = 0;
    this.query(x[i], y[i], this.cellSize, function (j) {
      if (j === i || !isMetal[el[j]] || nb[j]) return;
      if (sign && (mark[j] === tag) !== (sign > 0)) return;
      const p = ei * NZ + el[j], rc = pairRhoRc[p];
      const dx = x[j] - x[i], dy = y[j] - y[i], d2 = dx * dx + dy * dy;
      if (d2 >= rc * rc) return;
      const u = 1 - d2 / (rc * rc);
      r += u * u * u;
    });
    return r;
  };
  // a rate given per substep of DT_REF, for a substep of DT
  function perSub(r) { return 1 - Math.pow(1 - r, P.DT / P.DT_REF); }
  // ln of the real collisions one contact-step represents (LN_COLL is for a substep of DT_REF)
  function lnColl() { return P.LN_COLL + Math.log(P.DT / P.DT_REF); }

  // the embedding energy F(rho) of a free metal atom and its derivative
  // (smootherstep in the density: its curvature too is continuous, which keeps the energy)
  function embed(z, r) {
    const t = (r - P.METAL_RHO_A) / (P.METAL_RHO_B - P.METAL_RHO_A);
    return t <= 0 ? 0 : t >= 1 ? -metalA[z] : -metalA[z] * t * t * t * (t * (6 * t - 15) + 10);
  }
  function embedD(z, r) {
    const t = (r - P.METAL_RHO_A) / (P.METAL_RHO_B - P.METAL_RHO_A);
    return t <= 0 || t >= 1 ? 0 : -metalA[z] * 30 * t * t * (1 - t) * (1 - t) / (P.METAL_RHO_B - P.METAL_RHO_A);
  }

  // Embedding energy that changes when the atoms of set move or change: their own, plus that of
  // their metal neighbours outside the set, measured relative to the neighbours' density without
  // the set (which the event cannot change), so differences are exact.
  S.metalEnergyNear = function (set, count) {
    const el = this.el, nb = this.nb, x = this.x, y = this.y, mark = this.mark;
    let any = false;
    for (let a = 0; a < count; a++) if (isMetal[el[set[a]]]) { any = true; break; }
    if (!any) return 0;
    const tag = ++this.markTag;
    for (let a = 0; a < count; a++) mark[set[a]] = tag;
    const near = this._fsNear || (this._fsNear = []);
    near.length = 0;
    const tag2 = ++this.markTag; // second stamp for 'already listed' neighbours
    const listed = this._fsListed || (this._fsListed = new Int32Array(this.mark.length));
    let U = 0;
    for (let a = 0; a < count; a++) {
      const i = set[a];
      if (!isMetal[el[i]]) continue;
      U += embed(el[i], this.metalRho(i));
      this.query(x[i], y[i], this.cellSize, function (j) {
        if (mark[j] !== tag && listed[j] !== tag2 && isMetal[el[j]] && !nb[j]) { listed[j] = tag2; near.push(j); }
      });
    }
    for (const j of near) {
      const rIn = this.metalRho(j, tag, 1);
      if (rIn === 0) continue;
      const rOut = this.metalRho(j, tag, -1);
      U += embed(el[j], rIn + rOut) - embed(el[j], rOut);
    }
    return U;
  };

  // Water's hydrogen-bond network angle, after the three-body term of the coarse-grained mW water
  // model (Molinero & Moore, 2009): neighbouring water oxygens prefer O-O-O angles of 120 degrees,
  // the honeycomb that hexagonal ice shows when seen along its c axis. Dense packing (60 degree
  // angles) and square packing (90) cost energy, which makes ice open and less dense than water.
  // U = K sum over centres i and neighbour pairs j<k of f(r_ij) f(r_ik) (cos theta_jik + 1/2)^2,
  // with f = 1 out to W3_INNER and then fading smoothly to zero at the range W3_RANGE.
  S.isWaterO = function (i) {
    if (this.el[i] !== 8 || this.nb[i] !== 2) return false;
    const b = i * P.MAXB;
    return this.el[this.bp[b]] === 1 && this.el[this.bp[b + 1]] === 1;
  };
  // inSet (optional): a stamp array and value; only angle terms involving a stamped atom count
  S.waterAngle = function (i, forces, inSet, tag) {
    const rc = P.W3_RANGE, K = P.W3_EPS * P.KB, x = this.x, y = this.y, fx = this.fx, fy = this.fy;
    const L = this._w3 || (this._w3 = []);
    L.length = 0;
    const self = this;
    this.query(x[i], y[i], rc, function (j) { if (j !== i && self.isWaterO(j)) L.push(j); });
    let U = 0;
    // weights: 1 out to an inner radius (so hydrogen-bonded neighbours feel no radial push from
    // the angle term), then a smooth fade to zero
    const r1 = P.W3_INNER, iw = 1 / (rc - r1);
    function wt(r) { if (r <= r1) return 1; if (r >= rc) return 0; const g = (r - r1) * iw, h = 1 - g * g; return h * h; }
    function dwt(r) { if (r <= r1 || r >= rc) return 0; const g = (r - r1) * iw; return -4 * g * (1 - g * g) * iw; }
    for (let p = 0; p < L.length; p++) {
      const j = L[p], ax = x[j] - x[i], ay = y[j] - y[i];
      const ra = Math.sqrt(ax * ax + ay * ay) || 1e-6;
      const fa = wt(ra), dfa = dwt(ra);
      if (!fa) continue;
      for (let r = p + 1; r < L.length; r++) {
        const k = L[r];
        if (inSet && inSet[i] !== tag && inSet[j] !== tag && inSet[k] !== tag) continue;
        const bx = x[k] - x[i], by = y[k] - y[i];
        const rb = Math.sqrt(bx * bx + by * by) || 1e-6;
        const fb = wt(rb), dfb = dwt(rb);
        if (!fb) continue;
        const c = (ax * bx + ay * by) / (ra * rb), dc = c + 0.5;
        U += K * fa * fb * dc * dc;
        if (!forces) continue;
        const kc = 2 * K * fa * fb * dc; // dE/dc
        const ja = K * dc * dc * fb * dfa / ra, kb = K * dc * dc * fa * dfb / rb; // radial parts
        const jx = kc * (bx / (ra * rb) - c * ax / (ra * ra)) + ja * ax;
        const jy = kc * (by / (ra * rb) - c * ay / (ra * ra)) + ja * ay;
        const kx = kc * (ax / (ra * rb) - c * bx / (rb * rb)) + kb * bx;
        const ky = kc * (ay / (ra * rb) - c * by / (rb * rb)) + kb * by;
        fx[j] -= jx; fy[j] -= jy; fx[k] -= kx; fy[k] -= ky;
        fx[i] += jx + kx; fy[i] += jy + ky;
      }
    }
    return U;
  };
  // The network-angle energy of all angle terms that involve a set of atoms (for event
  // accounting: terms not touching the set cannot change, so they are left out on both sides).
  S.waterAngleEnergyNear = function (set, count) {
    const stamp = this._w3s || (this._w3s = new Uint32Array(P.MAX));
    const mark = this._w3m || (this._w3m = new Uint32Array(P.MAX));
    this._w3t = (this._w3t || 0) + 1;
    const t = this._w3t, centres = this._w3c || (this._w3c = []), self = this;
    centres.length = 0;
    let any = false;
    for (let a = 0; a < count; a++) { mark[set[a]] = t; if (this.isWaterO(set[a])) any = true; }
    if (!any) return 0;
    for (let a = 0; a < count; a++) {
      const i = set[a];
      this.query(this.x[i], this.y[i], P.W3_RANGE + 8, function (j) {
        if (stamp[j] !== t && self.isWaterO(j)) { stamp[j] = t; centres.push(j); }
      });
    }
    let U = 0;
    for (let k = 0; k < centres.length; k++) U += this.waterAngle(centres[k], false, mark, t);
    return U;
  };

  // Dense-water cohesion (see W2_EPS): each water oxygen's neighbour count rho is a sum of
  // w(r) = 1 out to W2_INNER, fading to 0 at W2_RANGE; its energy is F(rho) = -C smoothstep.
  function w2(r) { if (r <= P.W2_INNER) return 1; if (r >= P.W2_RANGE) return 0; const g = (r - P.W2_INNER) / (P.W2_RANGE - P.W2_INNER), h = 1 - g * g; return h * h; }
  function dw2(r) { if (r <= P.W2_INNER || r >= P.W2_RANGE) return 0; const iq = 1 / (P.W2_RANGE - P.W2_INNER), g = (r - P.W2_INNER) * iq; return -4 * g * (1 - g * g) * iq; }
  function w2F(r) {
    const t = (r - P.W2_RHO_A) / (P.W2_RHO_B - P.W2_RHO_A);
    return t <= 0 ? 0 : t >= 1 ? -P.W2_EPS * P.KB : -P.W2_EPS * P.KB * t * t * (3 - 2 * t);
  }
  function w2dF(r) {
    const t = (r - P.W2_RHO_A) / (P.W2_RHO_B - P.W2_RHO_A);
    return t <= 0 || t >= 1 ? 0 : -P.W2_EPS * P.KB * 6 * t * (1 - t) / (P.W2_RHO_B - P.W2_RHO_A);
  }
  S.waterDenseForces = function () {
    const n = this.n, x = this.x, y = this.y, fx = this.fx, fy = this.fy, self = this;
    const rho = this.wRho || (this.wRho = new Float64Array(P.MAX));
    const pi = this.wPI || (this.wPI = []), pj = this.wPJ || (this.wPJ = []), pd = this.wPD || (this.wPD = []);
    pi.length = pj.length = pd.length = 0;
    const isW = this.wIs || (this.wIs = new Uint8Array(P.MAX));
    for (let i = 0; i < n; i++) { isW[i] = this.isWaterO(i) ? 1 : 0; rho[i] = 0; }
    const R = P.W2_RANGE;
    for (let i = 0; i < n; i++) {
      if (!isW[i]) continue;
      this.query(x[i], y[i], R, function (j) {
        if (j <= i || !isW[j]) return;
        const dx = x[j] - x[i], dy = y[j] - y[i], d2 = dx * dx + dy * dy;
        if (d2 >= R * R) return;
        const d = Math.sqrt(d2) || 1e-6, w = w2(d);
        rho[i] += w; rho[j] += w;
        const dw = dw2(d);
        if (dw) { pi.push(i); pj.push(j); pd.push(dw / d, dx, dy); }
      });
    }
    for (let k = 0; k < pi.length; k++) {
      const i = pi[k], j = pj[k];
      const c = (w2dF(rho[i]) + w2dF(rho[j])) * pd[3 * k]; // dU/dr / r
      const gx = c * pd[3 * k + 1], gy = c * pd[3 * k + 2];
      fx[j] -= gx; fy[j] -= gy; fx[i] += gx; fy[i] += gy;
    }
  };
  // water-oxygen density of i: all neighbours, or only stamped (sign > 0) / unstamped (sign < 0)
  S.waterRho = function (i, mark, tag, sign) {
    if (!this.isWaterO(i)) return 0;
    const x = this.x, y = this.y, self = this, R = P.W2_RANGE;
    let r = 0;
    this.query(x[i], y[i], R, function (j) {
      if (j === i || !self.isWaterO(j)) return;
      if (sign && (mark[j] === tag) !== (sign > 0)) return;
      const dx = x[j] - x[i], dy = y[j] - y[i], d2 = dx * dx + dy * dy;
      if (d2 < R * R) r += w2(Math.sqrt(d2));
    });
    return r;
  };
  // Dense-water energy that a change of the atoms in set can affect (exact differences, as for
  // metals): their own, plus their water neighbours' relative to the density without the set.
  S.waterDenseEnergyNear = function (set, count) {
    let any = false;
    for (let a = 0; a < count; a++) if (this.el[set[a]] === 8) { any = true; break; }
    if (!any) return 0;
    const mark = this.mark, tag = ++this.markTag, tag2 = ++this.markTag, self = this;
    const listed = this._wListed || (this._wListed = new Int32Array(P.MAX));
    const near = this._wNear || (this._wNear = []);
    near.length = 0;
    for (let a = 0; a < count; a++) mark[set[a]] = tag;
    let U = 0;
    for (let a = 0; a < count; a++) {
      const i = set[a];
      if (!this.isWaterO(i)) continue;
      U += w2F(this.waterRho(i));
      this.query(this.x[i], this.y[i], P.W2_RANGE, function (j) {
        if (mark[j] !== tag && listed[j] !== tag2 && self.isWaterO(j)) { listed[j] = tag2; near.push(j); }
      });
    }
    for (const j of near) {
      const rIn = this.waterRho(j, mark, tag, 1);
      if (rIn === 0) continue;
      const rOut = this.waterRho(j, mark, tag, -1);
      U += w2F(rIn + rOut) - w2F(rOut);
    }
    return U;
  };

  // Bending (angle spring) energy around centre c, mirroring angleForces().
  S.angleEnergy = function (c) {
    const m = this.nb[c];
    if (m < 2) return 0;
    const MB = P.MAXB, base = c * MB, x = this.x, y = this.y, el = this.el, K = P.K_ANGLE;
    const self = this;
    function springU(a, b, target, k) {
      const d = Math.hypot(x[b] - x[a], y[b] - y[a]);
      const ma = mass[el[a]], mb = mass[el[b]];
      return 0.5 * (k || K) * (ma * mb / (ma + mb)) * (d - target) * (d - target);
    }
    if (m === 2) {
      const theta = angleFor(el[c], this.used[c]);
      const la = this.bl[base], lb = this.bl[base + 1];
      return springU(this.bp[base], this.bp[base + 1], Math.sqrt(la * la + lb * lb - 2 * la * lb * Math.cos(theta)), P.K_ANGLE_2);
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
      // candidates were collected before this substep's events; never close a 3-membered ring
      if (this.nb[i] && this.nb[j] && this.isOneThree(i, j)) continue;
      const ei = el[i], ej = el[j], p = ei * NZ + ej;
      const cls = C.pairClass[p];
      if (!cls) continue;
      // an ion does not bond to its own kind (no O-O peroxide chains on an oxide, no Cl-Cl on a salt)
      if (ei === ej && (this.ion[i] || this.ion[j])) continue;
      if (ei === 8 && ej === 8 && !this.oxygenChainOK(i, j)) continue;
      // freshly added solid or liquid still settling in (graphite assembling from added atoms, held
      // at room temperature by a strong thermostat that its rigid bonds make it read hot) does not
      // react with other substances yet (it still bonds within itself)
      if (ei !== ej && (this.settling(i) || this.settling(j))) continue;
      const T = 0.5 * (this.tLoc[i] + this.tLoc[j]);
      let RT = P.R_KJ * Math.max(T, 1);
      const E1 = C.pairEnergy[p * 4 + 1];
      const canI = this.canAccept(i, ej) && !this.loneC(i), canJ = this.canAccept(j, ei) && !this.loneC(j);
      if (canI && canJ) {
        // radical + radical: nearly barrierless recombination
        if (rng() < Math.exp(Math.min(0, lnColl() - 4 / RT))) {
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
        // on the surface of a lump of metal the rate goes with the surface's temperature (the
        // average one already holds the fast tail of the gas; a fast molecule's own reading on top
        // would count that tail twice)
        if (isMetal[el[a]] && this.condF[a] > 0.5) RT = P.R_KJ * Math.max(this.tLoc[a], 1);
        if (isMetal[el[a]]) {
          if (this.dissociativeAdsorption(a, b, RT)) continue;
          // a metal surface does not half-open O2 or Cl2 into a peroxide or superoxide (M-O-O.) or
          // knock a lone O atom off it: diatomics split across two metal atoms (above). A lone metal
          // atom in the gas does pluck one atom off a halogen (Na + Cl2 -> NaCl + Cl, the harpoon
          // reaction), but not off spin-paired O2
          if (this.nb[b] === 1 && el[this.bp[b * P.MAXB]] === el[b] && !isMetal[el[b]] && (el[b] === 8 || this.hasMetalNeighbour(a))) continue;
        }
        const k = this.weakestBond(b, a);
        if (k < 0 || el[k] === el[a]) continue; // X. + X-Y -> X-Y + X. changes nothing
        const step = this._weakE;
        // after losing a unit, b must be able to bond with a (hypervalence rule)
        if (this.used[b] > normVal[el[b]] && !(en[el[a]] > en[el[b]] + 0.3)) continue;
        const dH = step - E1;
        // (a free valence on a solid, such as an edge of graphite or a sulfur chain, is a surface
        // like a metal's: it meets the spin barrier too, so graphite needs heat to burn)
        // An edge of a covalent network also has to chemisorb the oxygen: graphite oxidises with
        // about 170 kJ/mol of activation, so it needs red heat to burn
        const spin = this._weakSpin && (isMetal[el[a]] || this.netF[a] || this.isSolidAtom(a))
          ? P.SPIN_BARRIER + (this.netF[a] && !isMetal[el[a]] ? P.NETWORK_O2_BARRIER : 0) : 0;
        // alkali and alkaline-earth metals react by handing over an electron, with almost no
        // barrier beyond the reaction enthalpy (sodium in water, potassium in air)
        const alpha = isMetal[el[a]] && en[el[a]] < P.ET_EN ? P.ALPHA_ET : P.ALPHA;
        // (and other metals meet the same halide-layer barrier as in dissociativeAdsorption)
        const passive = isMetal[el[a]] && en[el[a]] >= P.ET_EN && HALOGENS.has(el[b]) ? P.PASSIVE_BARRIER : 0;
        const Ea = alpha * step + Math.max(0, dH) + spin + passive;
        if (rng() < Math.exp(Math.min(0, lnColl() - Ea / RT)) && this.canAfford(a, b, dH)) {
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
        if (ord >= 2 && !(el[a] === 8 && el[b] === 8 && this.isBonded(b, k) && el[k] === 8 && !this.oxygenChainOK(a, b))) {
          const full = C.bondEnergy(el[b], el[k], ord);
          const dHt = full - E1;
          const Eat = P.ALPHA * full + Math.max(0, dHt) + spin;
          if (Eat / RT < lnColl() + 12 && rng() < Math.exp(Math.min(0, lnColl() - Eat / RT)) && this.canAfford(a, b, dHt)) {
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
      if (k !== l && !this.isBonded(k, l) && !this.isOneThree(k, l)) {
        const pkl = el[k] * NZ + el[l];
        const dx = this.x[l] - this.x[k], dy = this.y[l] - this.y[k];
        const near = dx * dx + dy * dy < 2.25 * pairSigma[pkl] * pairSigma[pkl];
        const ooOK = !(el[k] === 8 && el[l] === 8) || this.oxygenChainOK(k, l, i, j);
        if (C.pairClass[pkl] && near && ooOK) { formSum += C.pairEnergy[pkl * 4 + 1]; formKL = true; }
      }
      const dH = breakSum - formSum;
      // AB + CD -> AC + B. + D. makes one bond and two radicals: an abstraction by a molecule, with
      // the abstraction barrier on the weaker bond (H2 + F2 -> HF + H + F starts the dark explosion
      // of hydrogen and fluorine; with chlorine or oxygen the same step is far uphill). Swapping
      // both pairs at once (formKL) is the symmetry-forbidden four-centre route
      const Ea = (formKL ? P.ALPHA_4C * breakSum : P.ALPHA * Math.min(sI, sJ)) + Math.max(0, dH) + spin4;
      if (Ea / RT > lnColl() + 12) continue;
      if (rng() < Math.exp(Math.min(0, lnColl() - Ea / RT)) && this.canAfford(i, j, dH)) {
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

  // The carbon of carbon monoxide (C triple-bonded to O, nothing else): its spare valence is a
  // lone pair, not an unpaired electron, so CO is a closed-shell molecule. It does not stick to
  // other CO molecules; radicals still add to it (O + CO -> CO2, OH + CO -> HOCO -> CO2 + H).
  S.loneC = function (i) {
    if (this.el[i] !== 6 || this.nb[i] !== 1) return false;
    return this.el[this.bp[i * P.MAXB]] === 8 && this.bo[i * P.MAXB] === 3;
  };

  // May oxygens i and j bond? Not into an O-O-O chain on a molecule that holds other elements:
  // H2O3 and longer polyoxides (and organic HOOO- groups) fall apart far below room temperature,
  // so OH + H2O2 does not grow H2O3. Pure oxygen may, as far as ozone: O + O2 -> O3, but not
  // O2 + O2 -> O4.
  // (si, sj: partners i and j are losing in the same step, not counted)
  S.oxygenChainOK = function (i, j, si, sj) {
    const MB = P.MAXB, el = this.el;
    let hetero = false;
    const oN = [0, 0];
    for (let t = 0; t < 2; t++) {
      const a = t ? j : i;
      for (let k = 0; k < this.nb[a]; k++) {
        const b = this.bp[a * MB + k];
        if (b === (t ? sj : si)) continue;
        if (el[b] === 8) oN[t]++; else hetero = true;
      }
    }
    if (oN[0] + oN[1] && hetero) return false;
    return !(oN[0] && oN[1]) && oN[0] < 2 && oN[1] < 2;
  };

  S.settling = function (i) {
    const z = this.el[i], bp = E[z].bp || 1e4;
    // (a non-volatile element's atoms count even while knocked loose into the gas: carbon atoms
    // flung off graphite as it assembles)
    if (this.anneal[i] > 0 && (this.condF[i] > 0.5 || this.nb[i] >= 2 || this.netF[i] === 1 || bp > P.AMBIENT)) return true;
    // a lone atom of a non-volatile non-metal far below its boiling point (carbon has no vapour
    // until about 3900 K) is a stray from its lump: it may only rejoin it
    return this.nb[i] === 0 && !isMetal[z] && bp > 600 && this.condF[i] < 0.5 && this.tCell[this.cellOfT[i]] < 0.5 * bp;
  };

  S.hasMetalNeighbour = function (a) {
    const self = this, x = this.x, y = this.y, el = this.el, r = pairSigma[el[a] * NZ + el[a]] * 1.5;
    let found = false;
    this.query(x[a], y[a], r, function (j) {
      if (found || j === a || !isMetal[el[j]]) return;
      if ((x[j] - x[a]) ** 2 + (y[j] - y[a]) ** 2 < r * r) found = true;
    });
    return found;
  };

  // A diatomic molecule X2 (O2, Cl2, ...) hitting a metal surface splits across two metal atoms:
  // a takes b, a neighbouring metal atom m takes b's partner k (2 M + O2 -> 2 M-O). This is how
  // metals oxidise and halogenate, and it leaves no stray O atoms or peroxides behind. Returns
  // true when it happened.
  S.dissociativeAdsorption = function (a, b, RT) {
    if (this.nb[b] !== 1) return false;
    const MB = P.MAXB, k = this.bp[b * MB], el = this.el;
    if (el[k] !== el[b] || this.nb[k] !== 1 || isMetal[el[b]]) return false;
    const eX = el[b], ord = this.bondOrder(b, k);
    const pAB = el[a] * NZ + eX;
    if (!C.pairClass[pAB]) return false;
    // a second metal atom with a free valence: a surface neighbour of a, the one nearest to k (the
    // molecule tips over onto it as it splits)
    const self = this, x = this.x, y = this.y;
    let m = -1, best = Infinity;
    this.query(x[a], y[a], pairSigma[el[a] * NZ + el[a]] * 1.5, function (j) {
      if (j === a || j === b || j === k || !isMetal[el[j]] || !self.canAccept(j, eX)) return;
      const ra = pairSigma[el[a] * NZ + el[j]] * 1.5;
      if ((x[j] - x[a]) ** 2 + (y[j] - y[a]) ** 2 > ra * ra) return;
      const r = pairSigma[el[j] * NZ + eX] * 2.2, d2 = (x[j] - x[k]) ** 2 + (y[j] - y[k]) ** 2;
      if (d2 < r * r && d2 < best) { best = d2; m = j; }
    });
    if (m < 0) return false;
    const pMK = el[m] * NZ + eX;
    if (!C.pairClass[pMK]) return false;
    const eXX = C.bondEnergy(eX, eX, ord);
    const dH = eXX - C.pairEnergy[pAB * 4 + 1] - C.pairEnergy[pMK * 4 + 1];
    const spin = eX === 8 && ord === 2 ? P.SPIN_BARRIER : 0;
    // only the electropositive metals hand an electron straight to the molecule; on iron, copper or
    // aluminium the first layer of halide has to be broken through (iron burns in chlorine only once
    // heated to about 500 K; chlorine is stored in steel cylinders)
    const passive = en[el[a]] >= P.ET_EN && !spin ? P.PASSIVE_BARRIER : 0;
    const Ea = P.ALPHA_DISS * eXX + Math.max(0, dH) + spin + passive;
    if (!(this.rng() < Math.exp(Math.min(0, lnColl() - Ea / RT))) || !this.canAfford(a, b, dH)) return false;
    this.beginEvent(a, b, k, m);
    this.setBond(b, k, 0);
    this.addBond(a, b, 1, C.pairClass[pAB]);
    this.addBond(m, k, 1, C.pairClass[pMK]);
    this.endEvent(a, b);
    this.heat(a, b, -dH);
    this.logReaction();
    return true;
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
          if (rng() < Math.exp(Math.min(0, lnColl() - 8 / RT))) {
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

  // Looks up each small molecule's melting point (from the phase data), so frozen molecular
  // solids are recognised as solids. Molecules change only by reactions, so this runs rarely.
  S.updateMeltPoints = function () {
    const n = this.n, MB = P.MAXB, el = this.el, meltT = this.meltT;
    if (!CE.phase) return;
    const seen = this._mpSeen && this._mpSeen.length >= n ? this._mpSeen : (this._mpSeen = new Uint8Array(P.MAX));
    seen.fill(0, 0, n);
    const cache = this._mpCache || (this._mpCache = new Map());
    const mol = [];
    for (let i = 0; i < n; i++) {
      if (seen[i]) continue;
      mol.length = 0; mol.push(i); seen[i] = 1;
      for (let h = 0; h < mol.length; h++) {
        const a = mol[h];
        for (let k = 0; k < this.nb[a]; k++) {
          const b = this.bp[a * MB + k];
          if (!seen[b]) { seen[b] = 1; mol.push(b); }
        }
      }
      // a covalent network (graphite, silicon, sulfur chains): big and held by bonds, not ions
      const net = mol.length > 12 && !this.hasIonicBond(mol) ? 1 : 0;
      for (let k = 0; k < mol.length; k++) this.netF[mol[k]] = net;
      let T = 0, c = 1;
      if (mol.length > 1 && mol.length <= 12) {
        const counts = {};
        for (let k = 0; k < mol.length; k++) counts[el[mol[k]]] = (counts[el[mol[k]]] || 0) + 1;
        const sig = JSON.stringify(counts);
        let e = cache.get(sig);
        if (e === undefined) {
          const key = C.describe(counts, mol.length, false).key;
          const ov = Simulation.cohesionOverride;
          e = { T: CE.phase.meltingPoint(key), c: (ov && ov[key]) || COHESION[key] || 1 };
          cache.set(sig, e);
        }
        T = e.T; c = e.c;
      }
      for (let k = 0; k < mol.length; k++) meltT[mol[k]] = T;
      if (this.coh[i] !== c) {
        // the molecule became another substance (a reaction): its cohesion changes, and the
        // energy that releases or costs is settled as heat
        const set = this._cohSet || (this._cohSet = new Int32Array(12));
        for (let k = 0; k < mol.length; k++) set[k] = mol[k];
        const U0 = this.eventEnergy(set, mol.length);
        for (let k = 0; k < mol.length; k++) this.coh[mol[k]] = c;
        const dU = this.eventEnergy(set, mol.length) - U0;
        if (dU !== 0) this.heatRaw(mol[0], mol.length > 1 ? mol[1] : mol[0], -dU);
      }
    }
  };


  S.hasIonicBond = function (mol) {
    const MB = P.MAXB;
    for (const a of mol) for (let b = 0; b < this.nb[a]; b++) if (this.bc[a * MB + b] === 2) return true;
    return false;
  };

  // Solid-state annealing. Like reactions, molecular rearrangement in a solid is time-compressed:
  // every simulated moment stands for a great many real vibrations, during which molecules in a
  // real crystal turn into their best orientation. Small molecules that are packed in a liquid
  // or solid try a small rigid rotation, accepted with the Boltzmann (Metropolis) rule at the
  // local temperature. That only speeds up reaching equilibrium (a cooled liquid freezes into a
  // crystal instead of a glass); the energy change is paid from or given to the local heat.
  S.solidAnneal = function () {
    const n = this.n, frac = P.ANNEAL_TRIALS * P.DT / P.DT_REF;
    if (!frac || !n) return;
    const x = this.x, y = this.y, el = this.el, rng = this.rng;
    const set = this._anSet || (this._anSet = new Int32Array(6));
    const ox = this._anX || (this._anX = new Float64Array(12));
    const trials = Math.max(1, Math.round(n * frac));
    // fresh neighbour lists: a missed neighbour would make some moves look free, and the
    // Metropolis rule would pick exactly those, creating energy
    this.buildGrid();
    for (let t = 0; t < trials; t++) {
      const r = Math.floor(rng() * n);
      // only frozen molecular solids: packed, and colder than their documented melting point
      if (!this.frozenAt(r)) continue;
      const m = this.gatherMolecule(r, set, 0, 6);
      if (m < 0) continue;
      // a small turn and shift of the whole molecule about its centre of mass
      let M = 0, cx = 0, cy = 0;
      for (let a = 0; a < m; a++) { const mi = mass[el[set[a]]]; cx += mi * x[set[a]]; cy += mi * y[set[a]]; M += mi; }
      cx /= M; cy /= M;
      const th = (rng() * 2 - 1) * P.ANNEAL_ANGLE;
      const sx = (rng() * 2 - 1) * P.ANNEAL_SHIFT, sy = (rng() * 2 - 1) * P.ANNEAL_SHIFT;
      this.annealMove(set, m, cx, cy, th, sx, sy, ox, r);
    }
  };

  // Is atom i part of a frozen molecular solid (packed, away from the glass, below its melting point)?
  S.frozenAt = function (i) {
    if (!this.nb[i] || this.grainOf[i] >= 0 || this.condF[i] < 0.5 || !(this.tLoc[i] < this.meltT[i])) return false;
    const m = 40;
    return this.x[i] > m && this.x[i] < P.W - m && this.y[i] > m && this.y[i] < P.H - m;
  };

  // Collects the molecule of atom r into set[from..] (at most `cap` atoms in all); returns the
  // new end index, or -1 when the molecule is too big or not all of it is frozen.
  S.gatherMolecule = function (r, set, from, cap) {
    const MB = P.MAXB;
    let m = from;
    set[m++] = r;
    for (let h = from; h < m; h++) {
      const a = set[h];
      if (!this.frozenAt(a) && this.nb[a] > 0 && this.el[a] !== 1) return -1;
      for (let k = 0; k < this.nb[a]; k++) {
        const b = this.bp[a * MB + k];
        let seen = false;
        for (let c = from; c < m; c++) if (set[c] === b) { seen = true; break; }
        if (seen) continue;
        if (m >= cap || m - from >= 6) return -1;
        set[m++] = b;
      }
    }
    return m;
  };

  // Tries a rigid turn (th about cx, cy) plus shift of the atoms in set; Metropolis at the local
  // temperature, the energy change settled with the local heat. Returns true if accepted.
  S.annealMove = function (set, m, cx, cy, th, sx, sy, ox, r) {
    const x = this.x, y = this.y, px = this.px, py = this.py, el = this.el;
    const U0 = this.eventEnergy(set, m);
    const c = Math.cos(th), sn = Math.sin(th);
    let dG = 0;
    for (let a = 0; a < m; a++) {
      const i = set[a];
      ox[2 * a] = x[i]; ox[2 * a + 1] = y[i];
      const rx = x[i] - cx, ry = y[i] - cy;
      x[i] = cx + sx + rx * c - ry * sn; y[i] = cy + sy + rx * sn + ry * c;
      dG -= mass[el[i]] * this.gAcc[i] * (y[i] - ox[2 * a + 1]); // gravity's potential energy
    }
    const dU = this.eventEnergy(set, m) - U0 + dG;
    const RT = P.KB * Math.max(this.tLoc[r], 1);
    if (dU > 0 && this.rng() >= Math.exp(-dU / RT)) {
      for (let a = 0; a < m; a++) { const i = set[a]; x[i] = ox[2 * a]; y[i] = ox[2 * a + 1]; }
      return false;
    }
    // accepted: carry each atom's velocity round with it (about the molecule's own drift: turning
    // the drift too would turn the motion of a moving lump into heat), and settle the energy change
    // locally
    let M = 0, ux = 0, uy = 0;
    for (let a = 0; a < m; a++) {
      const i = set[a], mi = mass[el[i]];
      M += mi; ux += mi * (ox[2 * a] - px[i]); uy += mi * (ox[2 * a + 1] - py[i]);
    }
    ux /= M; uy /= M;
    for (let a = 0; a < m; a++) {
      const i = set[a], vx = ox[2 * a] - px[i] - ux, vy = ox[2 * a + 1] - py[i] - uy;
      px[i] = x[i] - (ux + vx * c - vy * sn); py[i] = y[i] - (uy + vx * sn + vy * c);
    }
    if (dU !== 0) this.heatRaw(set[0], set[1], -dU);
    return true;
  };


  // ---------------------------------------------------------------------------------------
  // Crystal engine
  // ---------------------------------------------------------------------------------------
  // A real crystal settles into its lattice over countless molecular vibrations, far more than
  // the time-compressed simulation can show one by one; left alone, water cooled on screen
  // freezes into a disordered network (a glass). So frozen material is gathered into crystal
  // grains. Each grain holds the ideal lattice of its substance: the honeycomb of ice (hexagonal
  // ice seen along its c axis, every edge carrying one hydrogen, alternate oxygens with a
  // dangling H) or the close-packed layer of a metal. Members are pulled toward the best-fitting
  // placement of that template (shape matching: the template may translate and rotate, so a
  // block of ice still falls and tips). Because the fit is re-optimised every step the pull is a
  // true potential energy, U = k/2 sum m |x - goal|^2, and every join and leave is settled with
  // the local heat, so energy stays conserved. Material colder than its documented melting point
  // joins a neighbouring grain at a free lattice site (or starts one, most readily on the glass,
  // aligned with it); a member warmer than the melting point leaves (melts).

  // 1: a water oxygen (ice), 2: a free metal atom, 0: nothing that crystallises here
  S.crystalKind = function (i) {
    if (this.isWaterO(i)) return 1;
    if (isMetal[this.el[i]] && this.nb[i] === 0) return 2;
    return 0;
  };
  S.crystalMeltT = function (i, kind) { return kind === 1 ? (this.meltT[i] || 273.15) : meltK[this.el[i]]; };
  S.crystalSpacing = function (kind, z) { return kind === 1 ? P.ICE_A : pairSigma[z * NZ + z]; };

  // lattice site nearest to rest-frame point (rx, ry): { key, x, y, sub }
  S.crystalSite = function (g, rx0, ry0) {
    const a = g.a, phi = g.phi, ox = g.ox, oy = g.oy, rx = rx0 - ox, ry = ry0 - oy;
    let best = null, bd = Infinity;
    if (g.kind === 2) {
      // triangular: t1 = a (cos phi, sin phi), t2 = a (cos(phi+60), sin(phi+60))
      const t1x = a * Math.cos(phi), t1y = a * Math.sin(phi), t2x = a * Math.cos(phi + Math.PI / 3), t2y = a * Math.sin(phi + Math.PI / 3);
      const det = t1x * t2y - t1y * t2x, f1 = (rx * t2y - ry * t2x) / det, f2 = (t1x * ry - t1y * rx) / det;
      for (let d1 = 0; d1 <= 1; d1++) for (let d2 = 0; d2 <= 1; d2++) {
        const n1 = Math.floor(f1) + d1, n2 = Math.floor(f2) + d2;
        const sx = n1 * t1x + n2 * t2x, sy = n1 * t1y + n2 * t2y, dd = (sx - rx) * (sx - rx) + (sy - ry) * (sy - ry);
        if (dd < bd) { bd = dd; best = { key: n1 + ',' + n2, x: sx + ox, y: sy + oy, sub: 0 }; }
      }
      return best;
    }
    // honeycomb: A sites at n1 b1 + n2 b2 with b1 = a (e0 - e1), b2 = a (e0 - e2); B = A + a e0
    const e = g.e, b1x = a * (e[0] - e[2]), b1y = a * (e[1] - e[3]), b2x = a * (e[0] - e[4]), b2y = a * (e[1] - e[5]);
    const det = b1x * b2y - b1y * b2x;
    for (let sub = 0; sub <= 1; sub++) {
      const qx = rx - sub * a * e[0], qy = ry - sub * a * e[1];
      const f1 = (qx * b2y - qy * b2x) / det, f2 = (b1x * qy - b1y * qx) / det;
      for (let d1 = 0; d1 <= 1; d1++) for (let d2 = 0; d2 <= 1; d2++) {
        const n1 = Math.floor(f1) + d1, n2 = Math.floor(f2) + d2;
        const sx = n1 * b1x + n2 * b2x + sub * a * e[0], sy = n1 * b1y + n2 * b2y + sub * a * e[1];
        const dd = (sx - rx) * (sx - rx) + (sy - ry) * (sy - ry);
        if (dd < bd) { bd = dd; best = { key: n1 + ',' + n2 + ',' + sub, x: sx + ox, y: sy + oy, sub: sub }; }
      }
    }
    return best;
  };

  // Best rigid placement of grain g (centre, template centre, angle) and its spring energy.
  S.crystalFit = function (g) {
    const x = this.x, y = this.y;
    let M = 0, cx = 0, cy = 0, sx = 0, sy = 0;
    for (let k = 0; k < g.mem.length; k++) {
      const u = g.mem[k];
      for (let q = 0; q < u.atoms.length; q++) {
        const i = u.atoms[q], m = mass[this.el[i]];
        M += m; cx += m * x[i]; cy += m * y[i]; sx += m * u.rest[2 * q]; sy += m * u.rest[2 * q + 1];
      }
    }
    cx /= M; cy /= M; sx /= M; sy /= M;
    let num = 0, den = 0;
    for (let k = 0; k < g.mem.length; k++) {
      const u = g.mem[k];
      for (let q = 0; q < u.atoms.length; q++) {
        const i = u.atoms[q], m = mass[this.el[i]];
        const ax = u.rest[2 * q] - sx, ay = u.rest[2 * q + 1] - sy, bx = x[i] - cx, by = y[i] - cy;
        num += m * (ax * by - ay * bx); den += m * (ax * bx + ay * by);
      }
    }
    g.cx = cx; g.cy = cy; g.sx = sx; g.sy = sy; g.th = Math.atan2(num, den);
    return g;
  };
  S.crystalEnergy = function (g, forces) {
    if (!g.mem.length) return 0;
    this.crystalFit(g);
    const x = this.x, y = this.y, c = Math.cos(g.th), sn = Math.sin(g.th);
    const k = g.kind === 2 ? P.CRYSTAL_K_METAL : P.CRYSTAL_K;
    let U = 0;
    for (let q0 = 0; q0 < g.mem.length; q0++) {
      const u = g.mem[q0];
      for (let q = 0; q < u.atoms.length; q++) {
        const i = u.atoms[q], m = mass[this.el[i]];
        const ax = u.rest[2 * q] - g.sx, ay = u.rest[2 * q + 1] - g.sy;
        const gx = g.cx + ax * c - ay * sn, gy = g.cy + ax * sn + ay * c;
        const dx = x[i] - gx, dy = y[i] - gy, d2 = dx * dx + dy * dy;
        U += 0.5 * k * m * d2;
        if (forces) { this.fx[i] -= k * m * dx; this.fy[i] -= k * m * dy; }
      }
    }
    return U;
  };
  S.crystalForces = function () {
    for (let k = 0; k < this.grains.length; k++) this.crystalEnergy(this.grains[k], true);
  };

  // world position -> grain template frame
  S.crystalToRest = function (g, wx, wy) {
    const c = Math.cos(g.th), sn = Math.sin(g.th), dx = wx - g.cx, dy = wy - g.cy;
    return [g.sx + dx * c + dy * sn, g.sy - dx * sn + dy * c];
  };

  // Adds atom i (a water oxygen or a metal atom) to grain g at site, if the strain of the new
  // spring is at most maxE; the spring then draws it onto the site over the next moments, pushing
  // neighbours aside like a real freezing front. The jump in spring energy is settled with the
  // local heat (the latent heat itself comes from the bonds the ordered layer makes).
  S.crystalJoin = function (g, i, site, maxE) {
    const mem = { atoms: [i], rest: [site.x, site.y], key: site.key };
    const U0 = this.crystalEnergy(g, false);
    g.mem.push(mem); g.sites.set(site.key, mem);
    const U1 = this.crystalEnergy(g, false);
    if (U1 - U0 > maxE) { g.mem.pop(); g.sites.delete(site.key); this.crystalFit(g); return false; }
    this.grainOf[i] = this.grains.indexOf(g);
    if (U1 !== U0) this.heatRaw(i, this.crystalPartner(i), -(U1 - U0));
    return true;
  };
  // how close (in lattice spacings) an atom must come to a site to join a grain, and how far a
  // member may stray from it before it leaves: springless grains (metals) only record order, so
  // an atom joins only when it really sits on a site (a melt next to a crystal is not crystalline)
  S.crystalSlack = function (g, join) {
    if (g.kind === 2 && !P.CRYSTAL_K_METAL) return join ? 0.25 : 0.5;
    return join ? 0.45 : 0.6;
  };
  S.crystalPartner = function (i) { return this.nb[i] ? this.bp[i * P.MAXB] : this.crystalNeighbour(i); };
  S.crystalNeighbour = function (i) {
    let best = i, bd = Infinity;
    const self = this;
    this.query(this.x[i], this.y[i], 60, function (j, d) { if (j !== i && d < bd) { bd = d; best = j; } });
    return best;
  };

  // Removes the unit containing atom i from its grain (melting, reaction, removal).
  S.crystalEject = function (i) {
    const gid = this.grainOf[i];
    if (gid < 0) return;
    const g = this.grains[gid];
    let k = -1;
    for (let q = 0; q < g.mem.length; q++) if (g.mem[q].atoms.indexOf(i) >= 0) { k = q; break; }
    if (k < 0) { this.grainOf[i] = -1; return; }
    const u = g.mem[k];
    this._crystalBusy = true;
    const U0 = this.crystalEnergy(g, false);
    g.mem.splice(k, 1); g.sites.delete(u.key);
    for (let q = 0; q < u.atoms.length; q++) this.grainOf[u.atoms[q]] = -1;
    const U1 = this.crystalEnergy(g, false);
    if (U1 !== U0) this.heatRaw(u.atoms[0], this.crystalPartner(u.atoms[0]), U0 - U1);
    this._crystalBusy = false;
  };

  // particle moved from index `from` to `to` (removeParticle's swap)
  S.crystalRenumber = function (from, to) {
    const g = this.grains[this.grainOf[from]];
    if (!g) return;
    for (let q = 0; q < g.mem.length; q++) {
      const k = g.mem[q].atoms.indexOf(from);
      if (k >= 0) { g.mem[q].atoms[k] = to; return; }
    }
  };

  // Joining, nucleation and melting, every few substeps.
  // Grains of a metal that touch are merged into one grain (bookkeeping: what holds the lump
  // together is the metal's own cohesion). The smaller grain's lattice sites are carried over
  // into the larger grain's frame (keeping their own orientation), and any change in spring
  // energy is settled as heat.
  S.crystalMerge = function () {
    const grains = this.grains, x = this.x, y = this.y, self = this;
    if (grains.length < 2) return;
    const touch = new Map();
    for (let gi = 0; gi < grains.length; gi++) {
      const g = grains[gi];
      if (g.kind !== 2 || !g.mem.length) continue;
      const r = g.a * 1.25;
      for (const u of g.mem) {
        const i = u.atoms[0];
        this.query(x[i], y[i], r, function (j) {
          const gj = self.grainOf[j];
          if (gj < 0 || gj === gi || grains[gj].kind !== 2 || grains[gj].z !== g.z) return;
          if ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 > r * r) return;
          const key = gi < gj ? gi + ':' + gj : gj + ':' + gi;
          touch.set(key, (touch.get(key) || 0) + 1);
        });
      }
    }
    const merged = new Set();
    for (const [key, cnt] of touch) {
      if (cnt < 6) continue; // a real boundary, not a stray touch
      let [ga, gb] = key.split(':').map(Number);
      if (merged.has(ga) || merged.has(gb)) continue;
      let A = grains[ga], B = grains[gb];
      if (B.mem.length > A.mem.length) { const t = A; A = B; B = t; const ti = ga; ga = gb; gb = ti; }
      this._crystalBusy = true;
      const U0 = this.crystalEnergy(A, false) + this.crystalEnergy(B, false);
      const c = Math.cos(B.th), sn = Math.sin(B.th);
      this._mergeId = (this._mergeId || 0) + 1;
      for (const u of B.mem) {
        const rest = [];
        for (let q = 0; q < u.atoms.length; q++) {
          const ax = u.rest[2 * q] - B.sx, ay = u.rest[2 * q + 1] - B.sy;
          const r = this.crystalToRest(A, B.cx + ax * c - ay * sn, B.cy + ax * sn + ay * c);
          rest.push(r[0], r[1]);
        }
        const mu = { atoms: u.atoms, rest: rest, key: 'm' + this._mergeId + ':' + u.key };
        A.mem.push(mu); A.sites.set(mu.key, mu);
        for (const i of u.atoms) this.grainOf[i] = ga;
      }
      B.mem = []; B.sites = new Map();
      const U1 = this.crystalEnergy(A, false);
      const i0 = A.mem[0].atoms[0];
      if (U1 !== U0) this.heatRaw(i0, this.crystalPartner(i0), U0 - U1);
      this._crystalBusy = false;
      merged.add(ga); merged.add(gb);
    }
  };

  S.crystalUpdate = function () {
    const n = this.n, x = this.x, y = this.y, rng = this.rng, self = this;
    if (!P.CRYSTAL_K || !n) return;
    if ((this._crystalCalls = (this._crystalCalls || 0) + 1) % 16 === 0) this.crystalMerge();
    // melting and dislodged members leave
    for (let gi = 0; gi < this.grains.length; gi++) {
      const g = this.grains[gi];
      if (!g.mem.length) continue;
      this.crystalFit(g);
      const c = Math.cos(g.th), sn = Math.sin(g.th);
      for (let q = g.mem.length - 1; q >= 0; q--) {
        const u = g.mem[q], i = u.atoms[0];
        const ax = u.rest[0] - g.sx, ay = u.rest[1] - g.sy;
        const off = Math.hypot(x[i] - (g.cx + ax * c - ay * sn), y[i] - (g.cy + ax * sn + ay * c));
        // a member melts when it is hot itself, not merely when the gas around it is (a lump
        // swung through vapour, or a cold crystal under a flame's hot gas, stays solid)
        if (Math.min(this.tLoc[i], this.tSelf[i]) > this.crystalMeltT(i, g.kind) + P.CRYSTAL_HYST || off > this.crystalSlack(g) * g.a) {
          this.crystalEject(i);
          this.crystalFit(g);
        }
      }
    }
    // tidy away empty grains (renumbering their members' grain ids)
    if (this.grains.some(function (g) { return !g.mem.length; })) {
      const keep = this.grains.filter(function (g) { return g.mem.length; });
      this.grains = keep;
      this.grainOf.fill(-1, 0, n);
      for (let gi = 0; gi < keep.length; gi++) for (let q = 0; q < keep[gi].mem.length; q++) for (const i of keep[gi].mem[q].atoms) this.grainOf[i] = gi;
    }
    // candidates: cold crystallisable units not in a grain
    const tries = Math.max(1, Math.round(n * P.CRYSTAL_TRIES * P.DT / P.DT_REF));
    this.buildGrid();
    for (let t = 0; t < tries; t++) {
      const i = Math.floor(rng() * n);
      if (this.grainOf[i] >= 0) continue;
      const kind = this.crystalKind(i);
      if (!kind || this.condF[i] < 0.5 || !(this.tLoc[i] < this.crystalMeltT(i, kind))) continue;
      const a = this.crystalSpacing(kind, this.el[i]);
      // nearest member of a grain of the same substance
      let best = -1, bd = 1.7 * a;
      this.query(x[i], y[i], 1.7 * a, function (j, d) {
        if (j === i || self.grainOf[j] < 0 || d >= bd) return;
        const g = self.grains[self.grainOf[j]];
        if (g.kind !== kind || (kind === 2 && g.z !== self.el[i]) || g.mem[0] === undefined) return;
        if (kind === 1 && !self.isWaterO(j)) return;
        best = j; bd = d;
      });
      if (best >= 0) {
        const g = this.grains[this.grainOf[best]];
        if (g.mem.length >= P.CRYSTAL_MAX_MEMBERS) continue;
        this.crystalFit(g);
        const r = this.crystalToRest(g, x[i], y[i]);
        const site = this.crystalSite(g, r[0], r[1]);
        if (g.sites.has(site.key) || Math.hypot(site.x - r[0], site.y - r[1]) > this.crystalSlack(g, true) * a) continue;
        this.crystalJoin(g, i, site, (g.kind === 2 ? P.CRYSTAL_JOIN_MAX_METAL : P.CRYSTAL_JOIN_MAX) * P.KB);
        continue;
      }
      // nobody to join: nucleate a new grain (readily on the glass, rarely in the bulk)
      if (this.grains.length >= P.CRYSTAL_MAX_GRAINS) continue;
      let crowd = false;
      this.query(x[i], y[i], 2.5 * a, function (j) { if (self.grainOf[j] >= 0) crowd = true; });
      if (crowd) continue;
      const W = P.W, H = P.H, wall = Math.min(x[i], W - x[i], H - y[i]);
      const onGlass = wall < P.CRYSTAL_GLASS_RANGE;
      if (rng() > (onGlass ? P.CRYSTAL_NUC_GLASS : P.CRYSTAL_NUC_BULK)) continue;
      // orientation: a lattice row along the nearest glass, or the molecule's own frame
      let phi;
      const floor = H - y[i] <= Math.min(x[i], W - x[i]);
      if (kind === 1) {
        if (onGlass) phi = floor ? -Math.PI / 2 : 0; // one bond across the wall: zig-zag rows along it
        else { const h1 = this.bp[i * P.MAXB], h2 = this.bp[i * P.MAXB + 1]; phi = Math.atan2(y[h1] + y[h2] - 2 * y[i], x[h1] + x[h2] - 2 * x[i]) + Math.PI; }
      } else phi = onGlass ? (floor ? 0 : Math.PI / 2) : rng() * Math.PI / 3;
      // the template frame is the world frame at the moment of nucleation, the lattice origin on
      // the nucleus
      const g = { kind: kind, z: this.el[i], a: a, phi: phi, e: [], mem: [], sites: new Map(), ox: x[i], oy: y[i], cx: x[i], cy: y[i], sx: x[i], sy: y[i], th: 0 };
      for (let k = 0; k < 3; k++) g.e.push(Math.cos(phi + k * 2 * Math.PI / 3), Math.sin(phi + k * 2 * Math.PI / 3));
      this.grains.push(g);
      const site = { key: kind === 1 ? '0,0,0' : '0,0', x: x[i], y: y[i], sub: 0 };
      const mem = { atoms: [i], rest: [site.x, site.y], key: site.key };
      g.mem.push(mem); g.sites.set(site.key, mem);
      this.grainOf[i] = this.grains.length - 1;
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
        if (dx * dx + dy * dy > 9 * rest * rest) {
          // the bond is broken, so its energy is paid from the surroundings like a dissociation
          const be = C.bondEnergy(el[i], el[j], this.bo[i * MB + k]);
          this.setBond(i, j, 0);
          this.heat(i, j, -be);
          k--;
        }
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
      if (this.grab) {
        if (this.substep % 4 === 1) this.updateCarried();
        this.grabForces();
      }
      this.integrate();
      this.solveBonds(P.SHAKE_ITERS);
      this.rampBonds();
      this.walls();
      this.gasConduction();
      this.react();
      if (this.substep % 4 === 0) {
        this.bondEvents(P.DT * 4);
        this.breakStretched();
        this.solidAnneal();
      }
      if (this.substep % 64 === 0) this.updateMeltPoints();
      if (this.substep % 4 === 2) this.crystalUpdate();
      this.feedGas(P.DT);
      this.time += P.DT;
      this.substep++;
      this.impulseTime += P.DT;
      this.reactionWindowTime += P.DT;
    }
    this.updateGauges(substeps * P.DT);
  };

  // ---------------------------------------------------------------------------------------
  // Grab: the user picks up the material under the cursor and drags it. The hand holds only the
  // patch under the cursor; the grabbed atoms share one uniform pull (a critically damped spring
  // from their centre of mass to the hand, capped), so their own thermal motion is left alone and
  // nothing is heated or frozen by holding it. Everything else follows only through the forces
  // its atoms exert on each other: a solid comes along as a lump because its atoms attract one
  // another strongly, a liquid stretches, necks and drips, a gas slips away.
  // ---------------------------------------------------------------------------------------
  S.grabStart = function (x, y, r) {
    this.buildGrid();
    const atoms = [];
    const self = this;
    this.query(x, y, r, function (j) {
      const dx = self.x[j] - x, dy = self.y[j] - y;
      if (dx * dx + dy * dy < r * r) atoms.push(j);
    });
    if (!atoms.length) { this.grab = null; return 0; }
    let cx = 0, cy = 0;
    for (const j of atoms) { cx += this.x[j]; cy += this.y[j]; }
    cx /= atoms.length; cy /= atoms.length;
    // the hand holds the patch where it was taken: target = hand + (centre - hand at the start)
    this.grab = { atoms: atoms, body: [], hx: x, hy: y, vhx: 0, vhy: 0, ox: cx - x, oy: cy - y, t: this.time, vx: 0, vy: 0, bvx: 0, bvy: 0 };
    this.held.fill(0);
    for (const j of atoms) this.held[j] = 1;
    const b = this.grabBody(); this.grab.ox = b[1] - x; this.grab.oy = b[2] - y;
    this.updateCarried();
    return atoms.length;
  };
  S.grabMove = function (x, y) {
    const g = this.grab;
    if (!g) return;
    const dtw = Math.max(1e-3, this.time - g.t);
    // the hand's own velocity (smoothed), so the spring damps motion relative to the hand
    g.vhx += 0.5 * ((x - g.hx) / dtw - g.vhx); g.vhy += 0.5 * ((y - g.hy) / dtw - g.vhy);
    g.hx = x; g.hy = y; g.t = this.time;
  };
  S.grabEnd = function () { this.grab = null; this.held.fill(0); this.carried.fill(0); };
  // The solid body the held patch belongs to: solid atoms in contact with it, and with each other.
  // (A liquid is not carried: it hangs from the hand by its own cohesion, and drips.)
  S.updateCarried = function () {
    const n = this.n, carried = this.carried, el = this.el, x = this.x, y = this.y, self = this;
    const before = this.grabBody();
    carried.fill(0, 0, n);
    const queue = [], body = [];
    for (const i of this.grab.atoms) if (this.isSolidAtom(i)) { carried[i] = 1; queue.push(i); }
    while (queue.length) {
      const i = queue.pop();
      body.push(i);
      this.query(x[i], y[i], this.cellSize, function (j, d) {
        if (carried[j] || d > 1.3 * pairSigma[el[i] * NZ + el[j]] || !self.isSolidAtom(j)) return;
        carried[j] = 1; queue.push(j);
      });
    }
    this.grab.body = body;
    // the hand holds the body where it is: a change in what hangs from it is not a pull
    const after = this.grabBody();
    if (before[0] && after[0]) { this.grab.ox += after[1] - before[1]; this.grab.oy += after[2] - before[2]; }
  };
  // The hand steers what it holds (critically damped, at most GRAB_MAX_ACCEL and GRAB_MAX_SPEED).
  // A solid it holds moves as a whole: the push is shared over the body's atoms like its
  // weight, as it effectively is for anything this small (a real hand accelerates a lump a few
  // nanometres across far too gently to strain its bonds; on this scale its motion is ~1e13
  // times harsher, and pulling a lump along by one patch would stretch it plastically and heat
  // it). The push depends only on the motion of the body as a whole, so it neither heats nor
  // cools it, nor stops its atoms moving; and once let go, or when it hits something, nothing but
  // its own bonds holds it together. A liquid has no such body: the hand drags the patch alone,
  // and the liquid comes along only as far as its own cohesion carries it.
  S.grabForces = function () {
    const g = this.grab, A = g.atoms, x = this.x, y = this.y, px = this.px, py = this.py, el = this.el;
    if (this.time - g.t > 0.05) { const f = Math.pow(0.9, P.DT / P.DT_REF); g.vhx *= f; g.vhy *= f; } // the hand has stopped
    let M = 0, vx = 0, vy = 0;
    for (let k = 0; k < A.length; k++) {
      const i = A[k], m = mass[el[i]];
      M += m; vx += m * (x[i] - px[i]); vy += m * (y[i] - py[i]);
    }
    if (!M) return;
    g.vx = vx / M; g.vy = vy / M; // the held patch's velocity, per substep (for the thermostats)
    const b = this.grabBody();
    g.bvx = b[3]; g.bvy = b[4]; // ...and the body's
    const bvx = b[3] / P.DT, bvy = b[4] / P.DT;
    const w = P.GRAB_OMEGA;
    let ax = w * w * (g.hx + g.ox - b[1]) - 2 * w * (bvx - g.vhx);
    let ay = w * w * (g.hy + g.oy - b[2]) - 2 * w * (bvy - g.vhy);
    const nvx = bvx + ax * P.DT, nvy = bvy + ay * P.DT, nv = Math.hypot(nvx, nvy);
    if (nv > P.GRAB_MAX_SPEED) {
      const f = P.GRAB_MAX_SPEED / nv;
      ax = (nvx * f - bvx) / P.DT; ay = (nvy * f - bvy) / P.DT;
    }
    const a = Math.hypot(ax, ay);
    if (a > P.GRAB_MAX_ACCEL) { ax *= P.GRAB_MAX_ACCEL / a; ay *= P.GRAB_MAX_ACCEL / a; }
    const set = g.body.length ? g.body : A;
    for (let k = 0; k < set.length; k++) {
      const i = set[k];
      if (i >= this.n) continue;
      const m = mass[el[i]];
      this.fx[i] += m * ax; this.fy[i] += m * ay;
    }
  };
  // mass, centre and velocity (per substep) of what the hand steers: the body, or the patch
  S.grabBody = function () {
    const g = this.grab, set = g.body.length ? g.body : g.atoms, el = this.el;
    let M = 0, cx = 0, cy = 0, vx = 0, vy = 0;
    for (const i of set) {
      if (i >= this.n) continue;
      const m = mass[el[i]];
      M += m; cx += m * this.x[i]; cy += m * this.y[i]; vx += m * (this.x[i] - this.px[i]); vy += m * (this.y[i] - this.py[i]);
    }
    return M ? [M, cx / M, cy / M, vx / M, vy / M] : [0, 0, 0, 0, 0];
  };
  // particle `from` moved into slot `to` (to is being removed and from is the last one)
  S.grabRenumber = function (to, from) {
    this.held[to] = this.held[from]; this.held[from] = 0;
    this.carried[to] = this.carried[from]; this.carried[from] = 0;
    for (const A of [this.grab.atoms, this.grab.body]) {
      const k = A.indexOf(to);
      if (k >= 0) { A[k] = A[A.length - 1]; A.pop(); }
      const m = A.indexOf(from);
      if (m >= 0) A[m] = to;
    }
    if (!this.grab.atoms.length) this.grabEnd();
  };


  S.updateGauges = function (elapsed) {
    if (this.vacuumOn && elapsed > 0) this.pump(elapsed);
    if (this._keCount > 0) {
      const T = this._keSum / this._keCount;
      this.temperature = this.n ? T : this.plateTemp;
      this._keSum = 0; this._keCount = 0;
    }
    // Pressure of the gas in the jar right now, from its contents: P = n R T / V, with n the
    // gas-phase molecules (not part of a liquid or solid) and T from their own motion.
    const atm = this.gasPressure();
    this.pressureAtm = this.n ? this.pressureAtm * 0.5 + atm * 0.5 : 0;
    if (this.reactionWindowTime >= 0.5) {
      this.reactionsPerSec = this.reactionWindow / this.reactionWindowTime;
      this.reactionWindow = 0; this.reactionWindowTime = 0;
    }
  };

  // Vacuum pump: removes gas-phase molecules from the whole jar at pumpRate (fraction per
  // simulated second), so the pressure falls exponentially like a real pumped vessel. Liquids
  // are untouched until they evaporate (so water boils away under vacuum).
  S.pump = function (elapsed) {
    const n = this.n, MB = P.MAXB;
    if (!n) return;
    const pRemove = 1 - Math.exp(-this.pumpRate * elapsed);
    const seen = this._pmSeen && this._pmSeen.length >= P.MAX ? this._pmSeen : (this._pmSeen = new Uint32Array(P.MAX));
    const stamp = this._pmStamp = (this._pmStamp || 0) + 1;
    const queue = this._pmQueue || (this._pmQueue = new Int32Array(P.MAX));
    const doomed = [];
    for (let i = 0; i < n; i++) {
      if (seen[i] === stamp) continue;
      let head = 0, tail = 0, cond = 0;
      queue[tail++] = i; seen[i] = stamp;
      while (head < tail) {
        const a = queue[head++];
        cond += this.condF[a];
        for (let k = 0; k < this.nb[a]; k++) {
          const b = this.bp[a * MB + k];
          if (seen[b] !== stamp) { seen[b] = stamp; queue[tail++] = b; }
        }
      }
      if (tail > 12 || cond / tail >= 0.5 || this.touchesCondensed(queue, tail)) continue;
      if (this.rng() < pRemove) for (let k = 0; k < tail; k++) doomed.push(queue[k]);
    }
    if (!doomed.length) return;
    doomed.sort(function (a, b) { return b - a; });
    for (let k = 0; k < doomed.length; k++) this.removeParticle(doomed[k]);
    this.removedMol += doomed.length * P.MOL_PER_PARTICLE;
  };

  // Is any of the atoms list[0..count) touching a liquid or solid (not bonded to it)? A surface atom
  // of a lump can have too few contacts to count as packed, yet it is no part of the vapour.
  S.touchesCondensed = function (list, count) {
    const x = this.x, y = this.y, el = this.el, self = this;
    let hit = false;
    for (let k = 0; k < count && !hit; k++) {
      const a = list[k];
      this.query(x[a], y[a], 1.25 * (radius[el[a]] + 26), function (j) {
        if (hit || j === a || (self.condF[j] < 0.5 && self.grainOf[j] < 0) || self.isBonded(a, j)) return;
        const r = 1.25 * pairSigma[el[a] * NZ + el[j]];
        if ((x[j] - x[a]) ** 2 + (y[j] - y[a]) ** 2 < r * r) hit = true;
      });
    }
    return hit;
  };

  // Instantaneous gas pressure (atm) from the jar's contents.
  S.gasPressure = function () {
    const n = this.n, MB = P.MAXB;
    if (!n) return 0;
    const seen = this._gpSeen && this._gpSeen.length >= P.MAX ? this._gpSeen : (this._gpSeen = new Uint32Array(P.MAX));
    const stamp = this._gpStamp = (this._gpStamp || 0) + 1;
    const queue = this._gpQueue || (this._gpQueue = new Int32Array(P.MAX));
    const x = this.x, y = this.y, px = this.px, py = this.py, el = this.el, dt = P.DT;
    let molecules = 0, ke = 0;
    for (let i = 0; i < n; i++) {
      if (seen[i] === stamp) continue;
      // walk the molecule containing i
      let head = 0, tail = 0, m = 0, mvx = 0, mvy = 0, cond = 0;
      queue[tail++] = i; seen[i] = stamp;
      while (head < tail) {
        const a = queue[head++];
        const ma = mass[el[a]];
        m += ma; mvx += ma * (x[a] - px[a]); mvy += ma * (y[a] - py[a]);
        cond += this.condF[a];
        for (let k = 0; k < this.nb[a]; k++) {
          const b = this.bp[a * MB + k];
          if (seen[b] !== stamp) { seen[b] = stamp; queue[tail++] = b; }
        }
      }
      if (tail > 12 || cond / tail >= 0.5) continue; // liquid, solid or network: not gas
      if (this.touchesCondensed(queue, tail)) continue; // on the surface of a liquid or solid
      molecules++;
      ke += (mvx * mvx + mvy * mvy) / (2 * m * dt * dt); // translational KE of the molecule
    }
    if (!molecules) return 0;
    const T = ke / (molecules * P.KB); // 2D: translational KE = kT per molecule
    return molecules * P.MOL_PER_PARTICLE * 0.082057 * T / P.JAR_VOLUME_L;
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
    const species = new Map(), tLoc = this.tLoc;
    const descCache = this._descCache || (this._descCache = new Map());
    mols.forEach(function (m) {
      const sig = JSON.stringify(m.counts) + (m.ionic ? 'i' : '');
      let d = descCache.get(sig);
      if (!d) { d = C.describe(m.counts, m.atoms.length, m.ionic); descCache.set(sig, d); }
      let s = species.get(d.key);
      if (!s) {
        s = { key: d.key, formula: d.formula, name: d.name, molecules: 0, atoms: 0, condensed: 0, color: speciesColor(d.key, m.counts), first: m.atoms[0], counts: m.counts, ionic: m.ionic, tSum: 0 };
        species.set(d.key, s);
      }
      s.molecules++;
      s.atoms += m.atoms.length;
      s.tSum += tLoc[m.atoms[0]];
      // networks/lattices are solids; small molecules are condensed when packed against others
      if (m.atoms.length > 12 || m.contacts / m.atoms.length >= 1) s.condensed++;
      m.species = s;
    });
    const molColor = this.molColor;
    mols.forEach(function (m) {
      const col = m.species.color;
      for (let k = 0; k < m.atoms.length; k++) molColor[m.atoms[k]] = col;
    });
    const list = Array.from(species.values());
    list.forEach(function (s) {
      // a "molecule" of a network/lattice counts its formula units; when the formula is no real
      // formula unit (graphite with a few oxygens on its edge reduces to C150O) it counts the main
      // element's atoms instead, so a lump does not vanish from the list when something sticks to it
      let units = s.molecules;
      if (s.key.slice(-2) === ')n') {
        const size = formulaSize(s.key);
        if (size <= 12) units = s.atoms / size;
        else {
          let total = 0, most = 0;
          for (const z in s.counts) { total += s.counts[z]; most = Math.max(most, s.counts[z]); }
          units = s.atoms * most / total;
        }
      }
      s.mol = units * P.MOL_PER_PARTICLE;
      s.phase = s.condensed / s.molecules > 0.5 ? 'condensed' : 'gas';
    });
    list.sort(function (a, b) { return b.mol - a.mol; });
    // per-atom species index and each species' real-world look (for the real-world view)
    const spIdx = this.spIdx || (this.spIdx = new Int32Array(P.MAX));
    list.forEach(function (s, k) {
      s.idx = k;
      s.T = s.tSum / s.molecules;
      s.look = CE.look ? CE.look.lookOf(s.key, s.counts, s.atoms / s.molecules, s.ionic, s.T) : null;
    });
    mols.forEach(function (m) {
      for (let k = 0; k < m.atoms.length; k++) spIdx[m.atoms[k]] = m.species.idx;
    });
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
