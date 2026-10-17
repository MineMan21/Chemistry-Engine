# Chemistry Engine

A particle chemistry sandbox that runs in the browser. Pick any of the 118 elements, drop it into a
jar, and watch atoms bond into molecules. Heat the jar on a hot plate, pump it down with a vacuum
pump, or pressurize it from a gas can.

There are three zoom levels over the same running simulation (switch with the **Atoms / Moles /
Real world** buttons at the top left, or just scroll):

- **Atoms:** ball-and-stick atoms with single, double, triple and ionic bonds, element symbols,
  partial charges (δ+/δ−) and ion charges. Clicking drops individual atoms.
- **Moles:** every atom becomes a dot (1 mmol each), coloured by the molecule it belongs to.
  Clicking adds material by the mole.
- **Real world:** the jar as you would see it: water as a translucent blue pool with a bright
  surface, metals as grey lumps, sulfur yellow, salt white, rust red-brown, and coloured gases
  (chlorine, bromine, iodine, NO₂) as a tinted haze, while colourless gases stay clear. Hot
  liquids and solids glow: dull red from about 800 K, orange to yellow-white when molten. Clicking
  adds moles here too.

Liquids fall under gravity, land in a pool and spread out; solids (metals below their melting
point, salts, covalent networks) fall as lumps and keep their shape until they melt. Glass attracts polar
molecules, so water climbs into a meniscus, thin films and drops cling to the walls, and vapour
condenses on the room-temperature glass of a hot jar.

## Running it

Open `index.html` in a modern browser. There is no build step and no dependency; it also works
straight from the file system.

Tests run under Node (no packages needed):

```sh
npm test
TEST_FILTER=iron npm test   # only the tests whose name contains "iron"
```

## Controls

| Action | How |
| --- | --- |
| Zoom | Mouse wheel, pinch, `+` / `−`, the Atoms / Moles / Real world buttons, or the buttons at the top right (`F` fits the jar) |
| Pan | Drag in the jar |
| Pick an element | Click it in the periodic table |
| Pick a molecule | Click a chip under **Molecules**: presets (water, CO₂, methane, ammonia, CO, HCl, salt, NaOH, H₂O₂, ethanol, benzene, sulfuric acid) or your own saved ones |
| Save a molecule | Choose **Inspect**, click any molecule in the jar, then press **Save this molecule** (or `S`). It appears under **My molecules**, kept in this browser; × deletes it |
| Add atoms or molecules | Atom view: click in the jar (1, 5, 20 or 100 per click) |
| Add moles | Moles or Real world view: click in the jar, or type an amount and press **Add** |
| Spark | Spark tool, then click: a short 4000 K burst over a small spot that also splits some bonds, like an igniter. Use it to light fuels such as H₂ + O₂ or CH₄ + O₂ |
| Chill | Chill tool, then click: blasts a small spot down to 80 K, to freeze or condense what is there |
| Inspect | Hover any atom; with the Inspect tool, click to pin a molecule |
| Grab | Grab tool, then press on material and drag: the hand takes hold of what is under the cursor. A solid (a lump of metal, salt, graphite or ice) is carried as a whole, its weight and the hand's push shared over all its atoms, so it is neither strained nor warmed by being moved; let go or thrown at the glass, only its own bonds hold it together. A liquid is dragged by the patch under the cursor alone, so it stretches, necks and drips; a gas slips away |
| Highlight | Click a row in **Contents** to highlight every molecule of that substance in the jar (others fade; in the real-world view each one gets a yellow dot). Click it again or press Esc to clear |
| Phase diagram | **Phase diagram** under the element or molecule card, or under Contents for a highlighted substance: where it is solid, liquid, gas or supercritical, with the triple and critical points and a ring showing the jar's current temperature and pressure. Hover to read values. Available for about 45 well-documented substances (water, CO₂, N₂, O₂, the halogens and noble gases, mercury, ammonia, methane, simple organics...) |
| Hot plate | Slider from 20 K to 4000 K (logarithmic) plus presets |
| Pressure | Vacuum pump toggle with a strength slider (2–95% of the gas per second); gas can toggle with gas choice and flow rate |
| Pause | Space |

The **Experiments** menu loads ready-made setups: hydrogen + oxygen, sodium + chlorine, a blob
of water dropped from mid-air, steam on cold glass, methane + oxygen, iron in hot oxygen, carbon +
oxygen, noble gases and liquid nitrogen, plus three textbook phase-change experiments that run
the hot plate through a program: freezing water into hexagonal ice, melting an iron lump at
2100 K and recasting it as it cools, and melting gallium just above room temperature. Touching
the plate control takes over from the program.

## Scale

The simulated jar is a 24 × 32 nm slice, and each particle statistically stands for **1 mmol** of atoms
(6.0 × 10²⁰), so **1 mol = 1000 particles**. That is what lets you add moles at a time and still
run in real time. The jar is treated as 24 L, so 1 mol of gas at room temperature reads about
1 atm on the gauge. The real-world view draws that sample as a jar 24 cm wide (its scale bar is
in cm), while the atom and mole views show nm and Å. Diatomic elements (H, N, O, F, Cl, Br, I) are added as X₂ by default, so
"1 mol of hydrogen" means 1 mol of H₂.

## How the chemistry works

Exact quantum chemistry for every combination of 118 elements is not computable in real time,
so the engine uses rules driven by real element data (`js/elements.js`: mass, Pauling
electronegativity, covalent radius, common valences, melting and boiling points).

- **Bond type.** Two metals do not form discrete bonds; they hold together by metallic
  cohesion (see Phases below). A metal with an electronegativity gap of at least 1.5 to its
  partner forms ionic bonds (the oxides, halides and sulfides of the metals: MgO, Fe₂O₃, CuO,
  NaCl). Everything else bonds covalently. Noble gases only bond to more electronegative partners
  when they have valences beyond 0 (Xe, Kr, Rn). Only strongly electropositive metals form
  hydrides (LiH, NaH, CaH₂, MgH₂); iron, copper and the other transition metals do not.
- **Ionic compounds are lattices.** Each ion–ion contact carries one unit of charge, so a doubly
  or triply charged ion bonds to several counter-ions and the compound grows as an extended
  lattice with its textbook ratio ((MgO)n, (Al₂O₃)n, iron's Fe₂O₃ and Fe₃O₄) instead of closing
  into small molecules. Ions touch at ionic radii (cations shrink, anions swell: Na⁺ 102 pm,
  Cl⁻ 181 pm) and feel each other's full charges over a longer range, so salts and oxides are
  solids.
- **Valence.** Each element has allowed valences. Going beyond the normal one (for example
  SF₆ or PCl₅) needs a more electronegative partner.
- **Bond energies.** Measured average bond enthalpies for about 100 common bonds (H–O 459,
  N≡N 942, C=O 799 kJ/mol …); everything else comes from Pauling's rule,
  D(A–B) = ½[D(A–A) + D(B–B)] + 96.5·Δχ². Bond lengths come from covalent radii, shortened for
  double and triple bonds.
- **Reactions** happen when atoms collide, with Arrhenius probabilities at the local
  temperature:
  - radical + radical → bond (almost no barrier),
  - radical abstraction (a free valence steals a bond unit, e.g. Na + Cl₂ → NaCl + Cl),
  - atom transfer in hot zones (H + O₂ → OH + O), which branches flame chains,
  - a diatomic molecule splitting across a metal surface (2 M + O₂ → 2 M–O, 2 Na + Cl₂ → 2 NaCl),
    which is how metals oxidise: no stray O atoms or peroxides,
  - alkali and alkaline-earth metals react by electron transfer, with almost no barrier
    beyond the reaction enthalpy (sodium in water),
  - four-centre exchange AB + CD → AC + BD (largely forbidden, so it needs high temperature),
    and its one-bond variant AB + CD → AC + B + D, an abstraction by a molecule, with the
    abstraction barrier (it starts the dark explosion of H₂ + F₂; with Cl₂ or O₂ it is far uphill),
  - thermal dissociation of the weakest bond unit, and pi bonds forming between radicals.

  Carbon monoxide is a closed-shell molecule (its carbon's spare valence is a lone pair, not an
  unpaired electron), so CO does not stick to other CO; radicals still add to it (O + CO → CO₂).
  Oxygen does not grow O–O–O chains on molecules that hold other elements (H₂O₃ and longer
  polyoxides fall apart far below room temperature); pure oxygen still makes ozone.

  Activation energy follows an Evans–Polanyi-style rule: a fraction of the energy of the bonds
  being broken plus any endothermicity. Triplet O₂ gets an extra spin barrier against
  closed-shell molecules and metal surfaces, which is why fuels and iron sit happily in oxygen
  until you heat or spark them, while radicals react with it readily. A halogen splitting on a
  metal that does not react by electron transfer (iron, copper, aluminium) meets an extra barrier
  for breaking through its own halide layer, so iron burns in chlorine only once heated (about
  500 K), while sodium reacts at once.
- **Energy is conserved.** Reaction enthalpy goes into (or comes out of) the atoms' kinetic
  energy, so exothermic reactions heat the jar and flames are hot. Changes in the atoms'
  non-bonded interactions during a reaction are accounted for too. The heat is shared with the
  molecules around the reacting pair (the third bodies that carry it off in a real
  recombination). When the colliding atoms cannot pay for an endothermic step, the energy is
  borrowed from the surrounding region.
- **Geometry.** VSEPR sets the angle around atoms with two bonds (water 104.5°, H₂S 92°, CO₂
  180°); atoms with more bonds spread their neighbours evenly. The inspector shows the real 3D
  shape and angle; the drawing is a 2D projection.
- **Phases** emerge from non-bonded forces: van der Waals attraction, partial-charge
  electrostatics, explicit hydrogen bonds (H on N/O/F to N/O/F), ionic Coulomb forces and
  metallic cohesion. These are calibrated so that water is liquid at room temperature and boils
  off above about 373 K, Br₂ is liquid, I₂ and metals are solid, and N₂, O₂, Cl₂ and Ar are gases
  that condense near their real boiling points. Hydrogen bonds are directional: strongest when the
  O–H points straight at the neighbouring O, which is what lets water freeze into an ordered,
  open lattice rather than a dense heap. Water oxygens also feel a gentle three-body term after
  the coarse-grained mW water model (Molinero & Moore, 2009) that leans toward 120° O–O–O angles,
  the hexagonal rings of ice Ih seen down its c axis; it is kept weak so liquid water stays a
  disordered, mobile jumble, and the ice honeycomb comes from the crystal engine below. Water also
  gains a many-body cohesion only when a molecule is surrounded by others (in the liquid or ice),
  never as a pair or a small clump, so steam stays a gas of separate molecules rather than drifting
  droplets (in 2D, molecules meet far more often than in 3D, so pairwise stickiness alone would
  cluster the vapour).
- **Boiling points are calibrated per substance.** The generic rules (van der Waals wells from
  each element, hydrogen bonds, partial charges) put some substances' boiling points far off, so
  documented molecular substances carry a cohesion factor on their van der Waals and hydrogen-
  bond attraction: N₂, O₂, F₂, Cl₂, Br₂, I₂, CO₂, CH₄, CO, HCl, NH₃, ethanol, benzene and H₂O₂.
  Each was bisected so that a fixed amount in the jar is half evaporated at the temperature where
  the real substance (by its documented vapour pressure) would be. When a reaction turns a
  molecule into another substance, the change in its attraction is settled as heat.
- **Metallic cohesion** has two parts. A pair well scaled to the metal's melting point sets when
  the close-packed solid melts. A many-body embedding term (in the spirit of Finnis–Sinclair and
  embedded-atom models) rewards each atom for sitting inside dense metal, however that metal is
  packed: crystal and melt gain the same, so it leaves the melting point alone, but a lone atom,
  a small cluster or a froth gains nothing. That is what keeps a liquid metal together far above
  its melting point; a pair well alone would boil a 2D metal almost as soon as it melted. The two
  together give each metal its real cohesive energy (iron 413 kJ/mol, lead 196, tungsten 859),
  so pulling a lump apart costs what it costs the real metal. Iron
  is solid at 1600 K and a liquid with almost no vapour at 2000 K, lead melts between 500 and
  650 K, mercury is liquid at room temperature and gallium at 340 K. Melting in 2D is gradual,
  and gallium and mercury freeze some way below their real melting points. The pair well alone
  decides where a metal melts: the many-body term acts only on atoms being pulled out of the
  metal (a corner atom, one leaving for the vapour, the last threads of a crack), never on the
  inside of a solid or a melt, so it cannot favour the crystal.
- **Crystal engine.** Freezing is a nucleation-and-growth process: below its melting point a
  substance that crystallises (water, free metal atoms) starts crystal grains, readily on the cold
  glass and floor and rarely in the bulk, and the grains grow as neighbours join. Each grain keeps
  an ideal lattice (the honeycomb of ice, the close-packed layer of a metal) fitted to its members
  by a best-fit rotation and shift (shape matching). In ice, gentle spring forces pull members
  toward their sites (every join, leave and spring is energy-accounted), and ice grains are
  spaced so that ice takes up about 7% more room than the liquid (real ice: 9%), so water swells
  as it freezes. Metal grains carry no springs: a metal orders and holds together by its own
  cohesion alone, and its grains only keep track of which atoms are crystalline. Salts and oxides
  are held by their ionic bonds, molecular solids by their own attractions, covalent networks
  such as graphite and silicon by their bonds. A member melts when it is hot itself, not merely
  when the gas around it is.
- **Crystals form on screen.** A real crystal reorganises over countless molecular vibrations
  that the time-compressed simulation cannot show one by one, so molecules in a frozen molecular
  solid (below its documented melting point) also make occasional small trial turns and shifts,
  accepted by the Boltzmann (Metropolis) rule at the local temperature with the energy change
  paid from local heat. That only speeds up reaching equilibrium: water cooled below 273 K, or
  frost growing from vapour, settles into ordered rows instead of a frozen glass. The same
  melting points mark frozen molecular solids as solids for the settling force, so ice is not
  squashed like a liquid.
- **Temperature** is the kinetic energy over the degrees of freedom actually left: each rigid
  bond removes one of an atom's two, and in a covalent network rings and stiff bond angles lock
  up about half of what remains, so a lump of graphite on a hot plate reads the plate's
  temperature (and can burn). Where liquid or solid shares a region with gas, the few gas
  molecules barely count, since they hold very little heat next to condensed matter.
- **Heat flow in gas.** Simulated gas molecules collide only a few times per simulated second,
  far too rarely to carry heat across the jar, so hot gas also passes heat to the cooler gas next
  to it (energy conserving). This is what lets a flame front run through a whole mixture.
- **Apparatus.** The hot plate exchanges heat with atoms near the bottom (plus gentle convection
  through the jar); the glass loses heat to a 293 K room unless the jar is insulated. The pressure
  gauge reads the gas in the jar directly each frame (P = nRT/V from the gas-phase molecules and
  their own motion), so it responds immediately when gas is added, heated, pumped out or removed. The vacuum pump
  draws gas out of the whole jar (about half of it per simulated second, so the pressure falls
  off exponentially); liquids stay until they evaporate. The gas can injects gas through the valve.
- **Gravity and glass.** At a few nanometres real gravity would do nothing, so the engine applies
  a strong settling force only to material that is clearly packed (liquids and solids); gas
  molecules feel almost none, so vapour still fills the jar. What the hand holds feels none:
  the hand carries its weight. The glass walls and floor attract
  nearby molecules (strongly for hydrogen-bonding ones like water, weakly for non-polar ones,
  barely for metals) and hold liquid touching them with friction, which gives the meniscus and
  clinging drops. Liquid touching the glass exchanges heat with it at room temperature.

### Approximations and known limits

- It is a 2D model: geometry is flattened and there are fewer degrees of freedom than in 3D,
  so heat capacities and some rates are only roughly right.
- Each particle is a statistical stand-in for 1 mmol, and every simulated collision represents
  many real ones (folded into the Arrhenius prefactor), so time is compressed and reaction
  speeds are qualitative. A sparked 2:1 H₂/O₂ mixture burns through in a second or two
  (a real one explodes in milliseconds); the flame peaks around 4000–5000 K, hotter than the
  real ~3000 K because 2D molecules have fewer ways to store the heat, then cools to the room.
- Resonance, aromaticity, dative bonds, coordination complexes, catalysis, solvation of ions,
  isomers with the same formula and radioactive decay are not modelled.
- Superheavy elements have estimated electronegativities (from the element above them) and
  predicted or unknown melting and boiling points.
- Metal cohesion is calibrated on the melting point, so metals boil in the right order but
  low-melting metals boil too early (gallium well below its real 2673 K), and a molten metal
  shows more vapour (a few tenths of an atmosphere for iron at 2100 K) than the real one.
- Evaporation rises too gently with temperature: in an insulated jar a small amount of water is
  about 80% liquid at 298 K (85% real), but still about two-thirds liquid at 373 K and half at
  450 K, where real water would be all vapour. Liquids boil over a wide range rather than at a
  sharp boiling point; a hot plate well above the boiling point drives them off.
- Reactions that the test suite checks and that work: the hydrogen flame, H₂ + Cl₂
  (spark) → HCl, H₂ + F₂ exploding without a spark, Na + Cl₂ → solid NaCl, metals burning to
  their oxides (MgO, Al₂O₃, iron's Fe₂O₃ and Fe₃O₄, CuO), iron in hot steam releasing hydrogen,
  iron burning in chlorine only when heated, gold untouched by hot oxygen, CO stable, hot
  hydrogen peroxide decomposing to water and oxygen, no ammonia from N₂ + H₂ without a catalyst.
  Still poorly reproduced: hydrocarbon flames stop at peroxyacids and polymers instead of going
  all the way to CO₂ (the high-temperature route through formaldehyde and CO is not modelled),
  sodium reacts with water only at its surface and slowly, and with no ions in solution there is
  no neutralisation and no NH₄Cl. Rusting at room temperature (electrochemical, over days) does
  not happen. See `Known_Issues.md` for the full list.
- Scale: moved at the speed of a real hand, nothing would visibly move in the time the simulation
  shows, so the hand moves things at a few times the thermal speed of their atoms, about 10^13
  times harsher than a real hand. That is why it carries a solid as a whole instead of pulling
  it along by the patch it touches (which would stretch and warm a lump that real handling never
  could). A solid thrown at the glass at the hand's top speed hits it at roughly the speed of a
  rifle bullet, so it warms (metals, salts, MgO and graphite stay whole; sulfur and iodine
  crumble). Silicon, almost fully rigid in 2D, reads too cold.
- Integration: a substep of 1 ms is short enough to follow a metal lattice's vibrations, so a
  lump left alone in an insulated jar keeps its temperature to a couple of kelvin per second.
- Hydrogen uses a heavier dynamical mass (3 u) to keep the integrator stable; all chemistry
  uses the real mass.

## Code layout

| File | Contents |
| --- | --- |
| `js/elements.js` | Data for all 118 elements, periodic table layout and categories |
| `js/chemistry.js` | Bond classes, bond energies and lengths, valence rules, VSEPR, formula naming |
| `js/sim.js` | The particle simulation: forces, constraints, reactions, thermostats, pump and can |
| `js/phase.js` | Triple/critical points and normal melting/boiling points for documented substances, and the phase-boundary curves built from them |
| `js/phasechart.js` | The phase diagram chart (SVG) and its table of key values |
| `js/render.js` | Canvas drawing of the apparatus, atom view and pixel bulk view |
| `js/ui.js` | Periodic table, element card, contents list and inspector |
| `js/main.js` | Input handling, controls, presets and the main loop |
| `tests/` | Headless tests (`npm test`) |

Tunable constants (time step, Arrhenius prefactors, calibration factors) sit at the top of
`js/sim.js` in the `P` object, each with a comment.
