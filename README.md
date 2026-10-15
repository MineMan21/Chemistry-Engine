# Chemistry Engine

A particle chemistry sandbox that runs in the browser. Pick any of the 118 elements, drop it into a
jar, and watch atoms bond into molecules. Heat the jar on a hot plate, pump it down with a vacuum
pump, or pressurize it from a gas can.

- **Zoomed in (atom view):** ball-and-stick atoms with single, double, triple and ionic bonds,
  element symbols, partial charges (δ+/δ−) and ion charges. Clicking drops individual atoms.
- **Zoomed out (bulk view):** every atom becomes a pixel, coloured by the molecule it belongs to,
  and the simulation keeps running exactly the same way. Clicking adds material by the mole.

## Running it

Open `index.html` in a modern browser. There is no build step and no dependency; it also works
straight from the file system.

Tests run under Node (no packages needed):

```sh
npm test
```

## Controls

| Action | How |
| --- | --- |
| Zoom | Mouse wheel, pinch, `+` / `−`, or the buttons at the top right (`F` fits the jar) |
| Pan | Drag in the jar |
| Pick an element | Click it in the periodic table |
| Add atoms | Atom view: click in the jar (1, 5, 20 or 100 per click) |
| Add moles | Bulk view: click in the jar, or type an amount and press **Add** |
| Spark | Spark tool, then click: a short 4000 K burst that also splits some bonds (makes radicals) |
| Chill | Chill tool, then click: cools a spot to 80 K |
| Inspect | Hover any atom; with the Inspect tool, click to pin a molecule. Clicking a row in **Contents** pins that species |
| Hot plate | Slider from 20 K to 4000 K (logarithmic) plus presets |
| Pressure | Vacuum pump toggle; gas can toggle with gas choice and flow rate |
| Pause | Space |

The **Experiments** menu loads ready-made setups: hydrogen + oxygen, sodium + chlorine, liquid
water, methane + oxygen, iron in hot oxygen, carbon + oxygen, noble gases and liquid nitrogen.

## Scale

The jar is a 32 × 40 nm slice, and each particle statistically stands for **1 mmol** of atoms
(6.0 × 10²⁰), so **1 mol = 1000 particles**. That is what lets you add moles at a time and still
run in real time. The jar is treated as 24 L, so 1 mol of gas at room temperature reads about
1 atm on the gauge. Diatomic elements (H, N, O, F, Cl, Br, I) are added as X₂ by default, so
"1 mol of hydrogen" means 1 mol of H₂.

## How the chemistry works

Exact quantum chemistry for every combination of 118 elements is not computable in real time,
so the engine uses rules driven by real element data (`js/elements.js`: mass, Pauling
electronegativity, covalent radius, common valences, melting and boiling points).

- **Bond type.** Two metals do not form discrete bonds; they attract each other with metallic
  cohesion scaled to the metal's boiling point. A metal with an electronegativity gap of at
  least 1.7 forms ionic bonds. Everything else bonds covalently. Noble gases only bond to more
  electronegative partners when they have valences beyond 0 (Xe, Kr, Rn).
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
  - four-centre exchange AB + CD → AC + BD (largely forbidden, so it needs high temperature),
  - thermal dissociation of the weakest bond unit, and pi bonds forming between radicals.

  Activation energy follows an Evans–Polanyi-style rule: a fraction of the energy of the bonds
  being broken plus any endothermicity. Triplet O₂ gets an extra spin barrier against
  closed-shell molecules and metal surfaces, which is why fuels and iron sit happily in oxygen
  until you heat or spark them, while radicals react with it readily.
- **Energy is conserved.** Reaction enthalpy goes into (or comes out of) the atoms' kinetic
  energy, so exothermic reactions heat the jar and flames are hot. Changes in the atoms'
  non-bonded interactions during a reaction are accounted for too. When the colliding atoms
  cannot pay for an endothermic step, the energy is borrowed from the surrounding region.
- **Geometry.** VSEPR sets the angle around atoms with two bonds (water 104.5°, H₂S 92°, CO₂
  180°); atoms with more bonds spread their neighbours evenly. The inspector shows the real 3D
  shape and angle; the drawing is a 2D projection.
- **Phases** emerge from non-bonded forces: van der Waals attraction, partial-charge
  electrostatics, explicit hydrogen bonds (H on N/O/F to N/O/F), ionic Coulomb forces and
  metallic cohesion. These are calibrated so that water is liquid at room temperature and boils
  off above about 373 K, Br₂ is liquid, I₂ and metals are solid, and N₂, O₂, Cl₂ and Ar are gases
  that condense near their real boiling points.
- **Apparatus.** The hot plate exchanges heat with atoms near the bottom (plus gentle convection
  through the jar); the glass loses heat to a 293 K room unless the jar is insulated. Pressure is
  measured from momentum transfer to the lid and upper walls, like a real gauge. The vacuum pump
  removes whole molecules that reach the valve; the gas can injects gas through it.

### Approximations and known limits

- It is a 2D model: geometry is flattened and there are fewer degrees of freedom than in 3D,
  so heat capacities and some rates are only roughly right.
- Each particle is a statistical stand-in for 1 mmol, and every simulated collision represents
  many real ones (folded into the Arrhenius prefactor), so time is compressed and reaction
  speeds are qualitative. Gas-phase combustion started by a spark spreads as a burn over tens of
  seconds rather than an instant explosion; a hot plate at 1500–2000 K burns a mixture faster.
- Resonance, aromaticity, dative bonds, coordination complexes, catalysis, solvation of ions,
  isomers with the same formula and radioactive decay are not modelled.
- Superheavy elements have estimated electronegativities (from the element above them) and
  predicted or unknown melting and boiling points.
- Hydrogen uses a heavier dynamical mass (3 u) to keep the integrator stable; all chemistry
  uses the real mass.

## Code layout

| File | Contents |
| --- | --- |
| `js/elements.js` | Data for all 118 elements, periodic table layout and categories |
| `js/chemistry.js` | Bond classes, bond energies and lengths, valence rules, VSEPR, formula naming |
| `js/sim.js` | The particle simulation: forces, constraints, reactions, thermostats, pump and can |
| `js/render.js` | Canvas drawing of the apparatus, atom view and pixel bulk view |
| `js/ui.js` | Periodic table, element card, contents list and inspector |
| `js/main.js` | Input handling, controls, presets and the main loop |
| `tests/` | Headless tests (`npm test`) |

Tunable constants (time step, Arrhenius prefactors, calibration factors) sit at the top of
`js/sim.js` in the `P` object, each with a comment.
