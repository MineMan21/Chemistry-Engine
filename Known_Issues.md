# Known Issues: real-world accuracy

This document looks at how closely the sandbox's experiments match real chemistry. To find
these issues, I ran every preset under **Experiments**, plus about 30 extra setups (alkali
metals in water, acid–base reactions, combustion of other fuels, metals in oxygen and chlorine,
phase changes of common substances). Each setup ran headlessly through `tests/load.js` for 15–120
simulated seconds, and some also ran in the browser UI. Each result was compared with what the
real experiment does.

Issues that the README already admits under *Approximations and known limits* are marked
**(README)**. They are listed anyway because a preset or the UI still promises the real
behaviour.

Severity: **High** means the outcome is qualitatively wrong (wrong products, a reaction that
shouldn't happen does, or one that should happen doesn't). **Medium** means the outcome is right
but numbers or timing are far off. **Low** means a cosmetic or data issue.

Each issue now carries a **Status** from the fix pass that followed this review: **Fixed**,
**Improved** (better, not yet right) or **Open**. Regression tests for the fixes are in
`tests/sim.test.js`.

---

## Status after the fix pass

| # | Issue | Status |
| --- | --- | --- |
| G | Solids break when grabbed or thrown; moving them heats them | **Fixed** except iodine and sulfur: see the table below |
| 1 | CO polymerises at room temperature | **Fixed** |
| 2 | Hydrocarbon flames never reach CO₂ | Open |
| 3 | Ethanol decomposes at room temperature | **Improved**: ~95% intact after 15 s |
| 4 | Gold oxidises | **Fixed** |
| 5 | Iron + chlorine react at room temperature | **Fixed** |
| 6 | Sodium and potassium barely react with water | Open |
| 7 | NH₃ + HCl makes radicals | **Improved**: no radicals; still no NH₄Cl |
| 8 | Neutralisation makes radicals | **Improved**: no radicals; still no reaction |
| 9 | H₂ + F₂ does not react in the dark | **Fixed** |
| 10 | Peroxide builds H₂O₃–H₂O₅ | **Fixed** |
| 11 | Sulfur does not burn | **Fixed** |
| 12 | Metal oxides as small gas molecules | **Improved**: rust is a solid scale |
| 13 | Carbon burns at room temperature | **Fixed** |
| 14 | No NO from hot N₂ + O₂ | Open |
| 15 | Network solids never warm up | **Fixed** for graphite; silicon reads low |
| 16 | Gauge reads vapour that isn't there | **Fixed**; Hg and I₂ remain too volatile |
| 17 | Liquids don't boil above their boiling point | Open |
| 18 | The jar runs far below the plate | Open |
| 19 | Condensed matter at impossible temperatures | Open |
| 20 | CO₂ polymerises when cold | **Fixed** |
| 21 | Sulfur makes the wrong molecules | Open |
| 22 | A network with a stray atom vanishes from Contents | **Fixed** |
| 23 | Presets add less than they ask for | **Fixed** |
| 24 | No glow from hot materials | **Fixed** |
| 25 | Data inconsistencies | **Fixed** (except NaH's bond class) |

### Root causes found during the fix pass

Several of the issues above traced back to a handful of engine bugs rather than to chemistry
rules:

- **Overlap flash on adding molecules.** `addMolecules` spaced molecules by half an atom
  radius, so long molecules (CO₂, ethanol) started inside their neighbours, and hand-drawn
  templates (ethanol's) put atoms inside each other. Both flashed to 4,000–10,000 K the moment
  they were added. That flash, not room-temperature chemistry, made CO₂ polymerise (#20), ethanol
  decompose (#3) and NH₃ + HCl and acid + base make radicals (#7, #8). Templates are now
  relaxed once before use and molecules are placed without overlaps. (The old "wet iron rusts at
  350 K" test also passed only because of this flash; it now tests the real steam–iron reaction.)
- **Temperature of rigid networks.** The thermometer assumed every atom keeps at least one
  of its two degrees of freedom. Graphite keeps about 0.3, so it read a quarter of its real
  temperature and could not burn (#15, #11, #13).
- **Annealing moves turned motion into heat.** Trial turns of molecules in a frozen solid
  also rotated their bulk velocity, so any moving molecular solid (a lifted ice cube) heated
  itself and melted.
- **Bulk motion counted as heat.** Surface atoms of a moving lump, and stationary atoms
  sharing a region with it, read its speed as temperature.

---

## Grabbing and throwing solids

Reported issue: many solids were too weak to survive being grabbed or thrown, and moving
something with the hand warmed it up.

Each lump below (0.3 mol, settled on the floor) was grabbed by a patch at its top, lifted
700 units, then thrown at the glass at the hand's top speed. Each entry is the share of the lump
left in its largest piece after the lift and after the throw.

| Solid | Before: lift / throw | After: lift / throw |
| --- | --- | --- |
| Iron, copper, gold | whole / whole | whole / whole |
| Aluminium | 8% / torn | whole / whole |
| Lead | 12% / torn | whole / whole |
| Gallium (270 K) | 8% / torn | whole / 97% (dented) |
| Graphite | 81% / whole | whole / whole |
| Silicon | 0% (stuck) | whole / whole |
| Salt (NaCl) | 8% / torn | whole / whole |
| Magnesium oxide | 0% / torn | whole / whole |
| Ice (200 K) | 12% / torn | whole / whole |
| Sulfur | 4% / torn | crumbles (its rings form hot, above its melting point) |
| Iodine | 2% / torn | still fails |

Heating: a 300-atom lump of iron, lead or gold, lifted and held in mid-air in an insulated jar,
now changes temperature by about 1–3 K per second. An untouched lump drifts about the same.
Before, a held lump of iron gained tens of kelvin per second.

What changed:

- **Cohesion comes only from the atoms' own forces.** Each metal's attraction now adds up to its
  real cohesive energy (iron 413 kJ/mol, lead 196, from Kittel). It is a pair well plus a
  many-body embedding term, which makes an atom inside dense metal more tightly bound than one
  at a surface, crack or cluster. Metal crystals no longer have template springs: springs to a
  lattice template are not forces between atoms, and they leaked energy as the template refit.
  Nothing slows atoms near the glass or the floor any more, and the heat exchange with the
  plate, glass and air works relative to the material's own flow, so it no longer acts as a
  drag that slowed falling lumps, most of all just above the hot plate. Only atoms that touch
  the glass or floor feel friction. A lump dropped from the top of the jar now lands at full
  speed, flashes hot for a moment (iron about 1200 K) and cools on the plate within a second.
- **Metal melting points come from the pair well alone.** Without the springs, the pair well was
  recalibrated so that iron is solid at 1600 K and liquid at 2000 K, and lead melts between 500
  and 650 K. The many-body term only acts on atoms being pulled out of the metal, so it cannot
  favour the crystal. A version that also acted inside the metal froze mercury at room
  temperature.
- **The time step is shorter** (1 ms instead of 1.5 ms). With real cohesion a metal lattice is
  stiff, and the longer step leaked energy into it: an iron lump heated by itself at about
  30 K per second. Per-substep rates are rescaled so that everything else happens at the same
  rate per simulated second.
- **The hand works at the scale of the simulation.** On screen, the hand moves things at a few
  times the thermal speed of their atoms, about 10^13 times harsher than a real hand. Pulling a
  lump along by the patch under the cursor therefore stretched it plastically and heated it,
  which real handling never does to something a few nanometres across. The hand now carries a
  solid as a whole: its weight and the hand's push are shared over the atoms touching each
  other, the way gravity acts. The push depends only on the motion of the body as a whole, so
  it neither heats, cools nor freezes the atoms, and the heat exchange with the jar no longer
  drags the lump back. Once the lump is let go, or when it hits the glass, nothing but its own
  bonds holds it together. A liquid is still dragged by the patch alone, so it stretches and
  drips.

Still open:

- **Iodine:** its model is far too volatile. At 293 K its solid has regions above its 387 K
  melting point, so it crumbles.
- **Metals melt gradually** (as 2D solids do), and gallium and mercury freeze some way below
  their real melting points (gallium is solid at 240 K and liquid at 340 K).
- **Impacts are violent at this scale:** the hand's top speed is near the speed of a rifle
  bullet in the real-world jar, so a lump thrown at the glass warms up noticeably.

---

## Summary of the built-in experiments

| Preset | Real world | What the sandbox does (after the fix pass) | Verdict |
| --- | --- | --- | --- |
| Hydrogen + oxygen (spark) | Explodes to water | Burns to water. Flame peaks around 4900 K (real ~3000 K) **(README)** | Good |
| Sodium + chlorine | Vigorous reaction, white NaCl smoke | Correct product; ~25% of the sodium reacted after 20 s | Partial |
| Drop a blob of water | Falls, pools, wets the glass | Falls, pools, wets the glass; now gets the full 1.4 mol | Good |
| Steam on cold glass (450 K plate) | Water boils off and condenses on the walls | Still nearly all liquid after 30 s **(README)** | Poor |
| Methane + oxygen (spark) | Burns completely to CO₂ + 2 H₂O | Burns to water, but carbon ends mostly in CH₂O₃ and polymers, little CO₂ **(README)** | Poor |
| Iron in hot oxygen | Black oxide scale grows on the iron | A solid oxide scale (between Fe₃O₄ and Fe₂O₃) | Good |
| Carbon + oxygen | Inert cold; burns when hot | Inert at room temperature; burns at 1500 K to surface oxides and some CO₂ | Partial |
| Noble gases | Inert; argon liquefies at 87 K | Inert; a 70 K plate only cools the jar to ~100 K | Partial |
| Liquid nitrogen at 60 K | N₂ condenses within seconds | Takes ~25 s; condensed at a reading above N₂'s critical point | Partial |
| Freeze water into ice | Freezes to hexagonal ice | Freezes correctly at 200 K | Good |
| Melt and recast iron | Melts at 1811 K, recrystallises when cooled | Melts, glows, and recrystallises correctly | Good |
| Gallium melts near room temp | Melts at 303 K, freezes below it | Melts and freezes correctly | Good |

Things that match reality well: H₂/O₂ and H₂/Cl₂ ignition by spark; H₂ + F₂ exploding without
one; N₂ + H₂ not reacting without a catalyst; H₂/O₂ not self-igniting at 650 K; methane
self-igniting at an 850 K plate but not at 700 K (real autoignition ~810–850 K); gold untouched
by hot oxygen; iron untouched by chlorine until heated; iron reacting with hot steam to release
hydrogen; mercury liquid at room temperature and freezing below 234 K; iron, gallium and ice
phase changes; evaporative cooling of water under vacuum; dry ice subliming.

---

## Reactions

### 1. Carbon monoxide polymerises on its own at room temperature — High
**Status: Fixed.** CO's carbon is now a closed-shell lone pair, not a radical. CO stays CO
indefinitely, and radicals still add to it (O + CO → CO₂).
- **Setup:** 0.4 mol CO, 293 K, no oxygen.
- **Sandbox (before):** within 15 s, 98% of the CO turned into C₂O₂, C₃O₃, C₄O₄, C₅O₅ and an
  (CO)ₙ network, releasing heat.
- **Real:** CO is a stable gas indefinitely. C₂O₂ (ethylenedione) is not a stable molecule.

### 2. Hydrocarbon and alcohol flames never reach CO₂ — High (README)
**Status: Open.** A sparked methane/oxygen mix burns and makes plenty of water, but the carbon
still ends mostly in CH₂O₃ and in carbon–oxygen polymers, with only ~0.01 mol CO₂. The
high-temperature route through formaldehyde and CO is not modelled.
- The preset text says "burn methane into CO₂ and H₂O", which still mostly doesn't happen.

### 3. Ethanol decomposes at room temperature — High
**Status: Improved.** The cause was the overlap flash on adding molecules (see Root causes).
Ethanol at 293 K is now ~95% intact after 15 s; a slow condensation into ether-like
(C₄H₁₀O)ₙ and water remains.
- **Real:** ethanol is stable at room temperature.

### 4. Gold oxidises in oxygen — High
**Status: Fixed.** Au–O now has the weak bond of gold's unstable oxide (Au₂O₃ gives off its
oxygen above ~430 K). Gold on a 900 K plate in oxygen stays untouched.

### 5. Iron and chlorine react violently at room temperature — High
**Status: Fixed.**
- Metals that don't react by electron transfer must break through their own halide layer
  (a passivation barrier).
- A gas reacting on a metal surface now goes with the surface's temperature.

Iron now sits unreacted in chlorine at room temperature and burns on a 700 K plate. FeCl₃
holds together as a solid up to its 579 K sublimation point.

### 6. Sodium and potassium barely react with water — High (README)
**Status: Open.** Sodium still reacts slowly, mostly at its surface, and makes NaH and Na₂O
alongside NaOH. The likely root cause is that ions in solution aren't modelled: the NaOH
crust cannot dissolve away to expose fresh metal.
- **Real:** both react violently, and potassium ignites the hydrogen.

### 7. Ammonia + hydrogen chloride makes radicals instead of NH₄Cl — High
**Status: Improved.** The radicals came from the overlap flash; NH₃ and HCl now simply
coexist. No NH₄Cl forms yet, because dative bonds and four-bonded nitrogen aren't modelled
**(README)**.
- **Real:** an instant white smoke of solid NH₄Cl.

### 8. Acid–base neutralisation produces radicals — Medium
**Status: Improved.** No more radicals, but HCl + NaOH in water now doesn't react at all,
since ions in solution aren't modelled **(README)**.

### 9. Hydrogen + fluorine does not react in the dark — High
**Status: Fixed.** A molecule pulling one atom off another (H₂ + F₂ → HF + H + F) now uses the
abstraction barrier. That starts the dark explosion of H₂ + F₂, while the same step is far
uphill for chlorine and oxygen, so H₂/Cl₂ and H₂/O₂ still need a spark.

### 10. Hydrogen peroxide forms impossible higher peroxides — Medium
**Status: Fixed.**
- Oxygen no longer forms O–O–O chains on molecules that contain other elements, and pure
  oxygen stops at ozone.
- Hot peroxide (750 K plate) decomposes quickly to water and O₂.
- At 450 K it mostly evaporates and decomposes slowly. That's reasonable for the pure vapour;
  the fast decomposition seen in practice is catalysed by surfaces and metal traces.

### 11. Sulfur does not burn — Medium (README)
**Status: Fixed.** Sulfur on a 700 K plate burns to SO₂ and SO₃ (plus some sulfur-oxide
polymer).

### 12. Metal oxides and chlorides appear as small gas molecules — Medium
**Status: Improved.** Iron heated in oxygen now grows a solid oxide scale: large lattices
between Fe₃O₄ and Fe₂O₃ in composition. Exact stoichiometry is not enforced.

### 13. Carbon may burn spontaneously at room temperature — Medium
**Status: Fixed.** The *Carbon + oxygen* preset is inert at 293 K in every seed tested.
- Free valences on a solid (graphite edges) meet the O₂ spin barrier plus a chemisorption
  barrier, giving graphite oxidation its real ~170 kJ/mol.
- Freshly added material doesn't react with other substances while it settles in.
- Stray carbon atoms knocked loose far below carbon's boiling point can only rejoin their lump.

### 14. Nitrogen and oxygen do not form NO when hot — Low
**Status: Open.**
- **Real:** NO forms above ~2000 K (Zeldovich mechanism), about 1% at equilibrium.

---

## Phases and thermodynamics

### 15. Covalent network solids never warm up — High
**Status: Fixed for graphite.** Temperature is now kinetic energy over the degrees of freedom
actually left (2 per atom minus one per rigid bond, with a measured correction for rings and
stiff angles in networks).
- Graphite on an 800 K plate reads close to the plate with its chemistry switched off (it was
  ~250 K). With chemistry on it runs warm, about 1000 K, because the bonds at its edges keep
  breaking and re-forming and that leaks a little heat into it.
- Silicon (four bonds per atom, almost fully rigid in 2D) still reads low: ~190 K on a 293 K
  plate (it was 20–120 K).

### 16. Pressure gauge reads vapour that isn't there — Medium
**Status: Fixed.** The gauge and the vacuum pump no longer count molecules touching a liquid
or solid surface.

| Sample (293 K) | Gauge before | Gauge after | Real |
| --- | --- | --- | --- |
| Gallium (280–340 K) | 0.05–0.07 atm | 0.001 atm | ~0 |
| Mercury | 0.10 atm | 0.013 atm | 0.0000017 atm |
| Iodine | 0.26 atm | 0.019 atm | 0.0003 atm |
| Iron at 2100 K | 0.2–0.3 atm | 0.06 atm | ~0.0001 atm |

What remains for mercury, iodine and iron is real free vapour in the model: those substances
are too volatile (see the README on metal cohesion).

### 17. Liquids don't boil above their boiling point — High (README)
**Status: Open.** The *Steam on cold glass* preset is still nearly all liquid after 30 s on a
450 K plate. Sulfuric acid and bromine were not re-checked.

### 18. The jar runs far below the plate — Medium
**Status: Open.**
- **Hot end:** a 3000 K plate holds an N₂/O₂ jar at ~2050 K.
- **Cold end:** 60–70 K plates hold the jar at 100–140 K, so argon never liquefies on the 70 K
  plate (the noble-gas preset asks you to try this).

### 19. Condensed matter at impossible temperatures — Medium
**Status: Open.**
- **Liquid nitrogen preset:** N₂ is reported as condensed while the jar reads 131–141 K,
  above nitrogen's 126 K critical temperature.
- **Water:** stays mostly liquid with the jar at 250–260 K for 30 s.

### 20. CO₂ polymerises when cold — High
**Status: Fixed.** The polymer came from the overlap flash on adding the molecules. Dry ice now
stays CO₂ at 150 K and sublimes cleanly when warmed.

### 21. Sulfur makes the wrong molecules — Low
**Status: Open.** Sulfur is still built from chains and fragments (S₄–S₇) rather than S₈ rings.

---

## Accounting and display

### 22. A network with a few stray atoms vanishes from Contents — Medium
**Status: Fixed.** When a network's reduced formula is no real formula unit (graphite with a
few oxygens reduces to C₁₅₀O), Contents counts moles of its main element, so the lump keeps its
~0.3 mol.

### 23. Presets silently add less than they ask for — Low
**Status: Fixed.** `addMolecules` grows its box until the requested amount fits (a caller can
pass `fixed: true` to keep a shaped box, like a film against the glass).

### 24. No glow from hot materials — Low
**Status: Fixed.** In the real-world view, liquids and solids glow from dull red around 800 K
(the Draper point) to orange and yellow-white when molten.

---

## Reference data

### 25. Small inconsistencies in the data — Low
**Status: Fixed**, except NaH's bond class.
- **Gallium boiling point:** the README now matches `elements.js` (2673 K).
- **Ethanol in `phase.js`:** the triple point is now 159 K, matching the melting point.
- **Carbon:** `elements.js` now gives the sublimation point (~3915 K) for both melting and
  boiling.
- **Helium:** `elements.js` no longer gives a melting point; helium does not solidify at 1 atm.
- **Bond classes:** AlCl₃ and FeCl₃ now hold together as solids up to their real
  sublimation points (453 K and 579 K), whatever their bond class. NaH and the other s-block
  hydrides are still covalent rather than ionic.
