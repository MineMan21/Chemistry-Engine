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

---

## Summary of the built-in experiments

| Preset | Real world | What the sandbox does | Verdict |
| --- | --- | --- | --- |
| Hydrogen + oxygen (spark) | Explodes to water | Burns to water. Flame peaks around 4900 K (real ~3000 K) **(README)** | Good |
| Sodium + chlorine | Vigorous reaction, white NaCl smoke | Correct product, but only ~20% of the sodium has reacted after 20 s. Some NaCl is listed as *gas* | Partial |
| Drop a blob of water | Falls, pools, wets the glass | Falls, pools, wets the glass. Gets 0.44 mol instead of the 1.4 mol requested (see #20) | Good |
| Steam on cold glass (450 K plate) | Water boils off and condenses on the walls | 92% stays liquid after 30 s **(README)** | Poor |
| Methane + oxygen (spark) | Burns completely to CO₂ + 2 H₂O | No CO₂ at all. Carbon ends in CH₄O₄, CH₂O₃ and other peroxide-like species; half the methane is left in the UI run **(README)** | Poor |
| Iron in hot oxygen | Black oxide scale (FeO/Fe₃O₄/Fe₂O₃) grows on the iron | Oxides form, but many are small *gas-phase* Fe–O molecules with odd formulas (Fe₂O, Fe₃O₂) | Partial |
| Carbon + oxygen | Nothing at room temperature; glows and burns to CO/CO₂ when hot | Graphite never heats up, so it can't burn. Occasionally burns spontaneously at room temperature instead (1 of 3 seeds) | Poor |
| Noble gases | Inert; argon liquefies at 87 K | Inert (correct). A 70 K plate only cools the jar to ~100 K, so argon never liquefies | Partial |
| Liquid nitrogen at 60 K | N₂ condenses within seconds | Takes ~25 s; N₂ reported as condensed while the jar reads 130–140 K, above N₂'s 126 K critical point | Partial |
| Freeze water into ice | Freezes to hexagonal ice | Freezes correctly at 200 K | Good |
| Melt and recast iron | Melts at 1811 K, recrystallises when cooled | Melts and recrystallises correctly | Good |
| Gallium melts near room temp | Melts at 303 K, freezes below it | Melts and freezes correctly | Good |

Things that match reality well: H₂/O₂ and H₂/Cl₂ ignition by spark; N₂ + H₂ not reacting
without a catalyst (with or without iron at 700 K); H₂/O₂ not self-igniting at 650 K; methane
self-igniting at an 850 K plate but not at 700 K (real autoignition ~810–850 K); mercury
liquid at room temperature and freezing below 234 K; iron, gallium and ice phase changes;
evaporative cooling of water under vacuum; CO₂ stable at room temperature.

---

## Reactions

### 1. Carbon monoxide polymerises on its own at room temperature — High
- **Setup:** 0.4 mol CO, 293 K, no oxygen.
- **Sandbox:** within 15 s, 98% of the CO turns into C₂O₂, C₃O₃, C₄O₄, C₅O₅ and an (CO)ₙ
  network, releasing heat (the jar briefly reaches 390 K).
- **Real:** CO is a stable gas indefinitely. C₂O₂ (ethylenedione) is not a stable molecule.
- **Likely cause:** CO's carbon has only 3 of its 4 valences used, so the engine treats it as a
  radical that bonds to other CO molecules with almost no barrier.
- **Knock-on effect:** CO + O₂ forms these chains before you even spark it. After the spark the
  main products are C₂O₂, C₂O₄, CO₃ and C₃O₅, with only a little CO₂.

### 2. Hydrocarbon and alcohol flames never reach CO₂ — High (README)
- **Methane + O₂, sparked:** water forms, but carbon stops at CH₂O₃, CH₂O₄, CH₄O₃, CH₄O₄ and
  C₂H₂O₄. Zero CO₂, in both the headless and browser runs.
- **Methane + O₂ on an 850 K plate:** same outcome.
- **Benzene + O₂, sparked:** mostly polymer networks ((C₃O)ₙ, (C₃H)ₙ) plus 0.005 mol CO₂.
- **Ethanol + O₂, sparked:** H₂O, H₂O₂ and CH₃O₃; no CO₂.
- The preset text says "burn methane into CO₂ and H₂O", which never happens.

### 3. Ethanol decomposes at room temperature — High
- **Setup:** pure ethanol, 293 K.
- **Sandbox:** within 5 s about 35% of it turns into water, H₂, ethane and polymer networks
  such as (C₃H₇O)ₙ and (C₄H₁₀O)ₙ.
- **Real:** ethanol is stable at room temperature.

### 4. Gold oxidises in oxygen — High
- **Setup:** Au + O₂ on a 900 K plate.
- **Sandbox:** gold forms AuO, Au₂O₂, Au₂O₃ and Au₃O₄ as gas molecules. The reaction is
  exothermic enough to push the jar to ~1900 K.
- **Real:** gold does not react with O₂ at any temperature. Au₂O₃ is endothermic and
  decomposes above about 430 K.
- **Likely cause:** with an Au–O electronegativity gap of 0.9, the Au–O bond is classed as
  covalent and gets a Pauling-rule bond energy, which makes oxidation favourable.

### 5. Iron and chlorine react violently at room temperature — High
- **Setup:** Fe + Cl₂ at 293 K.
- **Sandbox:** the jar heats itself to ~1300 K within 10 s.
- **Real:** dry chlorine barely attacks iron at room temperature; chlorine is shipped in steel
  cylinders. Iron burns in chlorine only once heated to roughly 500 K.
- **Products:** FeCl₃, FeCl₂ and FeCl are all listed as *gases*. Real FeCl₃ is a solid that
  sublimes near 590 K.

### 6. Sodium and potassium barely react with water — High (README)
- **Setup:** 0.2 mol of the metal dropped onto 0.56 mol of water.
- **Sodium:** after 15 s only ~3% has reacted.
- **Potassium:** essentially nothing happens.
- **Products:** where sodium does react, it makes Na₂O and NaH. The real products are NaOH +
  H₂.
- **Real:** both react violently, and potassium ignites the hydrogen. This is the textbook
  alkali-metal demonstration.

### 7. Ammonia + hydrogen chloride makes radicals instead of NH₄Cl — High
- **Setup:** NH₃ + HCl gases mixed at 293 K.
- **Sandbox:** no ammonium chloride. Instead H₂, Cl atoms, Cl₂, NH₂ and NH form at room
  temperature.
- **Real:** an instant white smoke of solid NH₄Cl.
- Ammonium can't form because dative bonds and N with 4 bonds aren't modelled **(README)**.
  Producing free radicals and H₂ from these two stable molecules is an extra problem.

### 8. Acid–base neutralisation produces radicals — Medium
- **Setup:** HCl + NaOH in water.
- **Sandbox:** some NaCl + H₂O forms, along with H₂, Cl atoms, H atoms, Cl₂ and Na₂O.
- **Real:** HCl + NaOH → NaCl(aq) + H₂O cleanly. No gas is evolved and no radicals form.
- Ions in solution aren't modelled **(README)**, but the radical side products are
  unrealistic either way.

### 9. Hydrogen + fluorine does not react in the dark — High
- **Setup:** H₂ + F₂ at 293 K, no spark.
- **Sandbox:** no reaction in 15 s.
- **Real:** reacts explosively on contact, even in the dark and at very low temperatures. This
  is the classic contrast with H₂ + Cl₂, which needs light or a spark.

### 10. Hydrogen peroxide forms impossible higher peroxides — Medium
- **Setup:** H₂O₂ on a 450 K plate.
- **Sandbox:** partial decomposition to H₂O + O₂, which is right. But it also builds H₂O₃,
  H₂O₄, H₂O₅ and HO₃ in quantity.
- **Real:** H₂O₂ → H₂O + ½ O₂. Polyoxides beyond H₂O₃ aren't observable.

### 11. Sulfur does not burn — Medium (README)
- **Setup:** sulfur + O₂ on a 700 K plate.
- **Sandbox:** no SO₂ in 30 s.
- **Real:** sulfur ignites in air at about 505 K with a blue flame, giving SO₂.

### 12. Metal oxides and chlorides appear as small gas molecules — Medium
- **Iron in O₂:** Fe₂O, Fe₃O₂, Fe₂O₂, Fe₄O₄ and FeO, many listed as gas.
- **Copper in O₂:** Cu₂O₂, Cu₃O₂ (gas).
- **Gold in O₂:** AuO (gas).
- **Magnesium in CO₂:** Mg₂O₂, CMgO₂, CMg₂.
- **Iron in Cl₂:** FeCl (gas).
- **Real:** these are non-volatile ionic solids with fixed stoichiometry (FeO/Fe₃O₄/Fe₂O₃, CuO,
  Cu₂O, MgO). They grow as a scale on the metal and do not float around as small molecules.
- In the real-world view, these show up as rust-coloured blobs hovering mid-jar.

### 13. Carbon may burn spontaneously at room temperature — Medium
- **Setup:** the *Carbon + oxygen* preset on a 293 K plate, run with 3 random seeds.
- **Sandbox:** one seed ran away. Within 5 s there were tens of thousands of reactions per
  second, the jar reached ~2100 K, and CO₃, CO₄ and CO₅ formed. The other two seeds were inert,
  which is correct.
- **Real:** graphite in O₂ is inert at room temperature; it ignites at about 700–1000 K.

### 14. Nitrogen and oxygen do not form NO when hot — Low
- **Setup:** N₂ + O₂ on a 3000 K plate. The jar only reaches ~2050 K (see #18).
- **Sandbox:** no NO in 30 s.
- **Real:** NO forms at these temperatures (Zeldovich mechanism, the reason engines and
  lightning make NOₓ). The equilibrium is about 1% NO near 2000–2500 K.

---

## Phases and thermodynamics

### 15. Covalent network solids never warm up — High
- **Setup:** graphite or silicon in an insulated jar on a 293 K plate, then on an 800 K plate.
- **Graphite:** the thermometer reads 80–115 K on the 293 K plate and ~250 K on the 800 K plate.
- **Silicon:** stays at 20–120 K throughout.
- **Comparison:** iron and salt under the same conditions track the plate within ~30 K.
- **Real:** a lump of graphite on a hot plate reaches the plate temperature.
- **Consequences:** this is probably why carbon and sulfur "hardly burn". It also distorts the
  jar thermometer whenever a network solid is present.

### 16. Pressure gauge reads vapour that isn't there — Medium
The gauge counts surface atoms of liquids and solids as gas, so non-volatile materials show a
substantial pressure:

| Sample (293 K unless noted) | Gauge | Real vapour pressure |
| --- | --- | --- |
| Gallium, 0.5 mol at ~300–340 K (no free atoms seen in the jar) | 0.05–0.07 atm | ~10⁻³⁰ atm |
| Mercury, 0.5 mol | 0.10 atm | 0.0000017 atm |
| Iodine, 0.4 mol | 0.26 atm | 0.0003 atm |
| Water, 0.56 mol | 0.06–0.085 atm | 0.023 atm |
| Iron at 2100 K | 0.2–0.3 atm | ~0.0001 atm |

Counting molecules that actually sit apart from the rest gives about 5–10× less than the
gauge, but that is still far above reality for Hg, I₂ and Fe.

### 17. Liquids don't boil above their boiling point — High (README)
- **Water:** on a 400 K plate (jar ~390 K) it stays 100% liquid after 30 s. The *Steam on cold
  glass* preset (450 K plate) is still 92% liquid after 30 s, so there is little steam to
  condense on the glass.
- **Water after the methane flame:** the water formed is drawn as liquid drops hanging in the
  middle of a 450 K jar. At that temperature 0.1 mol of water should be all vapour.
- **Sulfuric acid:** the opposite error. It is listed as *gas* at 293 K, although it is
  essentially non-volatile (boils at 610 K).
- **Bromine:** 0.4 mol at 293 K is listed as all gas. At the real vapour pressure of 0.23 atm,
  about 40% should be liquid.

### 18. The jar runs far below the plate — Medium
- **Hot end:** a 3000 K plate holds an N₂/O₂ jar at ~2050 K.
- **Cold end:** 60–70 K plates hold the jar at 100–140 K, so argon never liquefies on the 70 K
  plate (the noble-gas preset asks you to try this).
- Some heat loss through the glass is real, but the 900 K gap makes the plate setting a poor
  guide to what is happening.

### 19. Condensed matter at impossible temperatures — Medium
- **Liquid nitrogen preset:** N₂ is reported as condensed while the jar reads 131–141 K. That
  is above nitrogen's critical temperature (126 K), where no liquid can exist. Even at 95–98 K
  (with the jar insulated) a real liquid would boil at ~5 atm, while the gauge reads 0.2 atm.
- **Water:** stays mostly liquid with the jar at 250–260 K for 30 s (0–14% crystallised).
  Supercooling is real, but bulk water in a glass jar normally freezes well before this.
- **Salt:** with 0.3 mol Na + 0.15 mol Cl₂ heated to 1000–1400 K, unreacted sodium and
  chlorine coexist for over a minute. The jar temperature *drops* (958 K → 458 K) as the plate
  is turned up.

### 20. CO₂ polymerises when cold — High
- **Setup:** CO₂ on a 150 K plate, then warmed back to 293 K.
- **Sandbox:** half of it becomes a covalent (CO₂)ₙ network plus C₂O₄ and C₃O₆. The network is
  still there after warming back to 293 K.
- **Real:** dry ice is a molecular solid that sublimes at 195 K. Polymeric CO₂ only exists
  above ~40 GPa.
- The same networks appear in the magnesium + CO₂ experiment.

### 21. Sulfur makes the wrong molecules — Low
- **Sandbox:** solid sulfur at room temperature contains S₄, S₅, S₆ and S₇ fragments listed as
  *gas*, plus an (S)ₙ network.
- **Real:** sulfur is S₈ crown rings and its vapour pressure at 293 K is negligible.

---

## Accounting and display

### 22. A network with a few stray atoms vanishes from Contents — Medium
- **What happens:** when a single O atom bonds to a 0.3 mol graphite lump, the lump is listed
  as `(C₁₅₀O)ₙ 0.002 mol` (or `(C₂₉₆O₇)ₙ 0.001 mol`) instead of ~0.3 mol.
- **Cause:** `formulaSize` divides the atom count by the reduced formula. With no common
  divisor, the whole lump counts as one formula unit.

### 23. Presets silently add less than they ask for — Low
- **Cause:** `addMolecules` stops at the number of lattice sites in the spawn box.
- **Examples:**
  - *Drop a blob of water* requests 1400 molecules and gets 440 (0.44 mol).
  - *Freeze water* requests 900 and gets 476.
  - *Steam* requests 1200 and gets 924.
  - The NaCl molecule chip placed 56 of 300.
- No message tells the user.

### 24. No glow from hot materials — Low
- Iron at 1100–2100 K, molten salt and flames are drawn in their room-temperature colours.
  Only the hot plate glows.
- Real iron is visibly red above ~800 K and orange-white when molten.

---

## Reference data

### 25. Small inconsistencies in the data — Low
- **Gallium boiling point:** the README says 2477 K; `elements.js` has 2673 K.
- **Ethanol in `phase.js`:** the triple point is 150 K but the melting point is 159 K. With
  the default melting-line slope, the diagram puts ethanol's melting at 1 atm at ~150 K, not
  159 K.
- **Carbon:** `elements.js` gives a melting point of 3823 K and a boiling point of 4098 K.
  Graphite actually sublimes at ~3900 K at 1 atm without melting.
- **Helium:** `elements.js` gives a 1 K melting point. Helium does not solidify at 1 atm, as
  `phase.js` correctly notes.
- **Bond classes:**
  - AlCl₃ is classed as ionic (Al–Cl gap 1.55) although it is a molecular/covalent solid that
    sublimes at 453 K.
  - FeCl₃ is classed as covalent molecules, so it behaves as a gas.
  - NaH and other s-block hydrides fall below the 1.5 gap and are covalent rather than ionic.
