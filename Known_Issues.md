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
`tests/sim.test.js`. Any new problem found while working on the sandbox is added to the table
below as soon as it is found, even if it is not fixed.

---

## Status after the fix pass

| # | Issue | Status |
| --- | --- | --- |
| P | Reactive substances (sulfuric acid, benzene) react with themselves into polymers; sulfuric acid still did once something was dropped in | **Fixed**: see below and "Sulfuric acid with something dropped in" |
| L | Flames give off no light | **Fixed**: see below |
| W | 2:1 H₂/O₂ burns to junk; sparking or heating water makes new molecules | **Fixed** except a few lingering radicals: see below |
| D | Dropped iron (and tungsten) dents; metals too malleable | **Fixed**: directional bonding, and much gentler gravity for falling solids; see below |
| B | Solids float on liquids (iron on water) | **Fixed** for iron, gold and lead; liquid mercury now sinks partway: see below |
| E | Elements dropped on their own misbehave (silicon, sulfur, tellurium flash to 10⁴ K; silicon hangs in mid-air; sodium splats) | **Fixed** for 86 of the 96 solid elements: see "Every element dropped" below |
| E1 | Rubidium, caesium and francium (1 mol) are as wide as the jar, drift partly as vapour and do not land within 60 s | Open: see "Every element dropped" |
| E2 | Potassium and gallium, just below their melting points, splash on landing (22–30% of contacts rearranged) | Open |
| E3 | Iodine, astatine, tennessine and oganesson fall apart: crystals held only by dispersion forces are too weak in 2D | Open |
| E4 | Silicon and germanium read ~100 K at room temperature; sulfur, selenium and tellurium read 500–570 K after landing (temperature of nearly rigid networks and chains) | Open |
| E5 | Liquid mercury only partly sinks in water (its settling force is capped at ~1/4 of water's); solids neither float nor sink properly in a liquid metal | Open |
| E6 | Falling liquids and solids still warm on landing (water ~150–220 K, solids 50–150 K for a moment), because settling gravity is ~10¹¹–10¹³ times real | Open: inherent to the exaggerated gravity |
| E7 | Molten gallium at 340 K (37 K above its melting point) keeps 10–15% of its atoms in ordered patches; whether a run reads below 10% depends on the random seed | Open |
| G | Solids break when grabbed or thrown; moving them heats them | **Fixed** except iodine and sulfur: see the table below |
| 1 | CO polymerises at room temperature | **Fixed** |
| 2 | Hydrocarbon flames never reach CO₂ | Open |
| 3 | Ethanol decomposes at room temperature | **Improved**: ~95% intact after 15 s |
| 4 | Gold oxidises | **Fixed** |
| 5 | Iron + chlorine react at room temperature | **Fixed** |
| 6 | Sodium and potassium barely react with water | **Fixed**: they (and lithium, calcium, barium) now give their hydroxide and H₂ (2 Na + 2 H₂O → 2 NaOH + H₂), no hydride; see issue 6 below. Barium still leaves a little BaH₂ (issue 45) |
| 7 | NH₃ + HCl makes radicals | **Improved**: no radicals; still no NH₄Cl |
| 8 | Neutralisation makes radicals | **Improved**: no radicals; still no reaction |
| 9 | H₂ + F₂ does not react in the dark | **Fixed** |
| 10 | Peroxide builds H₂O₃–H₂O₅ | **Fixed** |
| 11 | Sulfur does not burn | **Fixed** |
| 12 | Metal oxides as small gas molecules | **Improved**: rust is a solid scale |
| 13 | Carbon burns at room temperature | **Fixed** |
| 14 | No NO from hot N₂ + O₂ | Open |
| 15 | Network solids never warm up | **Fixed**: graphite warms (it reads ~1.5x hot, issue 29); silicon reads low |
| 16 | Gauge reads vapour that isn't there | **Fixed**; Hg and I₂ remain too volatile |
| 17 | Liquids don't boil above their boiling point | **Fixed (partly)**: a liquid hotter than its boiling point at the jar's pressure (from the real vapour curve) loses cohesion and boils off, taking its latent heat. On the test's 520 K plate water is 25–28% liquid (was 37–39%); an insulated pool on a 400 K plate barely changed (~54% liquid after 30 s, was 55–63%; see below). See "Boiling, condensing and the plate" |
| 18 | The jar runs far below the plate | **Explained / fixed at the cold end**: an insulated jar reaches the plate (3000 K plate: ~3000–3200 K); with the glass losing heat to the 293 K room the gas settles between the two (~2050 K), as a real jar would. On a 70 K plate argon now condenses (58% after 30 s, was 17%): its atoms were bound too weakly |
| 19 | Condensed matter at impossible temperatures | **Fixed / explained**: no liquid can last above its critical temperature or boiling point now. The liquid nitrogen's own temperature is ~77 K; the 131–141 K reading is the jar's average, mostly the warmer gas over it. Water on a 250 K plate does freeze (62% ice at 15 s, over 80% from 30 s); the old reading counted ice as condensed liquid |
| 20 | CO₂ polymerises when cold | **Fixed** |
| 21 | Sulfur makes the wrong molecules | Open |
| 22 | A network with a stray atom vanishes from Contents | **Fixed** |
| 23 | Presets add less than they ask for | **Fixed** |
| 24 | No glow from hot materials | **Fixed** |
| 25 | Data inconsistencies | **Fixed** (except NaH's bond class) |
| 26 | Iron in 0.6 mol of chlorine ignites by itself at room temperature after 6–12 s in a dry jar (to ~1200–1600 K, FeCl/FeCl₂/FeCl₃); issue 5's fix holds for small amounts only. With water in the jar it stays at a couple of Fe–Cl bonds | Open (found while checking the gas rendering) |
| 27 | About a third of 0.6 mol of chlorine counts as condensed at room temperature (clinging to the iron and the glass), though Cl₂ boils at 239 K | Open |
| 28 | Metal hardness follows the melting point only: at the same impact speed tungsten (real hardness 5.6x iron's) dents about as much as iron. Thrown by hand, heavy metals now leave the hand slower (equal energy per atom), so tungsten lands whole where iron dents | **Improved**: see "Throwing tungsten" below |
| 29 | Graphite on an 800 K hot plate (insulated jar) reads 900–1300 K depending on the run, mostly ~1200 K: about 1.5x hot. Its test passed only because its one random seed happened to read 936 K; it now averages three runs | Open |
| 30 | The Moles view lags more than the other two | **Fixed**: it rebuilt and uploaded a full-screen image every frame (about 5 million pixels on a high-resolution screen); it now draws its squares directly, grouped by colour (0.7 ms a frame, was 3.1 ms at 2x pixel density, plus the upload) |
| 31 | Sulfuric acid does not take part in acid–base or acid–metal reactions: with NH₃, NaOH or NaCl nothing happens (real: ammonium sulfate; sodium sulfate and water; HCl and NaHSO₄); copper stays inert, sodium gives Na₂O and NaOH instead of Na₂SO₄ + H₂. These go through protons and ions, which the model does not have (as issues 7 and 8) | **Improved**: metals above hydrogen now dissolve in dilute acid (issue 37); the rest open |
| 32 | A lone fast molecule in a nearly empty patch of vapour read thousands of kelvin, enough to break its bonds and set its neighbours off | **Fixed**: sparse cells lean harder on the jar's temperature (T_PRIOR 2 → 8) |
| 33 | A deep pool of liquid warms as it settles: 6 mol of water dropped as one full-width blob goes from 293 K to ~333 K in 4 s; filled as a mixture it starts at its settled density and warms less (water to 85% of the jar ~300 K, 95:5 brine to 85% ~312 K after 3 s), and its lower layers end up ~25% denser than the top, so the liquid is too compressible under the exaggerated settling gravity (a bed filling 30% of the jar stays within ~10 K) | Open (found while adding mixtures) |
| 34 | Liquids and solids take up ~130x their real volume next to the gas scale: the 2D jar full of water holds about 10 mol (a real 24 L jar: 1300 mol), while 1 mol of gas reads 1 atm as in 24 L | Open: inherent to the 2D slice and the 1 mmol particle |
| 35 | Covalent network elements (B, C, Si, P, S, Ge, As, Se, Sb, Te) cannot be part of a mixture: they are only laid down as bonded sheets, and as loose atoms they would bond in a frenzy | Open: they are left out of the mixture's element list |
| 36 | Iron dropped onto sulfuric acid and water (1:3) floats on it: after 16 s its centre is still above the liquid's surface. Packed flat, iron comes out only about twice as dense as water (really 7.9 times) and the acid mixture 1.15 times, so the small difference barely pushed the iron through the thick, viscous liquid | **Fixed**: an immersed metal is pulled in proportion to its real density (new table of metal densities): iron goes under and settles on the floor; it sinks no faster than SINK_VMAX (issue 43) |
| 37 | Iron does not react with dilute sulfuric acid at room temperature (real: dissolves, giving H₂ and FeSO₄). The only route was a metal taking the acid's OH oxygen and setting the H free as an H atom, ~150 kJ/mol of activation | **Fixed**: metals with a negative standard potential dissolve in one step, two protons (one may come from a water molecule) leaving as H₂; needs water about the acid, so cold concentrated acid stays inert; copper does not react. See "Iron in dilute sulfuric acid" |
| 38 | The iron sulfate that forms does not dissolve away into the water (ions are not solvated), so it gathers around the iron and the reaction slows: a small lump in 1:3 acid loses about 15% in 8 s, then much more slowly | Open |
| 39 | Magnesium in dilute sulfuric acid: the reaction's heat (457 kJ/mol, the real value) speeds up magnesium's own radical reaction with water, giving Mg–H, MgO and H atoms as side products | Open |
| 40 | Material dropped inside a liquid (a click below its surface) has no free space where it lands, so its atoms are placed wherever they fit nearby: a 0.15 mol iron lump came out as a loose cloud 1400 units wide and tall instead of a lump | **Fixed**: a solid or liquid released inside material is released just above it instead, and falls in as a lump |
| 41 | Iron lying in cold concentrated sulfuric acid now and then breaks an acid's S–OH bond (HO₃S and OH radicals): 1–3 times in 24 s for 0.15 mol of iron. Seen since iron sinks into the acid instead of resting on it; none in 24 s before | Open |
| 42 | H₂ given off under the liquid stays among its molecules for a while (Contents lists it as liquid/solid) instead of bubbling straight out | **Explained, partly improved**: see "Gases in water and the hydrophobic effect". Gases barely dissolve from above (0–2 of 200 atoms), as in reality; gas made under water leaves only by diffusing out, slowly, because a few molecules cannot make a bubble (at this scale a real one could not exist either) |
| 43 | After issue 36's fix, a ball of iron dropped into water slammed in: pulled 3.6 times as hard as the water, it plunged at the bulk speed limit (900 units/s), flung the water up the glass and heated its leading face to 380–480 K, where water oxidised iron atoms off it; the ball could tear | **Fixed**: a metal in a liquid sinks no faster than 200 units/s (the liquid's drag; about how fast iron sank in water before): 0.4 mol of iron into water now peaks at ~310 K, with no iron oxidised and a smaller splash than before issue 36's change |
| 44 | A lump held by the hand and rammed into the glass at the hand's top speed heats to ~1900 K (iron, past its melting point; cell readings up to 2900 K) as it hits; lead to ~1150 K. Thrown, it bounces off and keeps most of its motion; held, all its motion (the hand gives up to ~1680 K worth per atom) turns to heat at once | **Fixed**: while the hand holds a solid that touches the glass (or for 4 s after throwing one), each part of it is kept to at most 80% of its melting point, the rest of the heat going into the hand and the glass; once against the glass the hand presses it on at no more than 30 units/s instead of ramming on. Iron now peaks around 900 K, lead 460 K, gold 470 K, and they still heat and dent on impact. See "Pushing a solid into the glass" |
| 45 | Barium in water still leaves a little BaH₂ and BaH (~0.006 mol from 0.06 mol of barium in 10 s): the hot metal takes up some of the H₂ it gives off faster than water hydrolyses the hydride | Open (found while fixing issue 6) |
| 46 | Hydrogen peroxide vapour decomposes far too slowly: 0.3 mol at 850 K in an insulated jar is still ~90% H₂O₂ after 12 s (real: gone in well under a second; at 600 K real vapour lasts hours, which the model gets right). Liquid peroxide on a hot plate decomposes quickly, but since issue 17's fix it boils off first (as real dilute peroxide would), so the hot-peroxide test now runs in an insulated jar | Open (found while fixing issue 17) |
| 47 | The real-world view built its picture (splatting, blurring and shading a fine grid) for the whole jar every frame, even when zoomed in or panned so that most of the jar was off screen | **Fixed**: it builds only the part on screen plus the blur's reach. See "Drawing only what is on screen" |
| 48 | Zoomed far into the atom view, atoms near the screen's edge popped in and out, and a bond between two atoms just off screen could vanish while crossing it: atoms were drawn only within a fixed 60 px of the edge, while a big atom is up to ~250 px across and a bond up to ~560 px long at full zoom | **Fixed**: the margin grows with the zoom |
| 49 | Between the mole and real-world views (scale 0.19–0.25) both are drawn and cross-faded, the slowest view: ~27 ms a frame with ~8,900 atoms in a 1400×900 window (Chromium, software rendering), most of it the mole view's thousands of half-transparent squares | Open (minor: a short range of zoom, passed through) |

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

## Substances reacting with themselves; light

Reported issues: very reactive substances, such as sulfuric acid, reacted with themselves into
huge polymer chains; and a hydrogen flame gave off no visible light.

A survey of every molecule in the palette, alone in the jar at room temperature for 15 s:

| Substance | Before | After |
| --- | --- | --- |
| Sulfuric acid | 99% condensed into S–O–S chains and water | 99% intact (91% at 600 K, where real sulfuric acid starts to break down) |
| Benzene | heated itself to ~9000 K and broke into CH₄ and carbon chains | 100% intact (also at 600 K) |
| Water, CO₂, CH₄, NH₃, CO, HCl, NaCl, NaOH, H₂O₂, ethanol | stable | stable (ethanol: 1% pairs up into an ether-like dimer) |

The elements (N₂, O₂, the halogens, alkali and alkaline-earth metals) stay as they are. Sulfur,
phosphorus, selenium, boron, silicon and carbon form networks, which their real solids are.
Mixtures (water with sulfuric acid, NaOH, ethanol, NH₃, HCl or CO₂; ethanol or H₂O₂ with
sulfuric acid) do not react either.

Causes and fixes:

- **Four-centre exchanges between whole molecules were too easy.** In the reaction AB + CD →
  AC + BD, two stable molecules swap partners in one step. Its barrier was 22% of the bonds it
  breaks, out of reach at room temperature but not at the brief local hot spots of a liquid's
  surface. Real reactions of this kind are symmetry forbidden; they go through radicals, ions or
  catalysts. Between molecules of six or more atoms the barrier is now 40%. Between two small
  molecules it stays at 22%, because there the exchange stands in for the radical chains of
  H₂ + Cl₂, of a methane flame or of decomposing H₂O₂, which this coarse model cannot sustain by
  itself.
- **Benzene had no aromatic stability.** Its "double bonds" added across like an alkene's,
  each releasing 90 kJ/mol, and the heat ran the reaction away. Benzene-like rings are now
  recognised, and breaking any bond of one costs the ring's 150 kJ/mol resonance energy.
- **Light.** Reactions that make a bond now give off a flash coloured by the bond (OH violet-blue,
  CH and C₂ blue, NO₂ yellow-green, metals their own colour), and hot metal atoms in a gas glow
  in their flame-test colours.
- **The spark lasts a moment of simulated time** rather than ten screen frames. On a slow
  computer it used to be too brief to light anything once it only heated. A freshly added gas of
  whole molecules is no longer held at room temperature while it "settles in" (meant for liquids,
  solids and assembling networks), which had smothered flames in the first two seconds, and the
  products of a reaction stop settling. The hydrogen/oxygen preset starts premixed.

Still open:

- **Benzene above ~900 K** breaks down to methane and carbon far too violently (real benzene
  pyrolyses slowly from about 1000 K).
- **Light is not physically calibrated:** brightness is set so a flame is clearly visible. A
  real hydrogen flame is much fainter in a lit room.

---

## Hydrogen flame, sparked water and heated water

Reported issue: a 2:1 H₂/O₂ mix burned to "a bunch of junk" rather than water, sparking water
split it into other molecules, and heating water made new molecules.

Causes and fixes:

- **Bonds formed from a distance snapped shut.** A reaction candidate is picked when two atoms
  touch, but a fast hot H atom can have moved on by the time the reaction fires. The new bond
  started at most 2.5× its rest length, so the constraint solver yanked the atoms across the gap
  in one step and flung them off with energy from nowhere. In an insulated jar, steam at about
  2000 K heated itself to 19,000 K. Reactions now need the pair to still be in contact. A new
  bond starts at the atoms' real distance and eases in.
- **Stretched-bond check broke bonds that were easing in.** It snapped any bond longer than 3×
  its rest length, including ones just made a little apart. Each snap charged the bond's full
  energy, and the atoms bonded again, around and around. The energy owed piled up to the whole
  heat of the flame, which kept the flame cool and then chilled the jar below room
  temperature. It now judges stretch against the length the bond is easing from.
- **The spark only heats.** It used to break bonds directly as well.

Results:

| Experiment | Before | After |
| --- | --- | --- |
| Sparked 2:1 H₂/O₂ | peaked ~11,000 K, still ~6000 K after 40 s; 81% water with H, OH, O, H₂O₂, HO₂ | peaks ~3000 K, back to ~300 K; ~91% water, rest mostly H₂ with traces of O₂, H₂O₂, OH |
| Water on a 3000 K plate, then cooled | steam reached 7200 K and stayed above 4000 K after cooling | steam reaches ~1900 K with ~3% split into H and OH; ~2% left as H and OH after cooling |
| Sparked water | ~10% split into H, OH, H₂O₂, HO₂ | ~6% split, ending mostly as H₂ with traces of O₂ and H₂O₂ (as flash-heated water really does) |

Real water stays water up to about 2000 K. It splits by about a tenth near 3000 K and much more
near 4000 K, into H₂, O₂, OH, H and O, and those recombine to water on cooling.

Still open:

- **A few radicals linger:** in the sparse simulated gas, stray H atoms and OH meet too rarely
  to finish recombining, so 1–2% stay unreacted at room temperature. Real ones would also
  recombine on the glass.

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
- **Water is somewhat too volatile:** 0.4 mol of water in the jar at 300 K keeps about a fifth of
  itself as vapour; real water would keep about a tenth (its vapour pressure is 0.035 atm).
- **Metals melt gradually** (as 2D solids do), and gallium and mercury freeze some way below
  their real melting points (gallium is solid at 200 K, only partly ordered at 240 K, and liquid
  at 340 K).
- **Impacts are violent at this scale:** the hand's top speed is near the speed of a rifle
  bullet in the real-world jar, so a lump thrown at the glass warms up noticeably. The same goes
  for a dropped lump: even the gentle settling force on falling solids is about 4×10¹¹ times
  real gravity, so a lump dropped from the top of the jar lands at a few hundred m/s in
  real-world terms and warms by 50–150 K for a moment (a real one warms by a hundredth of a
  kelvin). A drop of water falling from the top lands at about 640 units/s and warms by
  150–220 K for a moment (at the old liquid gravity of 650, by about 390 K). Landing gently enough for that to vanish would take minutes of simulated time per fall.

### Denting of dropped metal

Reported issue: iron dropped in the jar dented and was more malleable than iron.

What changed: metals now have directional bonding. Each atom prefers its neighbours at the
close-packed crystal's 60° angles, so the lattice resists shear and slip as well as compression.
The pair well was lowered so that iron, lead, mercury and gallium still melt where they should;
iron's melting is now sharper and closer to real (solid at 1650 K, liquid by about 2050 K).

| Iron lump (0.3 mol) dropped from high in the jar | Before | After |
| --- | --- | --- |
| Contacts rearranged on landing | 26% | ~15% |
| Peak temperature on impact | ~1220 K | ~1250 K |

That left about half the denting. The impact was the size of a rifle bullet's, and a real iron
bullet mushrooms too, so a melting-point-calibrated metal could not shrug it off entirely.

Second change: gravity. Settling gravity on a solid falling through gas is now 20 instead of 110
(liquids keep 650), so a lump lands at about 280 units/s instead of about 690. A falling lump
takes longer to land, about 17 simulated seconds from the top of the jar instead of 7.

| Iron lump (0.3 mol) dropped from the top of the jar | Before | After |
| --- | --- | --- |
| Contacts rearranged on landing | 19.6% | 1.6% |
| Shape (width x height) | 716x618 -> 1583x508 (spread flat) | 724x631 -> 721x612 |
| Peak temperature on impact | 1451 K | 448 K |

Weakly bound heavy metals (lead, tin, gold, zinc) already felt less settling gravity in
proportion to their cohesion per mass, so they would not be crushed under their own weight. They
now fall more slowly still: lead takes about 2.7 times as long as iron.

Third change (reported: iron still came apart, and even tungsten did): gravity on a falling solid
is now 6 (it was 20), and on liquids 150 (it was 650). Both are still far above real gravity:
about 4×10¹¹ and 10¹³ times. Real gravity in the jar's units is about 2×10⁻¹¹ units/s². So no
breakage seen here comes from gravity weaker than real; at real gravity nothing would ever land
hard enough to dent. The weakly bound metals' reduction no longer slows them below 40% of the
others while falling, so lead lands in about 32 s (at 14% it would take about a minute).

| 1 mol dropped from the top of the jar | Gravity 20 | Gravity 6 |
| --- | --- | --- |
| Iron: contacts rearranged / peak temperature | 3.5% / 479 K | 0.0% / 353 K |
| Tungsten: contacts rearranged / peak temperature | 4.6% / 757 K | 0.4% / 416 K |
| Lead (from the survey below) | | 5.9% (lead is soft) |

### Solids floating on liquids

Found while weakening gravity: liquids and solids felt different settling gravity (650 and 110).
As a result, an iron lump floated on water, sitting on the surface about half submerged, and
lead floated too. Lowering gravity for solids alone would have made this worse. Now a solid that
reaches into a liquid feels the liquid's settling force. The whole lump does, through every coarse
cell it fills that connects to the liquid. Whether it sinks is then decided by density.

| In 1500 water molecules | Before | After |
| --- | --- | --- |
| Iron | Floats at the surface | Sinks to the floor |
| Gold | Floats at the surface | Sinks to the floor |
| Lead | Floats at the surface | Sinks to the floor |
| Mercury (liquid) | Floats on top | Sinks below the surface, not yet to the floor after 16 s |

The weakly bound metals' reduced settling force is now a cap on the force itself (what keeps
lead and liquid mercury from being squeezed), not a fraction of the liquid's. With liquids at
150, lead in water now feels the same pull as the water; mercury keeps the pull it always had,
about 1/4 of water's now, so it only partly sinks. Molten metal does not count as a liquid around a
solid lump, so a solid does not float or sink properly in liquid metal (it falls as in gas).

### Throwing tungsten

Reported: tungsten thrown at the glass broke up just like iron, though it should be far tougher.

Cause: the hand threw everything at up to 300 units/s, so a tungsten atom (3.3 times iron's
mass) hit with 3.3 times the energy, while a metal's bonding, scaled from its melting point, is
only about twice iron's. Thrown that way, tungsten came out worse than iron: 4.5% of contacts
rearranged against 3.0%, and 956 K against 427 K. Giving tungsten's directional bonding a larger
share (from its real hardness, with the pair well lowered to keep its melting point) barely
helped, at 3.8–4.5%.

What changed: the hand now gives each atom at most the motion energy of an iron atom at 300
units/s, as an arm throws a heavier ball of the same size slower.

| 0.3 mol lump, lifted and thrown at the glass | Leaves the hand at | Contacts rearranged | Peak temperature |
| --- | --- | --- | --- |
| Iron | 300 units/s (unchanged) | 3.0% | 427 K |
| Tungsten | 164 units/s (was 300) | 0.0% (was 4.5%) | 287 K (was 956 K) |
| Gold | 158 units/s | 0.0% | 321 K |
| Lead | 154 units/s (was 300) | 0.0% (was 20.5%) | 302 K (was 1437 K) |
| Copper | 282 units/s | 1.4% | 475 K |
| Aluminium | 300 units/s | 1.4% | 354 K |

What's left: damage depends mostly on impact speed (every metal tried stays whole below about
220 units/s and starts to dent above about 250, tungsten included), not on how hard the metal
really is. Lead,
which is very soft, now survives a throw because it is thrown slower; at iron's speed it
dented badly (20.5%).

### Sulfuric acid with something dropped in

Reported: sulfuric acid sat fine on its own, but once something was dropped into it, it reacted
with that and with itself into huge polymer chains.

Checked against real chemistry: pure sulfuric acid holds only a trace of disulfuric acid, and
nothing that reacts with it turns it into long S–O–S chains. The chains were wrong.

| 300 H₂SO₄, 20 s after adding | Before: largest molecule, S–O–S bridges | After |
| --- | --- | --- |
| Nothing | 21 atoms, 5 (some runs ran away to ~190 reactions) | 7 atoms (all 300 intact), 0 |
| 200 water | 14 atoms, 2 | 7 atoms (all intact), 0 |
| 0.1 mol iron | 87 atoms, 46 | 7 atoms, 0 |
| 0.1 mol copper | 21 atoms, 9 | 11 atoms, 1 |
| 0.1 mol zinc | 15 atoms, 6 | 18 atoms, 2 |
| 0.1 mol sodium | 14 atoms, 4 | 8 atoms, 0 |
| On a 650 K plate | 75 atoms | 17 atoms (short polysulfuric acids, as near-boiling acid loses water) |

Causes and fixes:

- **Radicals attacked the acid's sulfur.** Once a metal knocked a hydrogen off a molecule, the
  •O–SO₃H radical bonded to the sulfur of the next molecule (allowed because oxygen is more
  electronegative than hypervalent sulfur) and passed the radical on, so the chain grew. Real
  radicals neither add to the polar S=O bonds of sulfate nor substitute at its crowded sulfur.
  Now a radical does not attack a hypervalent centre, does not add across the oxygen end of a
  hypervalent S=O, and does not displace a group at a saturated centre unless the bond broken is
  weaker than 250 kJ/mol (peroxide O–O, S–S, N–N: H + H₂O₂ → H₂O + OH still runs).
- **A metal broke the wrong bond.** A metal meeting the acid's OH oxygen broke the S–O bond
  (the weakest by homolytic energy), leaving metal hydroxides and SO₃H fragments. An acid's O–H
  gives way as a proton would, so the metal now takes the hydrogen's place.
- **A lone fast molecule read as a hot spot.** A molecule that had just evaporated into a nearly
  empty patch of vapour set its patch's temperature reading, up to ~3000 K. That broke its S–O
  bonds and set off the cascade that seeded the chains even in pure acid. Sparse cells now lean
  harder on the jar's temperature.

Hydrogen and methane flames burn as before (hydrogen: ~907 of 1000 H₂O; before: 902–919). Acid–
base and acid–metal reactions remain missing (issue 31).

### Iron in dilute sulfuric acid

Found from a 1:3 sulfuric acid to water mixture with iron dropped in: the iron floated and did not
react. Two separate causes.

**Floating.** In a liquid, a solid feels the liquid's settling force, so density decides. But in
the flat 2D packing iron comes out only about twice as dense as water (0.048 against 0.023 u per
unit²; really 7.9 times), and the acid mixture at 0.026 left iron only ~1.8 times as dense. The
acid mixture is also about 8–10 times as viscous as water (its molecules wander a third as far,
as in the real 64% acid), so the iron hung at the surface, its centre still above it after 16 s.
An immersed free metal atom is now pulled in proportion to its real density over its 2D one,
against liquids at their real density per LIQ_DENSITY (43.7, from water). Iron is pulled 3.6
times as hard, gold 2.6, lead 1.6, lithium 1.9 and sodium 1.4, which leaves those two lighter than water, as they really are.
Pulled that hard with nothing to slow it, a ball of iron then slammed into water (issue 43), so a
metal in a liquid now sinks no faster than SINK_VMAX, 200 units/s, standing in for the liquid's
drag.

**Not reacting.** Real iron dissolves through protons and electrons: 2 H⁺ + 2 e⁻ → H₂ on its
surface. The model has neither, and the only route was a metal atom taking the acid's OH oxygen
while the H left as a free atom, which costs about 150 kJ/mol. A metal now takes two oxygens at
once while their two protons leave as H₂ (Fe + H₂SO₄ → FeSO₄ + H₂). One proton must be the acid's
and in the water (a water molecule next to its OH: concentrated acid keeps its protons, so cold
concentrated acid does not attack iron, as in real steel tanks); the other may be another acid's
or a water molecule's on the metal (then the product is Fe(OH)(HSO₄), hydrated iron sulfate).
Only metals with a negative standard potential react, with an activation of 52 + 30·E° kJ/mol
(iron 39, zinc 29, tin 48). Metals with oxide films (titanium, zirconium, niobium, molybdenum,
hafnium, tantalum, tungsten) and lead (its sulfate is insoluble) do not. The heat given off is the
real 2FE° (iron 85 kJ/mol): the model's bond energies would give 280 for iron and 398 for zinc,
whose hot spots set off side reactions.

Along the way: a metal atom that had dissolved went on to pull oxygen off the sulfate (giving SO₂
and iron oxide), which takes hot concentrated acid; metals no longer take an oxygen off a
hypervalent centre. Two iron atoms each taking one of the oxygens were left with spare valences
that cross-linked the iron sulfate; the attacking atom now takes both.

| 0.15 mol metal, 16 s at 293 K | Before | After |
| --- | --- | --- |
| Iron, 1:3 acid:water | Floats; no reaction | Sinks; 0.018 mol H₂, iron sulfates |
| Iron, 1:10 | No reaction | 0.010 mol H₂ |
| Zinc, 1:10 | No reaction | 0.011 mol H₂, zinc sulfates (a few ZnOH and H side products) |
| Copper, 1:10 | No reaction | No reaction |
| Iron, concentrated acid | No reaction | No reaction |

The reaction slows as the iron sulfate builds up around the lump (issue 38): the model has no
ions in solution, so the salt does not dissolve away. A screening of the ions' charges by the
water around them was tried; it halved how strongly the sulfate held together but changed little,
the products being slow to diffuse away through the viscous acid, and was left out.

### Checked: water under a sunk lump, gas against the glass

Two things that look wrong and were measured; both are what matter does at this scale (the jar is
24 × 32 nm, an iron "ball" of 0.4 mol a few nm across).

**Water under iron.** A ball of iron resting on the floor of dry glass touches it (its lowest atoms
12–15 units up); in water its lowest point sits 26–31 units up, on one layer of water molecules.
The glass is polar and holds water strongly (hydration layers on silica take very high pressures
to squeeze out), so a monolayer under a nanometre ball is right. The rest of the water under it
is the wedge around the curve of a round ball, which slowly squeezes out (100 water molecules
under it at 6 s, about 55 at 18–24 s).

**Gas between a solid and the glass.** Thrown at the wall, iron reaches the glass (its face at
contact distance) with only about two N₂ molecules in the curved wedge beside the contact, and
bounces off. Held against the glass, the gas there is squeezed out within about a second. A lump
resting against the glass in 1 atm of N₂ has gas in the crevice beside it at the open jar's density
(3.0 against 2.7 atoms per 10⁴ units²), so none is packed in; molecules do linger in the corner,
where the wall and the floor both attract them, as they cling to bare glass (twice the open jar's
density in a layer on the wall). Gas collisions also slow a thrown lump a lot (from 300 to 170 units/s
in a second, stopped within ~500 units): a few-nanometre particle in a gas is in the free-molecule
regime, where that is right. What did turn up: a lump held and rammed into the glass gets very hot
(issue 44).

### Pushing a solid into the glass

Asked for: a solid pushed into the glass should heat up and reshape, as it would, but never melt.
At the hand's top speed a lump carries ~1700 K worth of motion per atom (for iron), and stopped dead
against the glass all of it became heat: iron reached 2900 K in places, lead 1150 K.

- **Cap.** While the hand holds a solid that is within reach of the glass, and for 4 s after it lets
  go of one (thrown), each coarse cell of the body is kept to at most 80% of its melting point (or
  the hot plate's temperature, if higher), by scaling its atoms' motion relative to that part's
  own motion back down. Melting points come from the phase data for molecules, the element for
  free metal atoms and networks, and 1000 K for salts and oxides with no data.
- **Sticky body.** As lead warmed, single atoms read past 601 K for a moment, stopped counting as
  solid and dropped out of the carried body, and with it out of the cap; the hand then dragged
  an ever-hotter patch of melt into the glass (4200 K). Atoms already carried now stay in the
  body unless their surroundings really are past the melting point. Cells the body barely
  reaches into (fewer than 4 atoms) are measured against the body's own motion instead of being
  skipped.
- **Press, not ram.** Once against the glass, the hand kept driving the lump on at full speed, so
  soft lead was squeezed out along the glass so fast that the heat each step outran the cap. A
  body in contact with the glass is now pushed on no faster than 30 units/s; it still hits at
  full speed and dents under the press.

| Pushed into the glass at full speed | Before | After | Melts at |
| --- | --- | --- | --- |
| Iron | 2910 K | ~920 K | 1811 K |
| Lead | 1153 K | ~460 K | 601 K |
| Gold | 1519 K | ~470 K | 1337 K |

(hottest coarse-cell reading of the lump)

### Gases in water and the hydrophobic effect

Issue 42 looked like a missing hydrophobic effect (water pushing non-polar molecules out). Measured:

- **From above, gases do not dissolve.** 0.1 mol of H₂ or N₂ over a pool of water: after 12 s,
  0 and 2 of 200 atoms are in the water, as for real gases (H₂: 1/50 of its gas-phase
  concentration). That was already so.
- **Oil and water already parted, weakly.** Benzene dealt evenly through water rises to the top and
  gathers, but in a small pool only 26% of each carbon's neighbours (from other molecules) were
  benzene. Water's hydrogen bonds did the work; the rest of the hydrophobic effect is an entropy of
  water's network that a pair potential lacks. The dispersion well between a water atom and a
  non-polar atom (partial charge under 0.5, no hydrogen bond, and not an N, O or F carrying an H:
  H₂, N₂, O₂, the C and H of hydrocarbons) is now cut to 30% (HYDROPHOBIC): 54% in the same test, 67% in a larger pool (37%
  before). Ethanol, held to water by its OH, stays mixed (36%). (Counting the O of H₂O₂ as
  non-polar, its charge being small, slowed hot peroxide's decomposition by a third; an N, O or F
  with an H on it now always counts as polar.)
- **A molecule far above its boiling point is a gas even when crowded.** Any molecule packed
  against others felt the liquids' settling force (150 times the gas's), so N₂ that touched a
  pool was pulled toward it. Now a molecule hotter than twice its boiling point (H₂, N₂, O₂, CH₄;
  not NH₃, HCl or CO₂, which dissolve well) never feels it.
- **H₂ made under water stays dissolved for minutes.** It is not held: a lone molecule in water
  only wanders, and leaves by diffusing to the surface, a few units a second at best, so from the
  bottom of a 480-unit pool it takes about a minute. Pulling it up as a bubble would rise (tried,
  up to 2.7 times the liquids' settling force) did nothing: collisions push it about ~50 times
  harder. Real gas escapes in bubbles, but a few molecules cannot make one, and at this scale (a
  pool a few nanometres deep) a real bubble could not exist either: its surface tension would
  squeeze it at hundreds of atmospheres. Contents lists such gas as liquid/solid because it is
  surrounded.

### Boiling, condensing and the plate

Issues 17–19 were measured together (0.5–1 mol, 30 s):

| | Before | After |
| --- | --- | --- |
| Water, insulated jar, 400 K plate | 55–63% liquid, 0.22–0.34 atm | ~54% liquid, ~0.35 atm (real: all vapour, 2.4 atm) |
| Water on a 520 K plate (the test) | 37–39% liquid | 25–28% |
| Steam on cold glass (450 K plate) | 60% liquid, 0.36 atm | ~51%, boiling and condensing on the glass |
| Argon on a 70 K plate | 17% condensed | 58% |
| Nitrogen on a 60 K plate | 57% condensed | 58% |
| Water at 300 K, no air | 72–74% liquid (the test) | 74–75% |
| H₂O₂ on a 750 K plate (the test), 15 s | 0.17 mol water, 0.06 O₂ | 0.11 water, 0.03 O₂: it boils off and its vapour decomposes too slowly (issue 46) |

- **Vapour pressure.** The flat model's liquids evaporate too little as they heat: water's vapour
  pressure was 0.08 atm at 335 K (real 0.16), 0.22 at 380 K (real 1.0), 0.67 at 470 K (real 8.8),
  and only a few molecules left a second, so a liquid far above its boiling point stayed. Now a
  molecule in a liquid or solid hotter than its boiling point at the jar's pressure (the real
  vapour or sublimation curve, from the phase data of ~45 substances; the critical point above
  the critical pressure) by more than 3% has its cohesion faded, over the next 5% of the boiling
  point, to 30%, a step every 64 substeps. The energy that takes is settled as heat, so boiling
  absorbs its latent heat and the plate must supply it. Below the boiling point nothing changes.
- **Which molecules.** Each step matters: faded in the gas, hot N₂ bounced off the cold plate
  (so a molecule out in the gas keeps its cohesion, and one that boiled off stays loose); judged by
  its region's temperature, a thin film on a cold plate read as warm as the gas over it (so a
  liquid's temperature is the lower of its region's and its own); judged by its own alone, random
  jiggling boiled room-temperature water; a pair of molecules meeting in the gas counted as
  liquid (so it takes a region holding packed matter); vapour below its boiling point in the gas
  gets its cohesion back (water boiling into an empty jar overshot), but not at a liquid's surface,
  which evaporation cools below the boiling point (that stopped boiling water boiling); and
  boiling points are taken at 0.2 atm or more (below that the model's own evaporation fills the
  jar: at 0.04 atm, benzene, whose real vapour pressure at room temperature is 0.1 atm, boiled
  through and through in a jar with no air, loosened, and mixed into the water it should float
  on; a 0.04 atm floor left the insulated pool at 46–48% liquid instead of ~54%).
- **Argon** held together too weakly (0.12 atm and a quarter condensed in an insulated jar at
  72 K, real 0.07–0.11 atm); lone argon atoms now get twice the cohesion (krypton and xenon already
  matched).
- **Why the insulated pool still lingers.** In the insulated 400 K jar the liquid sits at
  ~355–370 K, cooled by its own evaporation, and the vapour over it reads only ~0.3–0.4 atm
  where real water would build 2.4 atm and be all vapour: about half of what counts as liquid
  is droplets in the gas, clusters of molecules that boiled off part-loosened and met again.
  Letting droplets boil as well (a looser test for 'in a liquid') changed this by only a few
  percent while costing the 300 K water test and nitrogen's condensing, so it was left out.
- **The plate.** An insulated jar reaches the plate's temperature; otherwise the glass loses heat
  to the room and the gas settles in between, as in a real jar. Stronger plate coupling for liquids
  (5x) was tried and changed little.

### Drawing only what is on screen

Asked whether zooming in skips drawing what is not in view. Measured in Chromium (software
rendering) with ~8,900 atoms in a 1400×900 window:

- **Atom and mole views** already drew only the atoms on screen. Zooming in makes a frame
  cheaper: the atom view takes 8.3 ms at its widest, 4.9 ms at 2.2×, 2.7 ms at 6× and 1.1 ms
  fully zoomed in; the mole view 6.8 ms over the whole jar, 2.1 ms looking at a corner. Each frame
  still looks at every atom once to see whether it is on screen, which is cheap next to drawing.
  The edge margin was too small at high zoom (issue 48).
- **The real-world view** built the whole jar each frame (issue 47). It now builds the window on
  screen plus the blur's reach (about 100 screen px), snapped to whole cells so that a pan moves
  the eased picture along instead of restarting it. It looks the same: compared pixel by pixel
  with the old drawing at eight camera positions (whole jar, zoomed in on its middle, floor,
  corner and right wall, and the mole and atom views), no pixel differs by more than 2/255. The
  real-view part of a frame, before → after:

  | View | 1400×900 | 400×800 (phone) |
  | --- | --- | --- |
  | Whole jar (fit) | ~5.5 ms, unchanged (the same work: CPU profiles match function by function) | ~2–3 ms, unchanged |
  | Zoomed in on the floor of the jar | 12.3 → 8.0 ms | 10.0 → 4.6 ms |
  | Zoomed in on the top corner | 10.7 → 7.9 ms | 10.1 → 3.0 ms |

- **The simulation is not culled**: every atom keeps moving and reacting wherever the camera
  is. Drawing less leaves more of each frame for it (the main loop gives physics what drawing
  leaves over), so a zoomed-in view runs the jar faster.
- Seen while measuring: between the mole and real-world views (scale 0.19–0.25) both are drawn
  and cross-faded, which made it the slowest view (~27 ms a frame here). It is a short range of
  zoom, passed through rather than stayed in, so it was left (issue 49).

### Every element dropped

Reported: drop each element on its own and check that it moves and lands properly. A script
dropped 1 mol of each element, as the Add tool does, near the top of the jar at room temperature.
It recorded the landing time and impact speed, contacts rearranged on landing (solids), peak and
final temperature, and how much stayed condensed.

| | First run (gravity change only) | After the fixes below |
| --- | --- | --- |
| Solid elements that land whole (≤10% rearranged, no overheating) | 80 of 96 | 86 of 96 |
| Hottest element | oganesson 24,000 K, selenium 4609 K | germanium 1122 K for a moment on impact |
| Never landed | silicon, germanium, arsenic, antimony, caesium, francium | rubidium, caesium, francium (see below) |

What it found and what changed:

- **Covalent network elements flashed to 10⁴–10⁵ K** (boron, silicon, phosphorus, sulfur,
  germanium, arsenic, selenium, antimony, tellurium). They were laid down as a heap of loose
  atoms 1.4–1.8 bond lengths apart, collapsed onto one another, and bonded and re-bonded in a
  frenzy. An earlier fix (a reaction giving off heat ends the settling-in thermostat, so flames
  are not quenched) also switched that thermostat off for the atoms assembling. Now these elements
  go down already bonded in their 2D natural form: a square net (silicon, germanium), a
  honeycomb (boron, phosphorus, arsenic, antimony; carbon as graphite), or zigzag chains joined
  into long loops (sulfur, selenium, tellurium). Loops, because chain ends are radicals that pass
  bonds along without end. A small flat ring like S8 was tried too: bent far from these atoms'
  preferred 92°, it flashed hot, and its few contacts with the next rings could not hold them
  together. While newly added, a substance also does not react with itself.
- **Silicon hung in mid-air.** Energy owed for endothermic events was repaid by slowing the
  atoms of a cell toward rest, not by cooling them about their own motion, so it braked a
  falling lump. It is now paid out of motion about the flow of the packed matter in each cell. A
  gas is still cooled toward rest: a first version also cooled gas about its own mean motion,
  which cannot cool a lone fast atom, and iron burning in oxygen ran away to 10⁵ K.
- **A settling lump wandered.** The strong heat exchange that holds new material at room
  temperature gave random kicks that did not cancel. They walked a small lump of tellurium up to
  the ceiling at tens of units/s. Each cell's net kick is now removed.
- **Sodium, potassium and gallium splatted.** A few surface atoms jiggling past the melting
  point (both are close to room temperature) counted as melt. That gave them, and through the
  buoyancy rule the whole lump, liquid-strength gravity. A metal atom now switches between solid
  and liquid gravity only when the matter around it is clearly past the melting point (by 12%),
  and molten metal never counts as a liquid that wets a lump. Sodium went from 43% rearranged
  to 5%.
- **Sulfur floated like a gas** as S8 rings: an atom counted as condensed only if it or a bonded
  neighbour touched something. The whole of a molecule now counts.
- **Lumps did not fit.** 1 mol of a large-atom element is wider than the space above the drop
  point. What did not fit was scattered at random through the jar. A lump is now moved in from
  the glass until it fits.
- **Astatine reached up to 44,000 K, oganesson 24,000 K.** Astatine and tennessine came as loose
  atoms that paired up hot; they are now added as diatomic molecules, like iodine, packed side by
  side. Molecules of an element with no melting point data of their own (At2) take the
  element's, so they count as solid. Oganesson (a noble element that melts above room
  temperature) now counts as a solid too.

What's left (the 10 that do not land well):

- **Rubidium, caesium and francium** (melting points 300–312 K) and **potassium** (337 K) are
  barely below their melting points. They come in lumps as wide as the jar (their atoms are the
  largest), fall slowly, and do not land within 60 s (potassium just does). **Gallium** (303 K)
  half-melts and splashes.
- **Iodine, astatine, tennessine and oganesson**, held together only by dispersion forces, fall
  apart: in 2D a molecule has too few neighbours to hold the crystal together (iodine was already
  listed above).
- Sulfur, selenium and tellurium chains read 500–570 K after landing, falling slowly; silicon
  and germanium, almost fully rigid in 2D, read too cold.

---

## Summary of the built-in experiments

| Preset | Real world | What the sandbox does (after the fix pass) | Verdict |
| --- | --- | --- | --- |
| Hydrogen + oxygen (spark) | Explodes to water | Burns to about 90% water. Flame peaks around 3000 K (real ~3000 K) **(README)** | Good |
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
**Status: Fixed.** Sodium took water's hydrogen one atom at a time, a radical step with a high
cost, so it reacted slowly and the H atoms it set free stuck to the metal as NaH, while sodium
pulled the H off the NaOH it had made, leaving Na₂O. Calcium took water's H as CaH and freed OH.
The concerted step that dissolves metals in acid (issue 37) now takes both protons from water for
metals with a standard potential of −2.5 V or below (WATER_E0): 2 Na + 2 H₂O → 2 NaOH + H₂, with
the heat 2F(E° + 0.41) (444 kJ per H₂ for sodium). These metals no longer take a single H off an
O–H (from either end). A hydride of one of them next to water is hydrolysed at once
(CaH₂ + 2 H₂O → Ca(OH)₂ + 2 H₂), and H₂ splits onto a metal only with an extra 60 kJ/mol
(HYDRIDE_BARRIER): lithium, sodium and calcium take up hydrogen only at a few hundred degrees.

| 0.06 mol in water, 10 s | Before | After |
| --- | --- | --- |
| Sodium | NaOH 0.005, NaH 0.005, Na₂O 0.005, H₂ 0.005, H atoms | NaOH 0.018, H₂ 0.009 |
| Potassium | KH 0.007, K₂O 0.003, KOH 0.001 | KOH 0.016, H₂ 0.008 |
| Lithium | LiH 0.010, Li₂O 0.006, LiOH 0.004, H₂ 0.002, H atoms | LiOH 0.020, H₂ 0.010 |
| Calcium | Ca₂H₂O 0.008, CaH₂ 0.005, CaH₂O 0.004, H atoms, Ca(OH)₂ 0.002 | Ca(OH)₂ 0.014, H₂ 0.014, traces of CaHO/CaH₂ |
| Magnesium, iron | no reaction | no reaction |

The reaction still runs only where the metal touches water: its hydroxide does not dissolve away
(no ions in solution), and a simulated second is a few milliseconds of chemistry, so a lump takes
many simulated seconds to go.

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
**Status: Fixed (partly).** The *Steam on cold glass* preset was still nearly all liquid after
30 s on a 450 K plate. A liquid more than 3% above its boiling point at the jar's pressure (from
the real vapour curve) now loses cohesion and boils off, taking its latent heat from the plate:
the preset is ~51% liquid at 30 s (was 60%), water on the test's 520 K plate 25–28% (was 37–39%).
An insulated pool on a 400 K plate still lingers (~54% liquid; real water would all be vapour).
See "Boiling, condensing and the plate". Sulfuric acid and bromine were not re-checked.

### 18. The jar runs far below the plate — Medium
**Status: Explained; cold end fixed.**
- **Hot end:** a 3000 K plate holds an N₂/O₂ jar at ~2050 K. An insulated jar reaches the plate
  (~3000–3200 K): the gap is the glass losing heat to the 293 K room, as a real jar's would.
- **Cold end:** 60–70 K plates held the jar at 100–140 K, so argon never liquefied on the 70 K
  plate (the noble-gas preset asks you to try this). Argon's atoms held together too weakly; with
  twice the cohesion it is 58% condensed after 30 s (was 17%). The jar's reading stays above the
  plate because it averages in the gas.

### 19. Condensed matter at impossible temperatures — Medium
**Status: Fixed / explained.**
- **Liquid nitrogen preset:** N₂ was reported as condensed while the jar read 131–141 K, above
  nitrogen's 126 K critical temperature. The liquid's own temperature is ~77 K; the reading is
  the jar's average, mostly the warmer gas over it. No liquid can now last above its critical
  temperature, or above its boiling point.
- **Water:** stayed "mostly liquid" with the jar at 250–260 K for 30 s. It does freeze (62% ice at
  15 s, over 80% from 30 s); the old reading counted ice as condensed.

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
