// App wiring: input, controls, presets and the main loop.
(function (CE) {
  'use strict';

  const P = CE.Simulation.params;
  const E = CE.elements;
  const $ = function (id) { return document.getElementById(id); };

  const sim = new CE.Simulation();
  const canvas = $('view');
  const renderer = new CE.Renderer(canvas, sim);

  const state = {
    element: 1,
    mol: null, // selected molecule template (null = add the selected element)
    tool: 'drop',
    perClick: 1,
    paused: false,
    speed: 1,
    natural: true,
    pinned: -1, // atom whose molecule is pinned in the inspector
    focusKey: null, // species picked in Contents: all its molecules are highlighted
    hover: -1,
    mouse: null, // {x, y} in screen px
    msPerSub: 1,
    frame: 0
  };

  // ------------------------------------------------------------------ controls
  const cells = CE.UI.buildPeriodicTable($('ptable'), selectElement);
  function selectElement(z) {
    if (cells[state.element]) cells[state.element].classList.remove('selected');
    state.element = z;
    state.mol = null;
    cells[z].classList.add('selected');
    CE.UI.renderElementInfo($('el-info'), z, state.natural);
    addPhaseButton($('el-info'), { element: z });
    renderMoleculeLists();
    setTool('drop');
  }

  // ------------------------------------------------------------------ molecules
  let saved = CE.molecules.loadSaved();
  function selectMolecule(m) {
    if (cells[state.element]) cells[state.element].classList.remove('selected');
    state.mol = m;
    CE.UI.renderMoleculeInfo($('el-info'), m);
    addPhaseButton($('el-info'), { key: moleculeKey(m.atoms), formula: m.formula, name: m.name });
    renderMoleculeLists();
    setTool('drop');
  }
  function renderMoleculeLists() {
    const sel = state.mol ? state.mol.id : null;
    CE.UI.renderMoleculeChips($('mol-presets'), CE.molecules.PRESETS, sel, selectMolecule, null);
    CE.UI.renderMoleculeChips($('mol-saved'), saved, sel, selectMolecule, function (m) {
      saved = saved.filter(function (o) { return o.id !== m.id; });
      CE.molecules.storeSaved(saved);
      if (state.mol && state.mol.id === m.id) selectElement(state.element);
      renderMoleculeLists();
      toast('Deleted ' + m.name + ' from My molecules.');
    });
  }
  function savePinned() {
    const i = state.pinned;
    if (i < 0 || i >= sim.n) { toast('Pick Inspect and click a molecule in the jar first.'); return; }
    const tpl = sim.captureMolecule(i, 150);
    if (!tpl) { toast('That cluster is too big to save (150 atoms at most). Try a single molecule.'); return; }
    const info = sim.inspect(i);
    const sp = sim.species && sim.spIdx ? sim.species[sim.spIdx[i]] : null;
    const existing = saved.find(function (m) { return m.formula === info.formula && m.atoms.length === tpl.atoms.length; });
    if (existing) { selectMolecule(existing); toast(existing.name + ' is already in My molecules, so it is selected now.'); return; }
    const m = {
      id: 'saved-' + Date.now(), saved: true,
      name: info.name || info.formula, formula: info.formula,
      gas: !!(sp && sp.phase === 'gas' && tpl.atoms.length <= 12),
      atoms: tpl.atoms, bonds: tpl.bonds
    };
    saved.push(m);
    const ok = CE.molecules.storeSaved(saved);
    selectMolecule(m);
    toast('Saved ' + m.name + ' to My molecules' + (ok ? '' : ' for this visit (this browser blocks storage)') + '. Click in the jar to add it.');
  }
  $('save-pinned').addEventListener('click', savePinned);

  const TOOL_HINTS = {
    drop: 'Drop: click in the jar to add the selected element or molecule.',
    spark: 'Spark: click to set off a spark. It heats a small spot to 4000 K and splits some bonds, which ignites fuels such as H₂ + O₂ or CH₄ + O₂.',
    cool: 'Chill: click to blast a small spot down to 80 K, to freeze or condense what is there.',
    inspect: 'Inspect: click a molecule to pin its details. Then press Save this molecule to keep it in My molecules.',
    grab: 'Grab: press on material and drag to pull it around. A solid lifts as one lump, a liquid stretches and drips, a gas slips away.'
  };
  function setTool(t) {
    state.tool = t;
    document.querySelectorAll('#tools button').forEach(function (b) { b.classList.toggle('active', b.dataset.tool === t); });
    canvas.style.cursor = t === 'inspect' ? 'help' : t === 'grab' ? 'grab' : 'crosshair';
    $('tool-hint').textContent = TOOL_HINTS[t];
  }
  document.querySelectorAll('#tools button').forEach(function (b) {
    b.addEventListener('click', function () { setTool(b.dataset.tool); });
  });
  document.querySelectorAll('#per-click button').forEach(function (b) {
    b.addEventListener('click', function () {
      state.perClick = +b.dataset.n;
      document.querySelectorAll('#per-click button').forEach(function (o) { o.classList.toggle('active', o === b); });
    });
  });
  $('natural').addEventListener('change', function (e) {
    state.natural = e.target.checked;
    CE.UI.renderElementInfo($('el-info'), state.element, state.natural);
    addPhaseButton($('el-info'), { element: state.element });
  });
  $('add-btn').addEventListener('click', function () {
    addMolesAt(P.W / 2, P.H * 0.25);
  });
  $('scale-note').textContent = 'Each particle stands for ' + (P.MOL_PER_PARTICLE * 1000) +
    ' mmol (6.0×10²⁰ atoms), so 1 mol = ' + Math.round(1 / P.MOL_PER_PARTICLE) + ' particles. The jar holds 24 L.';

  // hot plate: logarithmic slider 20 K .. 4000 K
  const TMIN = 20, TMAX = 4000;
  function sliderToT(v) { return TMIN * Math.pow(TMAX / TMIN, v / 1000); }
  function tToSlider(T) { return Math.round(1000 * Math.log(T / TMIN) / Math.log(TMAX / TMIN)); }
  function setPlate(T) {
    sim.plateTemp = Math.max(TMIN, Math.min(TMAX, T));
    $('plate').value = tToSlider(sim.plateTemp);
    $('plate-out').textContent = Math.round(sim.plateTemp) + ' K (' + Math.round(sim.plateTemp - 273.15) + ' °C)';
  }
  $('plate').addEventListener('input', function (e) { program = null; setPlate(sliderToT(+e.target.value)); });
  document.querySelectorAll('#plate-presets button').forEach(function (b) {
    b.addEventListener('click', function () { program = null; setPlate(+b.dataset.t); });
  });

  // A plate program for experiments: stages of { T, sec, text } run in simulated time (so a slow
  // computer just takes longer). Touching the plate control cancels it.
  let program = null;
  function runProgram(stages) { program = { stages: stages, i: -1, until: 0 }; nextStage(); }
  function nextStage() {
    const p = program;
    if (++p.i >= p.stages.length) { program = null; return; }
    const st = p.stages[p.i];
    setPlate(st.T);
    if (st.text) toast(st.text);
    p.until = st.sec ? sim.time + st.sec : Infinity;
  }
  $('insulated').addEventListener('change', function (e) { sim.insulated = e.target.checked; });

  function refreshPressureButtons() {
    const v = $('vacuum'), c = $('can');
    v.textContent = 'Vacuum pump: ' + (sim.vacuumOn ? 'on' : 'off');
    v.classList.toggle('on', sim.vacuumOn); v.classList.toggle('vac', sim.vacuumOn);
    c.textContent = 'Gas can: ' + (sim.canOn ? 'open' : 'shut');
    c.classList.toggle('on', sim.canOn);
  }
  $('vacuum').addEventListener('click', function () {
    sim.vacuumOn = !sim.vacuumOn;
    if (sim.vacuumOn) sim.canOn = false;
    refreshPressureButtons();
  });
  $('can').addEventListener('click', function () {
    sim.canOn = !sim.canOn;
    if (sim.canOn) sim.vacuumOn = false;
    refreshPressureButtons();
  });
  // pump strength: logarithmic from a gentle 2%/s to a hard 95%/s of the gas removed per second
  function setPump(v) {
    const frac = 0.02 * Math.pow(0.95 / 0.02, v / 100);
    sim.pumpRate = -Math.log(1 - frac);
    $('pump-rate-out').textContent = Math.round(frac * 100) + '%/s';
  }
  $('pump-rate').addEventListener('input', function (e) { setPump(+e.target.value); });
  setPump(+$('pump-rate').value);
  $('can-gas').addEventListener('change', function (e) { sim.canGas = e.target.value; });
  $('can-rate').addEventListener('input', function (e) {
    sim.canRate = +e.target.value;
    $('can-rate-out').textContent = sim.canRate.toFixed(2) + ' mol/s';
  });

  function setPaused(p) {
    state.paused = p;
    $('play').textContent = p ? 'Resume' : 'Pause';
  }
  $('play').addEventListener('click', function () { setPaused(!state.paused); });
  $('speed').addEventListener('change', function (e) { state.speed = +e.target.value; });
  $('color-mode').addEventListener('change', function (e) { renderer.colorMode = e.target.value; });
  $('clear').addEventListener('click', function () {
    sim.clear();
    program = null;
    state.pinned = -1;
    state.focusKey = null;
    sim.analyze();
    refreshContents();
  });
  $('zoom-in').addEventListener('click', function () { renderer.zoomAt(renderer.width / 2, renderer.height / 2, 1.6); });
  $('zoom-out').addEventListener('click', function () { renderer.zoomAt(renderer.width / 2, renderer.height / 2, 1 / 1.6); });
  $('zoom-fit').addEventListener('click', function () { renderer.fit(); });
  // zoom-level switcher: jump to the atoms, moles or real-world view around the current spot
  document.querySelectorAll('#tiers button').forEach(function (b) {
    b.addEventListener('click', function () {
      const target = CE.Renderer.TIER_SCALE[b.dataset.tier];
      if (!target) { renderer.fit(); return; }
      const cx = Math.min(P.W - 50, Math.max(50, renderer.cam.x));
      let cy = Math.min(P.H - 50, Math.max(50, renderer.cam.y));
      if (b.dataset.tier === 'atoms' && renderer.realMix() > 0) {
        // coming from far out: look where the material is (the lowest stuff in the jar)
        cy = P.H - 0.4 * renderer.height / target;
      }
      renderer.cam.scale = target; renderer.cam.x = cx; renderer.cam.y = Math.max(cy, 50);
    });
  });

  // ------------------------------------------------------------------ presets
  const TEMPLATES = {
    H2O: { atoms: [[8, 0, 0], [1, 9.7, 0], [1, -2.43, 9.39]], bonds: [[0, 1, 1], [0, 2, 1]] },
    CH4: { atoms: [[6, 0, 0], [1, 10.7, 0], [1, -10.7, 0], [1, 0, 10.7], [1, 0, -10.7]], bonds: [[0, 1, 1], [0, 2, 1], [0, 3, 1], [0, 4, 1]] }
  };
  const PRESETS = {
    h2o2: { text: 'H₂ + O₂ is stable at room temperature. Use Spark (or heat the plate) to ignite it.', run: function () {
      setPlate(293); sim.addMoles(1, 1, P.W / 2, P.H * 0.35); sim.addMoles(8, 0.5, P.W / 2, P.H * 0.65); } },
    nacl: { text: 'Sodium metal reacts with chlorine gas: watch Na–Cl ionic bonds and salt crystals form.', run: function () {
      setPlate(293); sim.addMoles(11, 0.5, P.W / 2, P.H * 0.85); sim.addMoles(17, 0.3, P.W / 2, P.H * 0.3); } },
    water: { text: 'A blob of liquid water dropped from mid-air: it falls, splashes into a pool and wets the glass. Heat the plate past 373 K to boil it.', run: function () {
      setPlate(293); sim.addMolecules(TEMPLATES.H2O, 1400, P.W / 2, P.H * 0.3, { x: 330, y: 260 }); } },
    steam: { text: 'Water on a 450 K plate boils; the vapour condenses on the room-temperature glass and runs down the walls.', run: function () {
      setPlate(450); sim.addMolecules(TEMPLATES.H2O, 1200, P.W / 2, P.H - 170, { x: P.W / 2 - 30, y: 160 }); } },
    methane: { text: 'CH₄ + 2 O₂. Spark it to burn methane into CO₂ and H₂O.', run: function () {
      setPlate(293); sim.addMolecules(TEMPLATES.CH4, 250, P.W / 2, P.H / 2, { x: P.W / 2 - 60, y: P.H / 2 - 60 });
      sim.addMoles(8, 0.5, P.W / 2, P.H / 2); } },
    rust: { text: 'Iron is heated in oxygen; iron oxides form on the hot metal.', run: function () {
      setPlate(1100); sim.addMoles(26, 0.6, P.W / 2, P.H * 0.9); sim.addMoles(8, 0.6, P.W / 2, P.H * 0.35); } },
    carbon: { text: 'Carbon atoms link into a solid network. Heat the plate high to burn it to CO / CO₂.', run: function () {
      setPlate(293); sim.addMoles(6, 0.3, P.W / 2, P.H * 0.9); sim.addMoles(8, 0.6, P.W / 2, P.H * 0.35); } },
    noble: { text: 'Helium, neon and argon never bond. Cool them down: argon liquefies near 87 K.', run: function () {
      setPlate(293); sim.addMoles(2, 0.5, P.W / 2, P.H * 0.3); sim.addMoles(10, 0.5, P.W / 2, P.H * 0.5); sim.addMoles(18, 0.5, P.W / 2, P.H * 0.7); } },
    nitrogen: { text: 'Nitrogen below its 77 K boiling point condenses into a liquid at the bottom.', run: function () {
      setPlate(60); sim.addMoles(7, 1, P.W / 2, P.H * 0.5); } },
    ice: { text: 'Liquid water on a 200 K plate freezes into ice Ih. Zoom in on the solid to see its hexagonal rings.', run: function () {
      sim.addMolecules(TEMPLATES.H2O, 900, P.W / 2, P.H - 260, { x: 420, y: 220 });
      runProgram([
        { T: 293, sec: 3 },
        { T: 200, sec: 0, text: 'Plate at 200 K: the water freezes from the cold floor and glass up. Ice is an open honeycomb of hydrogen bonds, so the solid takes up about 9% more room than the liquid. Set the plate above 273 K to melt it again.' }
      ]); } },
    iron: { text: 'An iron lump on a 2100 K plate melts (iron melts at 1811 K); then the plate goes cold and the melt solidifies into close-packed crystal grains.', run: function () {
      sim.addMoles(26, 0.5, P.W / 2, P.H - 150);
      runProgram([
        { T: 293, sec: 2 },
        { T: 2100, sec: 14, text: 'Plate at 2100 K: the iron heats past 1811 K and melts into a liquid that flows across the plate.' },
        { T: 293, sec: 0, text: 'Plate back to 293 K: the molten iron cools, crystal grains nucleate on the cold plate and the metal hardens.' }
      ]); } },
    gallium: { text: 'Gallium melts at just 303 K (it melts in your hand). The 340 K plate melts it; then the plate cools it to 250 K and it freezes.', run: function () {
      sim.addMoles(31, 0.5, P.W / 2, P.H - 150);
      runProgram([
        { T: 293, sec: 3 },
        { T: 340, sec: 12, text: 'Plate at 340 K, just above gallium\'s 303 K melting point: the solid slumps into a liquid.' },
        { T: 250, sec: 0, text: 'Plate at 250 K: the liquid gallium freezes again.' }
      ]); } }
  };
  $('preset-go').addEventListener('click', function () {
    const p = PRESETS[$('preset').value];
    if (!p) return;
    sim.clear();
    program = null;
    setPlate(293);
    state.pinned = -1;
    state.focusKey = null;
    sim.vacuumOn = sim.canOn = false;
    refreshPressureButtons();
    p.run();
    sim.analyze();
    refreshContents();
    toast(p.text);
  });

  // ------------------------------------------------------------------ adding
  function insideJar(x, y) { return x > 0 && x < P.W && y > 0 && y < P.H; }

  // Adds `count` copies of the selected molecule around (x, y): packed like a liquid/solid, or
  // spread out for gases.
  function addMoleculesAt(x, y, count) {
    const m = state.mol;
    let size = 0;
    m.atoms.forEach(function (a) { size = Math.max(size, Math.hypot(a[1], a[2]) + E[a[0]].radius * 0.5); });
    const room = Math.floor((P.MAX - sim.n) / m.atoms.length);
    const n = Math.min(count, room);
    if (n <= 0) { toast('The jar is full (' + P.MAX + ' particles). Empty it or use the vacuum pump.'); return; }
    const spacing = m.gas ? 3 : 1;
    const a0 = size * 2.05 * spacing;
    let half = Math.sqrt(n / 0.87) * a0 / 2 + a0;
    // a square patch, stretched downward/upward when wider than the jar
    const hx = Math.min(P.W / 2 - 10, half), hy = Math.min(P.H / 2 - 10, half * half / hx);
    const cx = Math.min(P.W - hx - 5, Math.max(hx + 5, x)), cy = Math.min(P.H - hy - 5, Math.max(hy + 5, y));
    const before = sim.n;
    const placed = sim.addMolecules(m, n, cx, cy, { x: hx, y: hy }, P.AMBIENT, spacing);
    afterAdd(before, placed, count);
  }

  function addAtomsAt(x, y) {
    if (state.mol) { addMoleculesAt(x, y, state.perClick); return; }
    const z = state.element;
    const before = sim.n;
    const res = sim.addSubstance(z, state.perClick, x, y, { natural: state.natural, spacing: E[z].radius * (state.natural && CE.DIATOMIC.has(z) ? 3.4 : 2.3) });
    afterAdd(before, res.units, state.perClick);
  }

  function addMolesAt(x, y) {
    const mol = Math.max(0.001, Math.min(60, +$('mol-amount').value || 1));
    if (state.mol) { addMoleculesAt(x, y, Math.max(1, Math.round(mol / P.MOL_PER_PARTICLE))); return; }
    const z = state.element;
    const before = sim.n;
    const units = Math.max(1, Math.round(mol / P.MOL_PER_PARTICLE));
    const res = sim.addSubstance(z, units, x, y, { natural: state.natural });
    afterAdd(before, res.units, units);
  }

  function afterAdd(before, placed, requested) {
    if (sim.n >= P.MAX - 1) toast('The jar is full (' + P.MAX + ' particles). Empty it or use the vacuum pump.');
    else if (placed < requested) toast('Only room for ' + placed + ' of ' + requested + ' — the jar is crowded there.');
    if (sim.n !== before) { sim.analyze(); refreshContents(); }
  }

  let toastTimer = 0;
  function toast(msg) {
    const t = $('toast');
    t.textContent = msg; t.hidden = false;
    clearTimeout(toastTimer);
    toastTimer = setTimeout(function () { t.hidden = true; }, 5000);
  }

  // ------------------------------------------------------------------ pointer input
  let drag = null;
  const pointers = new Map();
  canvas.addEventListener('pointerdown', function (e) {
    canvas.setPointerCapture(e.pointerId);
    pointers.set(e.pointerId, { x: e.offsetX, y: e.offsetY });
    if (pointers.size === 2) { sim.grabEnd(); drag = { pinch: pinchInfo() }; return; }
    if (state.tool === 'grab' && e.button === 0) {
      // grab what is under the cursor (about 50 screen pixels, at least a few dozen atoms)
      const w = renderer.screenToWorld(e.offsetX, e.offsetY);
      if (sim.grabStart(w[0], w[1], Math.min(600, Math.max(120, 50 / renderer.cam.scale)))) {
        drag = { grab: true };
        canvas.style.cursor = 'grabbing';
        return;
      }
    }
    drag = { x: e.offsetX, y: e.offsetY, camX: renderer.cam.x, camY: renderer.cam.y, moved: false, button: e.button };
  });
  canvas.addEventListener('pointermove', function (e) {
    state.mouse = { x: e.offsetX, y: e.offsetY };
    if (pointers.has(e.pointerId)) pointers.set(e.pointerId, { x: e.offsetX, y: e.offsetY });
    if (!drag) return;
    if (drag.grab) {
      const w = renderer.screenToWorld(e.offsetX, e.offsetY);
      sim.grabMove(w[0], w[1]);
      return;
    }
    if (drag.pinch && pointers.size === 2) {
      const now = pinchInfo();
      renderer.zoomAt(now.cx, now.cy, now.d / drag.pinch.d);
      renderer.cam.x -= (now.cx - drag.pinch.cx) / renderer.cam.scale;
      renderer.cam.y -= (now.cy - drag.pinch.cy) / renderer.cam.scale;
      drag.pinch = now;
      return;
    }
    if (drag.pinch) return;
    const dx = e.offsetX - drag.x, dy = e.offsetY - drag.y;
    if (!drag.moved && dx * dx + dy * dy > 25) drag.moved = true;
    if (drag.moved) {
      renderer.cam.x = drag.camX - dx / renderer.cam.scale;
      renderer.cam.y = drag.camY - dy / renderer.cam.scale;
    }
  });
  function endPointer(e) {
    pointers.delete(e.pointerId);
    if (!drag) return;
    if (drag.grab) { sim.grabEnd(); drag = null; canvas.style.cursor = 'grab'; return; }
    if (drag.pinch) { if (pointers.size === 0) drag = null; return; }
    if (!drag.moved && drag.button === 0 && e.type === 'pointerup') click(e.offsetX, e.offsetY);
    drag = null;
  }
  canvas.addEventListener('pointerup', endPointer);
  canvas.addEventListener('pointercancel', endPointer);
  canvas.addEventListener('pointerleave', function () { state.mouse = null; });
  canvas.addEventListener('contextmenu', function (e) { e.preventDefault(); });
  canvas.addEventListener('wheel', function (e) {
    e.preventDefault();
    renderer.zoomAt(e.offsetX, e.offsetY, Math.exp(-e.deltaY * 0.0015));
  }, { passive: false });

  function pinchInfo() {
    const pts = Array.from(pointers.values());
    return { cx: (pts[0].x + pts[1].x) / 2, cy: (pts[0].y + pts[1].y) / 2, d: Math.hypot(pts[0].x - pts[1].x, pts[0].y - pts[1].y) || 1 };
  }

  function click(sx, sy) {
    const w = renderer.screenToWorld(sx, sy);
    if (state.tool === 'spark' || state.tool === 'cool') {
      if (!insideJar(w[0], w[1])) return;
      // a short burst (several frames) rather than a single instant, like a real igniter
      state.burst = { x: w[0], y: w[1], frames: 10, hot: state.tool === 'spark' };
      if (state.tool === 'spark') renderer.addSpark(w[0], w[1]);
      return;
    }
    if (state.tool === 'inspect') {
      sim.buildGrid();
      const i = sim.nearest(w[0], w[1], Math.max(30, 12 / renderer.cam.scale));
      state.pinned = i;
      return;
    }
    if (!insideJar(w[0], w[1])) return;
    if (renderer.isAtomView()) addAtomsAt(w[0], w[1]);
    else addMolesAt(w[0], w[1]);
  }

  window.addEventListener('keydown', function (e) {
    if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
    if (!$('phase-dialog').hidden) { if (e.key === 'Escape') closePhase(); return; }
    if (e.key === ' ') { setPaused(!state.paused); e.preventDefault(); }
    else if (e.key === '+' || e.key === '=') renderer.zoomAt(renderer.width / 2, renderer.height / 2, 1.4);
    else if (e.key === '-' || e.key === '_') renderer.zoomAt(renderer.width / 2, renderer.height / 2, 1 / 1.4);
    else if (e.key === 'f' || e.key === 'F') renderer.fit();
    else if ((e.key === 's' || e.key === 'S') && state.pinned >= 0) savePinned();
    else if (e.key === 'Escape' && state.focusKey) { state.focusKey = null; refreshContents(); }
  });

  // ------------------------------------------------------------------ panels
  function refreshContents() {
    // a highlighted species that is all gone (reacted, pumped out) stops being highlighted
    if (state.focusKey && !sim.species.some(function (s) { return s.key === state.focusKey; })) state.focusKey = null;
    CE.UI.renderContents($('contents'), sim.species, function (sp) {
      state.focusKey = state.focusKey === sp.key ? null : sp.key;
      refreshContents();
      if (state.focusKey) toast('Highlighting every ' + sp.formula + (sp.name ? ' (' + sp.name + ')' : '') + ' in the jar. Click it again to clear.');
    }, state.focusKey);
    const fb = $('phase-focus'), sp = state.focusKey && sim.species.find(function (x) { return x.key === state.focusKey; });
    fb.hidden = !sp;
    if (sp) {
      fb.textContent = 'Phase diagram of ' + sp.formula;
      fb.onclick = function () { openPhase({ key: sp.key, formula: sp.formula, name: sp.name }); };
    }
  }

  // ------------------------------------------------------------------ phase diagrams
  const phaseUI = { chart: null, key: null, returnFocus: null };
  function moleculeKey(atoms) {
    const counts = {};
    atoms.forEach(function (a) { counts[a[0]] = (counts[a[0]] || 0) + 1; });
    return CE.chem.describe(counts, atoms.length, false).key;
  }
  function addPhaseButton(container, what) {
    const b = document.createElement('button');
    b.className = 'btn phase-btn';
    b.textContent = 'Phase diagram';
    b.addEventListener('click', function () { openPhase(what); });
    container.appendChild(b);
  }
  // what: { element: z } or { key, formula, name }
  function openPhase(what) {
    let d, title, missing;
    if (what.element) {
      const e = CE.elements[what.element];
      d = CE.phase.forElement(what.element);
      const form = CE.DIATOMIC.has(what.element) ? CE.chem.pretty(e.symbol + '2') : e.symbol;
      title = e.name + ' (' + form + ')';
      missing = 'melts at ' + (e.mp ? Math.round(e.mp) + ' K' : 'an unknown temperature') + ' and boils at ' + (e.bp ? Math.round(e.bp) + ' K' : 'an unknown temperature') + ' at 1 atm';
    } else {
      d = CE.phase.forKey(what.key);
      title = (d ? d.name : what.name || what.formula) + ' (' + what.formula + ')';
    }
    const body = $('phase-body');
    body.innerHTML = '';
    $('phase-title').textContent = 'Phase diagram: ' + title;
    phaseUI.chart = null;
    phaseUI.key = d ? d.key : null;
    if (!d) {
      $('phase-sub').textContent = 'No documented phase diagram';
      const p = document.createElement('p');
      p.className = 'phase-note';
      p.textContent = 'Its triple and critical points are not well documented, so no diagram is drawn' +
        (missing ? '. What is known: it ' + missing + '.' : '.');
      body.appendChild(p);
    } else {
      $('phase-sub').textContent = 'Which phase it is in at each temperature and pressure. Hover the chart to read values; the ring marks the jar right now.';
      phaseUI.chart = CE.PhaseChart.render(body, d);
      body.appendChild(CE.PhaseChart.table(d));
      const note = document.createElement('p');
      note.className = 'phase-note';
      note.textContent = (d.note ? d.note + ' ' : '') + 'Reference data: NIST Chemistry WebBook and CRC Handbook. Lines between the documented points follow the Clausius–Clapeyron relation' +
        (d.est.length ? '; values marked (est.) are estimates.' : '.');
      body.appendChild(note);
      updatePhaseMarker();
    }
    const dlg = $('phase-dialog');
    if (dlg.hidden) phaseUI.returnFocus = document.activeElement;
    dlg.hidden = false;
    $('phase-close').focus();
  }
  function closePhase() {
    $('phase-dialog').hidden = true;
    phaseUI.chart = null;
    if (phaseUI.returnFocus && phaseUI.returnFocus.focus) phaseUI.returnFocus.focus();
  }
  // the jar's pressure, and the substance's own temperature when it is in the jar
  function updatePhaseMarker() {
    if (!phaseUI.chart) return;
    const sp = sim.species && sim.species.find(function (x) { return x.key === phaseUI.key; });
    phaseUI.chart.setJar(sp ? sp.T : sim.temperature, sim.n ? sim.pressureAtm : 0);
  }
  $('phase-close').addEventListener('click', closePhase);
  $('phase-dialog').addEventListener('click', function (e) { if (e.target === e.currentTarget) closePhase(); });

  // Per-atom mask of the highlighted species, from the species index of the last analysis.
  function updateFocus() {
    const k = state.focusKey ? sim.species.findIndex(function (s) { return s.key === state.focusKey; }) : -1;
    if (k < 0 || !sim.spIdx) { renderer.focus = null; return; }
    const mask = state.focusMask || (state.focusMask = new Uint8Array(P.MAX));
    const sp = sim.spIdx;
    for (let i = 0; i < sim.n; i++) mask[i] = sp[i] === k ? 1 : 0;
    renderer.focus = mask;
  }

  function updateModeUI() {
    const atomView = renderer.isAtomView(), tier = renderer.tier();
    $('mode-badge').textContent = tier === 'atoms' ? 'Click to drop atoms' : tier === 'moles' ? 'Each dot is 1 mmol · click to add moles' : 'As it looks to the eye · click to add moles';
    document.querySelectorAll('#tiers button').forEach(function (b) { b.classList.toggle('active', b.dataset.tier === tier); });
    $('add-atoms').classList.toggle('active', atomView);
    $('add-moles').classList.toggle('active', !atomView);
    $('ro-temp').textContent = Math.round(sim.temperature) + ' K';
    const p = sim.pressureAtm;
    $('ro-press').textContent = (p > 0 && p < 0.01 ? p.toExponential(1) : p < 100 ? p.toFixed(2) : p.toFixed(0)) + ' atm';
    $('ro-atoms').textContent = sim.n.toLocaleString();
    $('ro-rx').textContent = Math.round(sim.reactionsPerSec).toLocaleString() + '/s';
  }

  function updateInspector() {
    const box = $('inspector');
    let i = state.pinned >= 0 && state.pinned < sim.n ? state.pinned : -1;
    const pinned = i >= 0;
    if (!pinned && state.mouse && !drag) {
      const w = renderer.screenToWorld(state.mouse.x, state.mouse.y);
      if (insideJar(w[0], w[1])) i = sim.nearest(w[0], w[1], Math.max(25, 8 / renderer.cam.scale));
    }
    $('save-pinned').hidden = !pinned;
    if (i < 0) { box.hidden = true; renderer.highlight = null; return; }
    const info = sim.inspect(i);
    if (!info) { box.hidden = true; renderer.highlight = null; return; }
    renderer.highlight = new Set(sim.moleculeOf(i, []));
    CE.UI.renderInspector(box, info, pinned);
    box.hidden = false;
    // place near the atom on screen, kept inside the stage
    const s = renderer.worldToScreen(sim.x[i], sim.y[i]);
    const bw = box.offsetWidth, bh = box.offsetHeight;
    let left = s[0] + 18, top = s[1] + 18;
    if (left + bw > renderer.width - 8) left = s[0] - bw - 18;
    if (top + bh > renderer.height - 8) top = s[1] - bh - 18;
    box.style.left = Math.max(8, left) + 'px';
    box.style.top = Math.max(8, top) + 'px';
  }

  // ------------------------------------------------------------------ main loop
  let last = performance.now();
  function frame(now) {
    const dt = Math.min(0.1, (now - last) / 1000);
    last = now;
    if (state.burst && state.burst.frames > 0 && !state.paused) {
      const b = state.burst;
      sim.spark(b.x, b.y, b.hot ? 550 : 260, b.hot ? 4000 : 80);
      if (--b.frames === 0) state.burst = null;
    }
    if (program && sim.time >= program.until) nextStage();
    if (!state.paused && sim.n > 0) {
      // adaptive sub-stepping: spend ~9 ms per frame on physics (scaled by the speed setting)
      const budget = 12 * state.speed;
      const subs = Math.max(1, Math.min(24, Math.floor(budget / Math.max(0.05, state.msPerSub))));
      const t0 = performance.now();
      sim.step(subs);
      const per = (performance.now() - t0) / subs;
      state.msPerSub = state.msPerSub * 0.9 + per * 0.1;
    } else if (!state.paused) {
      sim.step(1); // keep gauges / gas can alive in an empty jar
    }
    state.frame++;
    if (state.frame % 10 === 0) {
      sim.analyze();
      if (state.frame % 30 === 0) { refreshContents(); updatePhaseMarker(); }
    }
    updateFocus();
    renderer.draw(dt);
    updateModeUI();
    updateInspector();
    requestAnimationFrame(frame);
  }

  function onResize() { renderer.resize(); }
  window.addEventListener('resize', onResize);

  // ------------------------------------------------------------------ boot
  renderer.resize();
  renderer.fit();
  selectElement(1);
  setPlate(293);
  refreshPressureButtons();
  sim.analyze();
  refreshContents();
  toast('Pick an element, then click in the jar. Zoom in (scroll) to see atoms bond; zoom out to work in moles.');
  requestAnimationFrame(frame);

  CE.app = { sim: sim, renderer: renderer, state: state };
})(window.CE = window.CE || {});
