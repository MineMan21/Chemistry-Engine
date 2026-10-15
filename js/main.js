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
    tool: 'drop',
    perClick: 1,
    paused: false,
    speed: 1,
    natural: true,
    pinned: -1, // atom whose molecule is pinned in the inspector
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
    cells[z].classList.add('selected');
    CE.UI.renderElementInfo($('el-info'), z, state.natural);
    setTool('drop');
  }

  function setTool(t) {
    state.tool = t;
    document.querySelectorAll('#tools button').forEach(function (b) { b.classList.toggle('active', b.dataset.tool === t); });
    canvas.style.cursor = t === 'inspect' ? 'help' : 'crosshair';
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
  $('plate').addEventListener('input', function (e) { setPlate(sliderToT(+e.target.value)); });
  document.querySelectorAll('#plate-presets button').forEach(function (b) {
    b.addEventListener('click', function () { setPlate(+b.dataset.t); });
  });
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
    state.pinned = -1;
    sim.analyze();
    refreshContents();
  });
  $('zoom-in').addEventListener('click', function () { renderer.zoomAt(renderer.width / 2, renderer.height / 2, 1.6); });
  $('zoom-out').addEventListener('click', function () { renderer.zoomAt(renderer.width / 2, renderer.height / 2, 1 / 1.6); });
  $('zoom-fit').addEventListener('click', function () { renderer.fit(); });

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
    water: { text: 'Liquid water held together by hydrogen bonds. Turn the plate up past 373 K to boil it.', run: function () {
      setPlate(293); sim.addMolecules(TEMPLATES.H2O, 900, P.W / 2, P.H - 200, { x: P.W / 2 - 30, y: 190 }); } },
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
      setPlate(60); sim.addMoles(7, 1, P.W / 2, P.H * 0.5); } }
  };
  $('preset-go').addEventListener('click', function () {
    const p = PRESETS[$('preset').value];
    if (!p) return;
    sim.clear();
    state.pinned = -1;
    sim.vacuumOn = sim.canOn = false;
    refreshPressureButtons();
    p.run();
    sim.analyze();
    refreshContents();
    toast(p.text);
  });

  // ------------------------------------------------------------------ adding
  function insideJar(x, y) { return x > 0 && x < P.W && y > 0 && y < P.H; }

  function addAtomsAt(x, y) {
    const z = state.element;
    const before = sim.n;
    const res = sim.addSubstance(z, state.perClick, x, y, { natural: state.natural, spacing: E[z].radius * (state.natural && CE.DIATOMIC.has(z) ? 3.4 : 2.3) });
    afterAdd(before, res.units, state.perClick);
  }

  function addMolesAt(x, y) {
    const mol = Math.max(0.001, Math.min(60, +$('mol-amount').value || 1));
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
    if (pointers.size === 2) { drag = { pinch: pinchInfo() }; return; }
    drag = { x: e.offsetX, y: e.offsetY, camX: renderer.cam.x, camY: renderer.cam.y, moved: false, button: e.button };
  });
  canvas.addEventListener('pointermove', function (e) {
    state.mouse = { x: e.offsetX, y: e.offsetY };
    if (pointers.has(e.pointerId)) pointers.set(e.pointerId, { x: e.offsetX, y: e.offsetY });
    if (!drag) return;
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
    if (e.key === ' ') { setPaused(!state.paused); e.preventDefault(); }
    else if (e.key === '+' || e.key === '=') renderer.zoomAt(renderer.width / 2, renderer.height / 2, 1.4);
    else if (e.key === '-' || e.key === '_') renderer.zoomAt(renderer.width / 2, renderer.height / 2, 1 / 1.4);
    else if (e.key === 'f' || e.key === 'F') renderer.fit();
  });

  // ------------------------------------------------------------------ panels
  function refreshContents() {
    CE.UI.renderContents($('contents'), sim.species, function (sp) {
      if (sp.first < sim.n) { state.pinned = sp.first; setTool('inspect'); }
    });
  }

  function updateModeUI() {
    const atomView = renderer.isAtomView();
    $('mode-badge').textContent = atomView ? 'Atom view · click to drop atoms' : 'Bulk view · click to add moles';
    $('add-atoms').classList.toggle('active', atomView);
    $('add-moles').classList.toggle('active', !atomView);
    $('ro-temp').textContent = Math.round(sim.temperature) + ' K';
    const p = sim.pressureAtm;
    $('ro-press').textContent = (p < 0.01 ? p.toExponential(1) : p < 100 ? p.toFixed(2) : p.toFixed(0)) + ' atm';
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
      sim.spark(b.x, b.y, b.hot ? 380 : 220, b.hot ? 4000 : 80);
      if (--b.frames === 0) state.burst = null;
    }
    if (!state.paused && sim.n > 0) {
      // adaptive sub-stepping: spend ~9 ms per frame on physics (scaled by the speed setting)
      const budget = 9 * state.speed;
      const subs = Math.max(1, Math.min(16, Math.floor(budget / Math.max(0.05, state.msPerSub))));
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
      if (state.frame % 30 === 0) refreshContents();
    }
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
