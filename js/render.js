// Canvas renderer: the apparatus (jar, hot plate, vacuum pump, pressure can) plus particles.
// Zoomed in: ball-and-stick atoms with bonds, symbols and charges.
// Zoomed out: every atom becomes a pixel coloured by its molecule (or element).
(function (CE) {
  'use strict';

  const E = CE.elements;
  const P = CE.Simulation.params;
  const NZ = CE.chem.NZ;

  // Zoom tiers (camera scale = screen px per world unit, 1 unit = 10 pm):
  //   atoms      scale >= ATOM_SCALE   ball-and-stick atoms and bonds
  //   moles      in between            one pixel per atom, coloured by molecule
  //   real world scale <= REAL_FULL    substances drawn as they look (liquid, solid, tinted gas)
  // The mole and real-world views cross-fade between REAL_FULL and REAL_START.
  const ATOM_SCALE = 0.64, REAL_START = 0.25, REAL_FULL = 0.19;
  const ATOM_VIEW_PX = 7 * ATOM_SCALE;
  // In the real-world view the jar is drawn as a 24 cm wide jar: 1 world unit = 0.1 mm.
  const REAL_UNITS_PER_CM = 100;

  function Renderer(canvas, sim) {
    this.canvas = canvas;
    this.ctx = canvas.getContext('2d');
    this.sim = sim;
    this.dpr = 1;
    this.width = 0; this.height = 0;
    this.cam = { x: P.W / 2, y: P.H / 2 - 300, scale: 0.1 };
    this.colorMode = 'molecule';
    this.hover = -1;
    this.highlight = null; // Set of atom indices to outline
    this.focus = null; // per-atom mask (1 = part of the species picked in Contents); others dim
    this.showCharges = true;
    this.pixelGroups = new Map(); // mole view: squares to draw, by colour
    this.pixelCss = new Map();
    this.time = 0;
    this.sparks = [];
    this.lightCanvas = document.createElement('canvas');
    this.lightCtx = this.lightCanvas.getContext('2d');
    this.elColorCss = E.map(function (e) { return e ? e.color : '#fff'; });
    this.elTextColor = E.map(function (e) {
      if (!e) return '#000';
      const r = parseInt(e.color.substr(1, 2), 16), g = parseInt(e.color.substr(3, 2), 16), b = parseInt(e.color.substr(5, 2), 16);
      return (0.299 * r + 0.587 * g + 0.114 * b) > 140 ? '#111' : '#fff';
    });
  }
  const R = Renderer.prototype;
  Renderer.ATOM_VIEW_PX = ATOM_VIEW_PX;

  R.resize = function () {
    const dpr = window.devicePixelRatio || 1;
    const rect = this.canvas.getBoundingClientRect();
    this.dpr = dpr;
    this.width = Math.max(1, rect.width);
    this.height = Math.max(1, rect.height);
    this.canvas.width = Math.round(this.width * dpr);
    this.canvas.height = Math.round(this.height * dpr);
  };

  // World bounds of the whole apparatus (for "fit").
  R.fit = function () {
    // apparatus spans x: -1100 .. W + 800, y: -1700 .. H + 420
    const x0 = -1150, x1 = P.W + 850, y0 = -1750, y1 = P.H + 450;
    this.cam.scale = Math.min(this.width / (x1 - x0), this.height / (y1 - y0), REAL_FULL) * 0.97;
    this.cam.x = (x0 + x1) / 2;
    this.cam.y = (y0 + y1) / 2;
  };
  R.minScale = function () { return Math.min(this.width / (P.W + 6000), this.height / (P.H + 6000)); };
  R.maxScale = function () { return 14; };

  R.worldToScreen = function (x, y) {
    return [(x - this.cam.x) * this.cam.scale + this.width / 2, (y - this.cam.y) * this.cam.scale + this.height / 2];
  };
  R.screenToWorld = function (sx, sy) {
    return [(sx - this.width / 2) / this.cam.scale + this.cam.x, (sy - this.height / 2) / this.cam.scale + this.cam.y];
  };
  R.zoomAt = function (sx, sy, factor) {
    const before = this.screenToWorld(sx, sy);
    this.cam.scale = Math.min(this.maxScale(), Math.max(this.minScale(), this.cam.scale * factor));
    const after = this.screenToWorld(sx, sy);
    this.cam.x += before[0] - after[0];
    this.cam.y += before[1] - after[1];
  };
  // Atom view when a carbon-sized atom would be at least a few pixels across.
  R.isAtomView = function () { return this.cam.scale >= ATOM_SCALE; };
  // 0 = mole (pixel) view, 1 = real-world view, in between while cross-fading
  R.realMix = function () { return Math.min(1, Math.max(0, (REAL_START - this.cam.scale) / (REAL_START - REAL_FULL))); };
  R.tier = function () { return this.isAtomView() ? 'atoms' : this.realMix() >= 0.5 ? 'real' : 'moles'; };
  Renderer.TIER_SCALE = { atoms: 2.2, moles: 0.38, real: 0 }; // 0 = fit the apparatus

  R.addSpark = function (x, y) { this.sparks.push({ x: x, y: y, t: 0 }); };

  R.draw = function (dt) {
    this.time += dt;
    const ctx = this.ctx, dpr = this.dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    const css = getComputedStyle(document.documentElement);
    this.theme = {
      bg: css.getPropertyValue('--canvas-bg').trim() || '#0d1117',
      glass: css.getPropertyValue('--glass').trim() || 'rgba(170,210,255,0.5)',
      metal: css.getPropertyValue('--metal').trim() || '#5b6573',
      text: css.getPropertyValue('--canvas-text').trim() || '#cfd8e3',
      interior: css.getPropertyValue('--jar-interior').trim() || 'rgba(120,170,230,0.05)'
    };
    ctx.fillStyle = this.theme.bg;
    ctx.fillRect(0, 0, this.width, this.height);

    ctx.save();
    ctx.translate(this.width / 2, this.height / 2);
    ctx.scale(this.cam.scale, this.cam.scale);
    ctx.translate(-this.cam.x, -this.cam.y);
    this.drawApparatusBack(ctx);
    ctx.restore();

    if (this.isAtomView()) this.drawAtoms(ctx);
    else {
      const mix = this.realMix();
      if (mix < 1) this.drawPixels(ctx, 1 - mix);
      if (mix > 0) this.drawReal(ctx, mix, dt);
      if (mix > 0 && this.focus) this.drawFocusDots(ctx, mix);
    }

    this.drawLight(ctx, dt);

    ctx.save();
    ctx.translate(this.width / 2, this.height / 2);
    ctx.scale(this.cam.scale, this.cam.scale);
    ctx.translate(-this.cam.x, -this.cam.y);
    this.drawApparatusFront(ctx);
    this.drawSparks(ctx, dt);
    this.drawGrab(ctx);
    ctx.restore();
    this.drawScaleBar(ctx);
  };

  function plateColor(T) {
    // black-body-ish glow for hot, frosty blue for cold
    if (T < 273) {
      const f = Math.min(1, (273 - T) / 250);
      return 'rgb(' + Math.round(70 - 30 * f) + ',' + Math.round(90 + 60 * f) + ',' + Math.round(110 + 120 * f) + ')';
    }
    const f = Math.min(1, Math.max(0, (T - 350) / 2500));
    const r = Math.round(70 + 185 * Math.min(1, f * 2.2));
    const g = Math.round(70 + 150 * Math.max(0, f * 1.6 - 0.45));
    const b = Math.round(75 + 150 * Math.max(0, f * 2 - 1.2));
    return 'rgb(' + r + ',' + g + ',' + b + ')';
  }

  R.drawApparatusBack = function (ctx) {
    const W = P.W, H = P.H, t = this.theme;
    // jar interior tint
    ctx.fillStyle = t.interior;
    ctx.fillRect(0, 0, W, H);
    // hot plate body
    const plateT = this.sim.plateTemp;
    ctx.fillStyle = t.metal;
    roundRect(ctx, -260, H + 40, W + 520, 340, 40); ctx.fill();
    ctx.fillStyle = plateColor(plateT);
    roundRect(ctx, -160, H + 40, W + 320, 90, 20); ctx.fill();
    if (plateT > 500) {
      const glow = ctx.createLinearGradient(0, H + 40, 0, H - 500);
      const a = Math.min(0.35, (plateT - 500) / 6000);
      glow.addColorStop(0, 'rgba(255,120,40,' + a + ')');
      glow.addColorStop(1, 'rgba(255,120,40,0)');
      ctx.fillStyle = glow;
      ctx.fillRect(0, H - 500, W, 540);
    } else if (plateT < 250) {
      const glow = ctx.createLinearGradient(0, H + 40, 0, H - 400);
      const a = Math.min(0.3, (250 - plateT) / 700);
      glow.addColorStop(0, 'rgba(140,200,255,' + a + ')');
      glow.addColorStop(1, 'rgba(140,200,255,0)');
      ctx.fillStyle = glow;
      ctx.fillRect(0, H - 400, W, 440);
    }
    // knob + label on the hot plate
    ctx.fillStyle = '#20252c';
    ctx.beginPath(); ctx.arc(W / 2, H + 260, 70, 0, Math.PI * 2); ctx.fill();
    const ang = -Math.PI * 0.75 + Math.min(1, plateT / 4000) * Math.PI * 1.5;
    ctx.strokeStyle = '#ffb347'; ctx.lineWidth = 18;
    ctx.beginPath(); ctx.moveTo(W / 2, H + 260); ctx.lineTo(W / 2 + Math.sin(ang) * 60, H + 260 - Math.cos(ang) * 60); ctx.stroke();
    ctx.fillStyle = t.text;
    ctx.font = '120px system-ui, sans-serif';
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    ctx.fillText('HOT PLATE  ' + Math.round(plateT) + ' K', W / 2 - 120, H + 260);
    ctx.textAlign = 'left';
    ctx.fillText('JAR  ' + Math.round(this.sim.temperature) + ' K', W / 2 + 120, H + 260);
  };

  R.drawApparatusFront = function (ctx) {
    const W = P.W, H = P.H, t = this.theme, sim = this.sim;
    // glass walls
    ctx.strokeStyle = t.glass;
    ctx.lineWidth = 36;
    ctx.beginPath();
    ctx.moveTo(-18, -18); ctx.lineTo(-18, H + 18); ctx.lineTo(W + 18, H + 18); ctx.lineTo(W + 18, -18);
    ctx.stroke();
    // lid with valve opening
    const vl = W / 2 - P.VALVE_HALF, vr = W / 2 + P.VALVE_HALF;
    ctx.fillStyle = t.metal;
    ctx.fillRect(-60, -110, vl + 60, 92);
    ctx.fillRect(vr, -110, W - vr + 60, 92);
    // valve manifold and pipes
    ctx.fillRect(vl - 40, -420, vr - vl + 80, 310);
    ctx.fillStyle = sim.vacuumOn ? '#4fc3f7' : sim.canOn ? '#9ccc65' : '#2a3038';
    ctx.fillRect(vl, -300, vr - vl, 190);
    // pipe to the vacuum pump (left) and pressure can (right)
    ctx.strokeStyle = t.metal; ctx.lineWidth = 70; ctx.lineCap = 'round';
    ctx.beginPath(); ctx.moveTo(vl - 20, -340); ctx.lineTo(-500, -340); ctx.lineTo(-500, -620); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(vr + 20, -340); ctx.lineTo(W + 500, -340); ctx.lineTo(W + 500, -620); ctx.stroke();
    ctx.lineCap = 'butt';
    // vacuum pump
    ctx.fillStyle = '#39424d';
    roundRect(ctx, -1050, -1250, 1100, 650, 60); ctx.fill();
    ctx.fillStyle = sim.vacuumOn ? '#4fc3f7' : '#59636f';
    ctx.beginPath(); ctx.arc(-500, -925, 210, 0, Math.PI * 2); ctx.fill();
    ctx.strokeStyle = '#1b2026'; ctx.lineWidth = 34;
    const spin = sim.vacuumOn ? this.time * 8 : 0;
    for (let k = 0; k < 3; k++) {
      const a = spin + k * Math.PI * 2 / 3;
      ctx.beginPath(); ctx.moveTo(-500, -925); ctx.lineTo(-500 + Math.cos(a) * 180, -925 + Math.sin(a) * 180); ctx.stroke();
    }
    ctx.fillStyle = t.text; ctx.font = '110px system-ui, sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'alphabetic';
    ctx.fillText('VACUUM ' + (sim.vacuumOn ? 'ON' : 'OFF'), -500, -1300);
    // pressure can
    ctx.fillStyle = '#7d8a99';
    roundRect(ctx, W + 250, -1500, 500, 880, 120); ctx.fill();
    ctx.fillStyle = sim.canOn ? '#9ccc65' : '#4b5560';
    ctx.fillRect(W + 250, -1150, 500, 160);
    ctx.fillStyle = '#111'; ctx.font = 'bold 120px system-ui, sans-serif';
    ctx.fillText(CE.chem.pretty(sim.canGas === 'air' ? 'AIR' : sim.canGas), W + 500, -1030);
    ctx.fillStyle = t.text; ctx.font = '110px system-ui, sans-serif';
    ctx.fillText('GAS CAN ' + (sim.canOn ? 'OPEN' : 'SHUT'), W + 500, -1560);
    // pressure gauge on the manifold
    const gx = W / 2, gy = -640;
    ctx.fillStyle = '#e8edf2'; ctx.beginPath(); ctx.arc(gx, gy, 200, 0, Math.PI * 2); ctx.fill();
    ctx.strokeStyle = '#39424d'; ctx.lineWidth = 26; ctx.stroke();
    const p = sim.pressureAtm;
    const frac = Math.min(1, Math.log10(1 + p * 9) / Math.log10(1 + 9 * 20));
    const ga = -Math.PI * 1.25 + frac * Math.PI * 1.5;
    ctx.strokeStyle = '#d32f2f'; ctx.lineWidth = 16;
    ctx.beginPath(); ctx.moveTo(gx, gy); ctx.lineTo(gx + Math.cos(ga) * 160, gy + Math.sin(ga) * 160); ctx.stroke();
    ctx.fillStyle = '#111'; ctx.font = 'bold 80px system-ui, sans-serif';
    ctx.fillText(p < 10 ? p.toFixed(2) : p.toFixed(0), gx, gy + 120);
    ctx.font = '60px system-ui, sans-serif';
    ctx.fillText('atm', gx, gy + 175);
    ctx.fillStyle = '#39424d'; ctx.fillRect(gx - 30, gy + 200, 60, 220);
    // gas jets
    if (sim.canOn || sim.vacuumOn) {
      ctx.fillStyle = sim.canOn ? 'rgba(160,220,120,0.6)' : 'rgba(80,190,250,0.6)';
      for (let k = 0; k < 14; k++) {
        const ph = (this.time * 2 + k / 14) % 1;
        const yy = sim.canOn ? -300 + ph * 380 : 80 - ph * 380;
        ctx.beginPath(); ctx.arc(vl + ((k * 53) % (vr - vl)), yy, 14, 0, Math.PI * 2); ctx.fill();
      }
    }
  };

  R.drawSparks = function (ctx, dt) {
    for (let k = this.sparks.length - 1; k >= 0; k--) {
      const s = this.sparks[k];
      s.t += dt;
      if (s.t > 0.6) { this.sparks.splice(k, 1); continue; }
      const a = 1 - s.t / 0.6;
      ctx.strokeStyle = 'rgba(255,220,120,' + a + ')';
      ctx.lineWidth = 6 / this.cam.scale;
      for (let r = 0; r < 8; r++) {
        const ang = r * Math.PI / 4 + s.t * 3;
        const r0 = 20 + s.t * 200, r1 = r0 + 60 + 40 * a;
        ctx.beginPath(); ctx.moveTo(s.x + Math.cos(ang) * r0, s.y + Math.sin(ang) * r0);
        ctx.lineTo(s.x + Math.cos(ang) * r1, s.y + Math.sin(ang) * r1); ctx.stroke();
      }
    }
  };

  // Zoomed out: one pixel (or small square) per atom.
  // Mole view: each atom a small square. Drawn as rectangles grouped by colour (a handful of
  // fill colours a frame), not as a full-screen image rebuilt and uploaded to the graphics card
  // every frame, which on a high-resolution screen made this the laggiest view.
  R.drawPixels = function (ctx, alpha) {
    const sim = this.sim, n = sim.n;
    const w = this.width * this.dpr, h = this.height * this.dpr;
    const s = this.cam.scale * this.dpr;
    const ox = this.width / 2 * this.dpr - this.cam.x * s, oy = this.height / 2 * this.dpr - this.cam.y * s;
    const colors = this.colorMode === 'element' ? null : sim.molColor;
    const elColor = CE.Simulation.elColor;
    const x = sim.x, y = sim.y, el = sim.el;
    const hl = this.highlight, focus = this.focus, groups = this.pixelGroups;
    if (groups.size > 512) groups.clear();
    groups.forEach(function (list) { list.length = 0; });
    for (let i = 0; i < n; i++) {
      const px = Math.floor(x[i] * s + ox), py = Math.floor(y[i] * s + oy);
      const size = Math.max(1, Math.round(E[el[i]].radius * 1.2 * s));
      if (px < -size || py < -size || px >= w + size || py >= h + size) continue;
      let col = colors ? colors[i] : elColor[el[i]];
      if (focus) col = focus[i] ? 0xff00ffff : (col & 0x00ffffff) | 0x30000000; // picked: yellow, rest faint
      if (hl && hl.has(i)) col = 0xff00ffff; // yellow
      let list = groups.get(col);
      if (!list) { list = []; groups.set(col, list); }
      list.push(px - (size >> 1), py - (size >> 1), size);
    }
    ctx.save();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.globalAlpha = alpha === undefined ? 1 : alpha;
    const css = this.pixelCss;
    groups.forEach(function (list, col) {
      if (!list.length) return;
      let fill = css.get(col);
      if (!fill) {
        // (colours are packed as in image data: alpha, blue, green, red from the high byte down)
        fill = 'rgba(' + (col & 255) + ',' + ((col >>> 8) & 255) + ',' + ((col >>> 16) & 255) + ',' + ((col >>> 24) / 255).toFixed(3) + ')';
        css.set(col, fill);
      }
      ctx.fillStyle = fill;
      for (let k = 0; k < list.length; k += 3) ctx.fillRect(list[k], list[k + 1], list[k + 2], list[k + 2]);
    });
    ctx.restore();
  };

  // Real-world view: atoms are splatted onto a coarse grid, blurred into continuous material and
  // shaded by what the substance looks like (translucent liquid with a bright surface, opaque
  // solid, tinted gas). The simulation underneath is unchanged. Everything that decides how an
  // atom is drawn varies smoothly (splatting, condensed and solid weights) and the grid is eased
  // over a fraction of a second, so the material does not flicker as hot atoms jostle.
  R.drawReal = function (ctx, alpha, dt) {
    const sim = this.sim, n = sim.n, W = P.W, H = P.H;
    const cw = Math.max(6, 2.6 / this.cam.scale); // grid cell in world units (~2.6 screen px)
    // gas is a soft haze: drawn from a grid GF times coarser, scaled up underneath the rest
    const GF = 4, cg4 = cw * GF;
    const R0 = Math.max(1, Math.round(30 / cw)), RG = Math.max(1, Math.round(160 / cg4));
    // Only the part of the jar on screen is built (zoomed in or panned, the rest is skipped), plus
    // a margin as wide as the blur reaches, so the picture at the screen's edge is the same as if
    // the whole jar were built. The window starts on a whole gas cell, so the fine and the gas
    // grids line up and a pan moves both by whole cells; it ends at the glass where the jar does
    // (as the whole-jar grid did: run on past it, the blur thinned liquid lying on the floor).
    const tl = this.screenToWorld(0, 0), br = this.screenToWorld(this.width, this.height);
    const pad = Math.max(2 * R0 + 2, (2 * RG + 2) * GF) * cw;
    const wx0 = Math.max(0, Math.floor((tl[0] - pad) / cg4)) * cg4, wy0 = Math.max(0, Math.floor((tl[1] - pad) / cg4)) * cg4;
    const wx1 = Math.min(W, br[0] + pad), wy1 = Math.min(H, br[1] + pad);
    if (wx1 <= wx0 || wy1 <= wy0) { this._rg = null; return; } // the jar is off screen
    const gw = Math.ceil((wx1 - wx0) / cw), gh = Math.ceil((wy1 - wy0) / cw), N = gw * gh;
    const ggw = Math.ceil((wx1 - wx0) / cg4), ggh = Math.ceil((wy1 - wy0) / cg4), GN = ggw * ggh;
    const mk = function (k, len) { const a = []; for (let q = 0; q < k; q++) a.push(new Float32Array(len)); return a; };
    let rg = this._rg;
    if (!rg || rg.cw !== cw) {
      // a new zoom level starts afresh
      rg = this._rg = { cw: cw, gw: gw, gh: gh, ggw: ggw, ggh: ggh, x0: wx0, y0: wy0, ch: mk(7, N), prev: mk(7, N),
        gch: mk(4, GN), gprev: mk(4, GN), fresh: true, tmp: new Float32Array(N), gtmp: new Float32Array(GN) };
      this.realCanvas = null;
    } else if (rg.x0 !== wx0 || rg.y0 !== wy0 || rg.gw !== gw || rg.gh !== gh) {
      // a pan keeps the eased picture, moved with the window (newly uncovered cells start empty,
      // in the margin off screen)
      const dx = Math.round((wx0 - rg.x0) / cw), dy = Math.round((wy0 - rg.y0) / cw);
      const gdx = Math.round((wx0 - rg.x0) / cg4), gdy = Math.round((wy0 - rg.y0) / cg4);
      const prev = mk(7, N), gprev = mk(4, GN);
      for (let k = 0; k < 7; k++) remapGrid(rg.prev[k], rg.gw, rg.gh, prev[k], gw, gh, dx, dy);
      for (let k = 0; k < 4; k++) remapGrid(rg.gprev[k], rg.ggw, rg.ggh, gprev[k], ggw, ggh, gdx, gdy);
      if (rg.gw !== gw || rg.gh !== gh) { rg.ch = mk(7, N); rg.gch = mk(4, GN); rg.tmp = new Float32Array(N); rg.gtmp = new Float32Array(GN); this.realCanvas = null; }
      Object.assign(rg, { gw: gw, gh: gh, ggw: ggw, ggh: ggh, x0: wx0, y0: wy0, prev: prev, gprev: gprev });
    }
    if (!this.realCanvas) {
      this.realCanvas = document.createElement('canvas');
      this.realCanvas.width = gw; this.realCanvas.height = gh;
      this.realCtx = this.realCanvas.getContext('2d');
      this.realImg = this.realCtx.createImageData(gw, gh);
      this.gasCanvas = document.createElement('canvas');
      this.gasCanvas.width = ggw; this.gasCanvas.height = ggh;
      this.gasCtx = this.gasCanvas.getContext('2d');
      this.gasImg = this.gasCtx.createImageData(ggw, ggh);
    }
    // channels: density, colour (3), opacity, solid, heat; gas: strength, colour (3)
    const ch = rg.ch, gch = rg.gch;
    for (let k = 0; k < 7; k++) ch[k].fill(0);
    for (let k = 0; k < 4; k++) gch[k].fill(0);
    const [dens, cr, cg, cb, cop, csol, heat] = ch;
    const [ga, gr, gg, gb] = gch;
    const tLoc = sim.tLoc;
    const species = sim.species || [], spIdx = sim.spIdx, x = sim.x, y = sim.y, el = sim.el;
    const area = 1 / (cw * cw);
    const condF = sim.condF;
    // (atoms off the window are skipped; at the glass, as before, an atom right at it is drawn in
    // the edge cell)
    const sx0 = wx0 > 0 ? wx0 : -Infinity, sy0 = wy0 > 0 ? wy0 : -Infinity;
    const sx1 = wx1 < W ? wx1 : Infinity, sy1 = wy1 < H ? wy1 : Infinity;
    let anyCond = false, anyGas = false;
    for (let i = 0; i < n; i++) {
      if (x[i] < sx0 || x[i] >= sx1 || y[i] < sy0 || y[i] >= sy1) continue; // off screen
      const sp = spIdx ? species[spIdx[i]] : null;
      const look = sp && sp.look;
      if (!look) continue;
      // bilinear splat over the four nearest cells, so an atom moving across a cell boundary
      // fades across instead of jumping
      let fx = (x[i] - wx0) / cw - 0.5, fy = (y[i] - wy0) / cw - 0.5;
      fx = fx < 0 ? 0 : fx > gw - 1.001 ? gw - 1.001 : fx;
      fy = fy < 0 ? 0 : fy > gh - 1.001 ? gh - 1.001 : fy;
      const gx = Math.floor(fx), gy = Math.floor(fy), tx = fx - gx, ty = fy - gy;
      const c = gy * gw + gx;
      const r = E[el[i]].radius;
      const w = Math.PI * r * r * area;
      // how much of the atom is drawn as liquid/solid: a smooth ramp in its time-smoothed
      // 'packed' flag (an instantaneous contact count would flicker), so a few gas molecules
      // that happen to touch stay invisible vapour
      const wc = smooth(0.35, 0.75, condF[i]);
      // solid or liquid: crystal-forming substances (water, metals) by the atom's own state, so
      // a melt near its melting point does not flip between looks with the average temperature
      const sol = sim.crystalKind(i) || (sim.nb[i] === 1 && sim.isWaterO(sim.bp[i * P.MAXB]))
        ? (sim.isSolidAtom(i) ? 1 : 0) : (look.solid ? 1 : 0);
      const wm = w * wc, wg = w * (1 - wc) * look.gas;
      if (wm > 0) {
        anyCond = true;
        for (let k = 0; k < 4; k++) {
          const b = k & 1 ? tx : 1 - tx, a = k & 2 ? ty : 1 - ty, f = a * b;
          if (f <= 0) continue;
          const cc = c + (k & 1) + (k & 2 ? gw : 0), m = wm * f;
          dens[cc] += m;
          cr[cc] += look.rgb[0] * m; cg[cc] += look.rgb[1] * m; cb[cc] += look.rgb[2] * m;
          cop[cc] += look.opacity * m; csol[cc] += sol * m;
          heat[cc] += tLoc[i] * m;
        }
      }
      if (wg > 0) {
        // (on the coarse grid each cell holds GF^2 fine cells' worth: per fine cell, as above)
        anyGas = true;
        let qx = (x[i] - wx0) / cg4 - 0.5, qy = (y[i] - wy0) / cg4 - 0.5;
        qx = qx < 0 ? 0 : qx > ggw - 1.001 ? ggw - 1.001 : qx;
        qy = qy < 0 ? 0 : qy > ggh - 1.001 ? ggh - 1.001 : qy;
        const hx = Math.floor(qx), hy = Math.floor(qy), ux = qx - hx, uy = qy - hy, c4 = hy * ggw + hx, g0 = wg / (GF * GF);
        for (let k = 0; k < 4; k++) {
          const f = (k & 2 ? uy : 1 - uy) * (k & 1 ? ux : 1 - ux);
          if (f <= 0) continue;
          const cc = c4 + (k & 1) + (k & 2 ? ggw : 0), g = g0 * f;
          ga[cc] += g; gr[cc] += look.rgb[0] * g; gg[cc] += look.rgb[1] * g; gb[cc] += look.rgb[2] * g;
        }
      }
    }
    // blur about one atom diameter, so packed atoms (metals sit ~36 units apart) read as one piece
    // (a picture with no liquid, solid or gas in it skips that blur: it is all zero)
    if (anyCond) for (let k = 0; k < 7; k++) { boxBlur(ch[k], gw, gh, R0, rg.tmp); boxBlur(ch[k], gw, gh, R0, rg.tmp); }
    if (anyGas) for (let k = 0; k < 4; k++) { boxBlur(gch[k], ggw, ggh, RG, rg.gtmp); boxBlur(gch[k], ggw, ggh, RG, rg.gtmp); }
    // ease toward this frame's picture (time constant ~0.12 s)
    const prev = rg.prev, gprev = rg.gprev;
    const ease = rg.fresh ? 1 : 1 - Math.exp(-(dt || 0.016) / 0.12);
    rg.fresh = false;
    for (let k = 0; k < 7; k++) {
      const a = ch[k], b = prev[k];
      for (let c = 0; c < N; c++) { b[c] += ease * (a[c] - b[c]); a[c] = b[c]; }
    }
    for (let k = 0; k < 4; k++) {
      const a = gch[k], b = gprev[k];
      for (let c = 0; c < GN; c++) { b[c] += ease * (a[c] - b[c]); a[c] = b[c]; }
    }
    // tinted gas
    const gimg = this.gasImg.data;
    for (let c = 0; c < GN; c++) {
      const o = c * 4;
      if (ga[c] > 1e-4) {
        gimg[o] = gr[c] / ga[c]; gimg[o + 1] = gg[c] / ga[c]; gimg[o + 2] = gb[c] / ga[c];
        gimg[o + 3] = Math.min(0.7, ga[c] * 1.6) * 255;
      } else gimg[o + 3] = 0;
    }
    const img = this.realImg.data;
    for (let gy = 0; gy < gh; gy++) {
      for (let gx = 0; gx < gw; gx++) {
        const c = gy * gw + gx, o = c * 4;
        let R = 0, G = 0, B = 0, A = 0;
        // liquids and solids (drawn over the gas)
        const d = dens[c];
        const cov = smooth(0.3, 0.55, d);
        if (cov > 0) {
          let r = cr[c] / d, g = cg[c] / d, b = cb[c] / d;
          // solid fraction blends the looks (opaque, matte) and (translucent, glossy) smoothly
          const sf = smooth(0.2, 0.8, csol[c] / d);
          let a = cov * (sf + (1 - sf) * Math.min(1, cop[c] / d));
          // light from above: brighten upward-facing surfaces, darken undersides
          const up = gy > 0 ? dens[c - gw] : 0, down = gy < gh - 1 ? dens[c + gw] : d;
          const edge = smooth(0.15, 0.5, d - up) - 0.6 * smooth(0.15, 0.5, d - down);
          const lift = (0.7 - 0.5 * sf) * edge;
          if (lift > 0) { r += (255 - r) * lift; g += (255 - g) * lift; b += (255 - b) * lift; a = Math.min(1, a + 0.45 * lift); }
          else { r *= 1 + lift * 0.35; g *= 1 + lift * 0.35; b *= 1 + lift * 0.35; }
          // hot matter glows: dull red from about 800 K (the Draper point), orange-white when molten
          const glow = incandescence(heat[c] / d);
          if (glow) {
            const k = glow[3];
            r += (glow[0] - r) * k; g += (glow[1] - g) * k; b += (glow[2] - b) * k;
            a += (1 - a) * k;
          }
          R = r; G = g; B = b; A = a;
        }
        img[o] = R; img[o + 1] = G; img[o + 2] = B; img[o + 3] = A * 255;
      }
    }
    this.realCtx.putImageData(this.realImg, 0, 0);
    if (anyGas) this.gasCtx.putImageData(this.gasImg, 0, 0);
    ctx.save();
    ctx.translate(this.width / 2, this.height / 2);
    ctx.scale(this.cam.scale, this.cam.scale);
    ctx.translate(-this.cam.x, -this.cam.y);
    ctx.globalAlpha = alpha;
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    // (clipped to the jar: the window's last cells can reach past the glass)
    ctx.beginPath(); ctx.rect(0, 0, W, H); ctx.clip();
    if (anyGas) ctx.drawImage(this.gasCanvas, wx0, wy0, ggw * cg4, ggh * cg4);
    ctx.drawImage(this.realCanvas, wx0, wy0, gw * cw, gh * cw);
    // glass reflections
    const glare = ctx.createLinearGradient(0, 0, W, 0);
    glare.addColorStop(0, 'rgba(255,255,255,0.10)');
    glare.addColorStop(0.08, 'rgba(255,255,255,0.02)');
    glare.addColorStop(0.85, 'rgba(255,255,255,0)');
    glare.addColorStop(0.93, 'rgba(255,255,255,0.07)');
    glare.addColorStop(1, 'rgba(255,255,255,0)');
    ctx.fillStyle = glare;
    ctx.fillRect(0, 0, W, H);
    ctx.restore();
  };

  // Light: flashes of chemiluminescence from reactions (drained from the simulation) and the
  // flame-test glow of hot metal atoms in the gas, gathered on a coarse grid that fades with a
  // fraction of a second's afterglow (as the eye sees a flame), blurred and added over the scene.
  // A bright flame also lights up the jar around it.
  const LC = 40, LIGHT_TAU = 0.25;
  R.drawLight = function (ctx, dt) {
    const sim = this.sim, gw = Math.ceil(P.W / LC), gh = Math.ceil(P.H / LC), N = gw * gh;
    let L = this._light;
    if (!L || L.N !== N) {
      L = this._light = { N: N, r: new Float32Array(N), g: new Float32Array(N), b: new Float32Array(N),
        br: new Float32Array(N), bg: new Float32Array(N), bb: new Float32Array(N), tmp: new Float32Array(N), img: null };
      this.lightCanvas.width = gw; this.lightCanvas.height = gh;
      L.img = this.lightCtx.createImageData(gw, gh);
    }
    const fade = Math.exp(-(dt || 0.016) / LIGHT_TAU);
    let any = false;
    for (let c = 0; c < N; c++) { L.r[c] *= fade; L.g[c] *= fade; L.b[c] *= fade; if (L.r[c] + L.g[c] + L.b[c] > 1e-3) any = true; }
    const cell = function (x, y) {
      const cx = Math.min(gw - 1, Math.max(0, Math.floor(x / LC))), cy = Math.min(gh - 1, Math.max(0, Math.floor(y / LC)));
      return cy * gw + cx;
    };
    // flashes from reactions (energy in kJ/mol; ~400 for a new O-H bond)
    const fl = sim.light;
    if (fl && fl.n) {
      for (let k = 0; k < fl.n; k++) {
        const o = 6 * k, c = cell(fl.buf[o], fl.buf[o + 1]), e = fl.buf[o + 2] / 4000;
        L.r[c] += e * fl.buf[o + 3] / 255; L.g[c] += e * fl.buf[o + 4] / 255; L.b[c] += e * fl.buf[o + 5] / 255;
      }
      fl.n = 0; any = true;
    }
    // hot metal atoms in the gas: their flame-test colours, rising from ~1000 K
    const flame = CE.Simulation.flameColor, el = sim.el, tLoc = sim.tLoc, condF = sim.condF, d = dt || 0.016;
    for (let i = 0; i < sim.n; i++) {
      const col = flame[el[i]];
      if (!col || tLoc[i] < 1000 || condF[i] > 0.5) continue;
      const t = Math.min(1, (tLoc[i] - 1000) / 1500), e = 1.2 * t * t * d, c = cell(sim.x[i], sim.y[i]);
      L.r[c] += e * col[0] / 255; L.g[c] += e * col[1] / 255; L.b[c] += e * col[2] / 255; any = true;
    }
    if (!any) return;
    // a soft copy, a few cells wide (a glow, not sparkles)
    L.br.set(L.r); L.bg.set(L.g); L.bb.set(L.b);
    for (const a of [L.br, L.bg, L.bb]) { boxBlur(a, gw, gh, 2, L.tmp); boxBlur(a, gw, gh, 2, L.tmp); }
    // on a dark background light simply adds; on a light one (a lit room, where a flame is hard
    // to see) it is drawn as a coloured haze of the same strength
    const bgc = this.theme.bg, dark = bgc[0] !== '#' || parseInt(bgc.substr(1, 2), 16) < 128;
    const img = L.img.data;
    let sr = 0, sg = 0, sb = 0;
    for (let c = 0; c < N; c++) {
      const r = 1 - Math.exp(-L.br[c]), g = 1 - Math.exp(-L.bg[c]), b = 1 - Math.exp(-L.bb[c]), o = 4 * c;
      if (dark) { img[o] = 255 * r; img[o + 1] = 255 * g; img[o + 2] = 255 * b; img[o + 3] = 255; }
      else {
        const m = Math.max(r, g, b, 1e-6);
        img[o] = 255 * r / m; img[o + 1] = 255 * g / m; img[o + 2] = 255 * b / m; img[o + 3] = 255 * Math.min(1, 2.5 * m);
      }
      sr += L.r[c]; sg += L.g[c]; sb += L.b[c];
    }
    this.lightCtx.putImageData(L.img, 0, 0);
    ctx.save();
    ctx.globalCompositeOperation = dark ? 'lighter' : 'source-over';
    // the flame lights up the whole jar a little
    const tot = sr + sg + sb;
    if (dark && tot > 0.05) {
      const a = Math.min(0.15, tot * 0.003), m = Math.max(sr, sg, sb) || 1;
      ctx.fillStyle = 'rgba(' + Math.round(255 * sr / m) + ',' + Math.round(255 * sg / m) + ',' + Math.round(255 * sb / m) + ',' + a.toFixed(3) + ')';
      ctx.fillRect(0, 0, this.width, this.height);
    }
    ctx.translate(this.width / 2, this.height / 2);
    ctx.scale(this.cam.scale, this.cam.scale);
    ctx.translate(-this.cam.x, -this.cam.y);
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(this.lightCanvas, 0, 0, gw * LC, gh * LC);
    ctx.restore();
  };

  // The grab tool's hand: a ring at the cursor and a line to the centre of what it holds, which
  // stretches when the material resists
  R.drawGrab = function (ctx) {
    const g = this.sim.grab;
    if (!g) return;
    const sim = this.sim, A = g.atoms;
    let cx = 0, cy = 0;
    for (let k = 0; k < A.length; k++) { cx += sim.x[A[k]]; cy += sim.y[A[k]]; }
    cx /= A.length; cy /= A.length;
    const hx = g.hx + g.ox, hy = g.hy + g.oy, px = 1 / this.cam.scale;
    ctx.save();
    ctx.strokeStyle = 'rgba(255, 200, 60, 0.9)';
    ctx.lineWidth = 2 * px;
    ctx.setLineDash([6 * px, 4 * px]);
    ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(hx, hy); ctx.stroke();
    ctx.setLineDash([]);
    ctx.beginPath(); ctx.arc(hx, hy, 9 * px, 0, Math.PI * 2); ctx.stroke();
    ctx.fillStyle = 'rgba(255, 200, 60, 0.9)';
    ctx.beginPath(); ctx.arc(cx, cy, 3.5 * px, 0, Math.PI * 2); ctx.fill();
    ctx.restore();
  };

  // Colour and strength of the thermal glow of condensed matter at temperature T (K): a
  // black-body ramp from dark red through orange to yellow-white, or null when it is too cool
  // to be seen glowing in a lit room.
  const GLOW = [[800, 110, 10, 0], [1000, 190, 35, 0], [1300, 255, 95, 10], [1700, 255, 165, 55],
    [2200, 255, 210, 130], [3200, 255, 240, 215]];
  function incandescence(T) {
    if (!(T > GLOW[0][0])) return null;
    let k = 1;
    while (k < GLOW.length - 1 && T > GLOW[k][0]) k++;
    const a = GLOW[k - 1], b = GLOW[k], f = Math.min(1, (T - a[0]) / (b[0] - a[0]));
    return [a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f, a[3] + (b[3] - a[3]) * f, smooth(800, 1500, T) * 0.9];
  }

  function smooth(a, b, v) {
    const t = Math.min(1, Math.max(0, (v - a) / (b - a)));
    return t * t * (3 - 2 * t);
  }

  // Copies grid src (sw x sh) into dst (dw x dh) moved by (dx, dy) cells: dst cell (i, j) takes
  // src cell (i + dx, j + dy), or 0 where that is outside src.
  function remapGrid(src, sw, sh, dst, dw, dh, dx, dy) {
    dst.fill(0);
    const i0 = Math.max(0, -dx), i1 = Math.min(dw, sw - dx);
    if (i1 <= i0) return;
    for (let j = Math.max(0, -dy), j1 = Math.min(dh, sh - dy); j < j1; j++) {
      const s = (j + dy) * sw + dx;
      dst.set(src.subarray(s + i0, s + i1), j * dw + i0);
    }
  }

  // In-place separable box blur (running sums), radius r cells. The vertical pass keeps one
  // running sum per column and goes row by row, so it reads and writes memory in order.
  let colSum = new Float32Array(0);
  function boxBlur(a, w, h, r, tmp) {
    const inv = 1 / (2 * r + 1);
    for (let yy = 0; yy < h; yy++) {
      const row = yy * w;
      let sum = 0;
      for (let k = -r; k <= r; k++) sum += a[row + Math.min(w - 1, Math.max(0, k))];
      for (let xx = 0; xx < w; xx++) {
        tmp[row + xx] = sum * inv;
        sum += a[row + Math.min(w - 1, xx + r + 1)] - a[row + Math.max(0, xx - r)];
      }
    }
    if (colSum.length < w) colSum = new Float32Array(w);
    const cs = colSum;
    cs.fill(0, 0, w);
    for (let k = -r; k <= r; k++) {
      const row = Math.min(h - 1, Math.max(0, k)) * w;
      for (let xx = 0; xx < w; xx++) cs[xx] += tmp[row + xx];
    }
    for (let yy = 0; yy < h; yy++) {
      const out = yy * w, add = Math.min(h - 1, yy + r + 1) * w, sub = Math.max(0, yy - r) * w;
      for (let xx = 0; xx < w; xx++) {
        a[out + xx] = cs[xx] * inv;
        cs[xx] += tmp[add + xx] - tmp[sub + xx];
      }
    }
  }

  // Zoomed in: ball-and-stick atoms with bonds.
  // Real-world view: the picked species is marked with a yellow dot per atom on top of the
  // material, so you can see where it is (e.g. which part of a pool is water).
  R.drawFocusDots = function (ctx, alpha) {
    const sim = this.sim, n = sim.n, focus = this.focus, s = this.cam.scale;
    const ox = this.width / 2 - this.cam.x * s, oy = this.height / 2 - this.cam.y * s;
    const d = 2.5;
    ctx.save();
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.globalAlpha = alpha;
    ctx.fillStyle = '#ffeb3b';
    ctx.beginPath();
    const vw = this.width + d, vh = this.height + d;
    for (let i = 0; i < n; i++) {
      if (!focus[i]) continue;
      const px = sim.x[i] * s + ox, py = sim.y[i] * s + oy;
      if (px < -d || py < -d || px > vw || py > vh) continue; // off screen
      ctx.rect(px - d / 2, py - d / 2, d, d);
    }
    ctx.fill();
    ctx.restore();
  };

  R.drawAtoms = function (ctx) {
    const sim = this.sim, n = sim.n, MB = P.MAXB;
    const s = this.cam.scale;
    const ox = this.width / 2 - this.cam.x * s, oy = this.height / 2 - this.cam.y * s;
    const x = sim.x, y = sim.y, el = sim.el;
    // Only atoms on screen are drawn: 2 = on screen (within the biggest atom's radius, with its
    // label and highlight ring, of the edge), 1 = off screen but near enough for a bond to it to
    // cross the screen (bonds longer than 40 units are not drawn), 0 = neither. (A fixed margin of
    // 60 px let big atoms and long bonds pop in and out at the edge when zoomed far in.)
    const margin = Math.max(60, 26 * 0.42 * s * 1.6), reach = margin + 40 * s;
    const visible = this._visible && this._visible.length >= n ? this._visible : (this._visible = new Uint8Array(P.MAX));
    const sx = this._sx && this._sx.length >= n ? this._sx : (this._sx = new Float32Array(P.MAX));
    const sy = this._sy && this._sy.length >= n ? this._sy : (this._sy = new Float32Array(P.MAX));
    const vw = this.width, vh = this.height;
    for (let i = 0; i < n; i++) {
      const a = x[i] * s + ox, b = y[i] * s + oy;
      sx[i] = a; sy[i] = b;
      visible[i] = a > -margin && b > -margin && a < vw + margin && b < vh + margin ? 2
        : a > -reach && b > -reach && a < vw + reach && b < vh + reach ? 1 : 0;
    }
    ctx.save();
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    const focus = this.focus;
    // bonds
    ctx.lineCap = 'round';
    const bondW = Math.max(1, 1.6 * s);
    for (let i = 0; i < n; i++) {
      if (!visible[i]) continue;
      for (let k = 0; k < sim.nb[i]; k++) {
        const j = sim.bp[i * MB + k];
        if (j < i && visible[j]) continue; // drawn from j
        if (visible[i] < 2 && visible[j] < 2 && !crossesView(sx[i], sy[i], sx[j], sy[j], vw, vh)) continue;
        const order = sim.bo[i * MB + k], ionic = sim.bc[i * MB + k] === 2;
        const dx = sx[j] - sx[i], dy = sy[j] - sy[i];
        const len = Math.sqrt(dx * dx + dy * dy) || 1;
        if (len > 40 * s) continue; // wrapped / stretched bond
        const nx = -dy / len, ny = dx / len;
        ctx.globalAlpha = focus && !focus[i] ? 0.15 : 1;
        ctx.strokeStyle = ionic ? 'rgba(255,214,102,0.75)' : 'rgba(200,210,225,0.9)';
        ctx.lineWidth = order > 1 ? bondW * 0.6 : bondW;
        if (ionic) ctx.setLineDash([bondW * 1.5, bondW * 1.5]);
        const gap = Math.max(3, bondW * 1.7);
        for (let o = 0; o < order; o++) {
          const off = (o - (order - 1) / 2) * gap;
          ctx.beginPath();
          ctx.moveTo(sx[i] + nx * off, sy[i] + ny * off);
          ctx.lineTo(sx[j] + nx * off, sy[j] + ny * off);
          ctx.stroke();
          if (ionic) break;
        }
        if (ionic) ctx.setLineDash([]);
      }
    }
    // atoms
    const hl = this.highlight;
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    for (let i = 0; i < n; i++) {
      if (visible[i] < 2) continue;
      const z = el[i];
      const r = E[z].radius * 0.42 * s;
      ctx.globalAlpha = focus && !focus[i] ? 0.15 : 1;
      ctx.fillStyle = this.colorMode === 'element' || r > 6 ? this.elColorCss[z] : cssFromPacked(sim.molColor[i]);
      ctx.beginPath(); ctx.arc(sx[i], sy[i], r, 0, Math.PI * 2); ctx.fill();
      if (r > 5) {
        ctx.strokeStyle = 'rgba(0,0,0,0.45)'; ctx.lineWidth = Math.max(0.5, r * 0.08); ctx.stroke();
      }
      if ((hl && hl.has(i)) || (focus && focus[i])) {
        ctx.strokeStyle = '#ffeb3b'; ctx.lineWidth = Math.max(1.5, r * 0.18);
        ctx.beginPath(); ctx.arc(sx[i], sy[i], r + Math.max(2, r * 0.2), 0, Math.PI * 2); ctx.stroke();
      }
      if (r > 8) {
        ctx.fillStyle = this.elTextColor[z];
        ctx.font = '600 ' + Math.round(r * 0.9) + 'px system-ui, sans-serif';
        ctx.fillText(E[z].symbol, sx[i], sy[i] + r * 0.04);
      }
      if (this.showCharges && r > 13) {
        const q = sim.q[i];
        if (Math.abs(q) > 0.15) {
          const ionic = Math.abs(q) >= 0.95 && hasIonic(sim, i);
          const label = ionic ? (Math.abs(Math.round(q)) > 1 ? Math.abs(Math.round(q)) : '') + (q > 0 ? '+' : '−') : (q > 0 ? 'δ+' : 'δ−');
          ctx.font = '600 ' + Math.round(r * 0.55) + 'px system-ui, sans-serif';
          ctx.fillStyle = q > 0 ? '#ff8a80' : '#82b1ff';
          ctx.fillText(label, sx[i] + r * 0.95, sy[i] - r * 0.95);
        }
      }
    }
    ctx.restore();
  };

  // Whether the segment (ax, ay)-(bx, by) passes over the w x h screen (its bounding box does:
  // close enough for a bond, a short segment)
  function crossesView(ax, ay, bx, by, w, h) {
    return Math.max(ax, bx) > 0 && Math.min(ax, bx) < w && Math.max(ay, by) > 0 && Math.min(ay, by) < h;
  }

  function hasIonic(sim, i) {
    for (let k = 0; k < sim.nb[i]; k++) if (sim.bc[i * P.MAXB + k] === 2) return true;
    return false;
  }

  const cssCache = new Map();
  function cssFromPacked(c) {
    let s = cssCache.get(c);
    if (!s) {
      s = 'rgb(' + (c & 255) + ',' + ((c >> 8) & 255) + ',' + ((c >> 16) & 255) + ')';
      cssCache.set(c, s);
    }
    return s;
  }

  R.drawScaleBar = function (ctx) {
    if (this.realMix() >= 0.5) {
      const pxPerCm = this.cam.scale * REAL_UNITS_PER_CM;
      const optsCm = [0.5, 1, 2, 5, 10, 20, 50];
      let Lc = optsCm[0];
      for (let k = 0; k < optsCm.length; k++) if (optsCm[k] * pxPerCm <= 140) Lc = optsCm[k];
      return this.drawBar(ctx, Lc * pxPerCm, Lc + ' cm');
    }
    // choose a round length in nm (1 world unit = 10 pm = 0.01 nm)
    const pxPerNm = this.cam.scale * 100;
    const opts = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100];
    let L = opts[0];
    for (let k = 0; k < opts.length; k++) if (opts[k] * pxPerNm <= 140) L = opts[k];
    this.drawBar(ctx, L * pxPerNm, L < 1 ? (L * 10) + ' Å' : L + ' nm');
  };

  R.drawBar = function (ctx, w, label) {
    const x0 = 16, y0 = this.height - 22;
    ctx.strokeStyle = this.theme.text; ctx.fillStyle = this.theme.text; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x0 + w, y0);
    ctx.moveTo(x0, y0 - 5); ctx.lineTo(x0, y0 + 5); ctx.moveTo(x0 + w, y0 - 5); ctx.lineTo(x0 + w, y0 + 5);
    ctx.stroke();
    ctx.font = '12px system-ui, sans-serif'; ctx.textAlign = 'left'; ctx.textBaseline = 'bottom';
    ctx.fillText(label, x0, y0 - 6);
  };

  function roundRect(ctx, x, y, w, h, r) {
    ctx.beginPath();
    ctx.moveTo(x + r, y);
    ctx.arcTo(x + w, y, x + w, y + h, r);
    ctx.arcTo(x + w, y + h, x, y + h, r);
    ctx.arcTo(x, y + h, x, y, r);
    ctx.arcTo(x, y, x + w, y, r);
    ctx.closePath();
  }

  CE.Renderer = Renderer;
})(window.CE = window.CE || {});
