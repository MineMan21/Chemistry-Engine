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
    this.showCharges = true;
    this.pixelCanvas = document.createElement('canvas');
    this.pixelCtx = this.pixelCanvas.getContext('2d');
    this.imageData = null;
    this.time = 0;
    this.sparks = [];
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
      if (mix > 0) this.drawReal(ctx, mix);
    }

    ctx.save();
    ctx.translate(this.width / 2, this.height / 2);
    ctx.scale(this.cam.scale, this.cam.scale);
    ctx.translate(-this.cam.x, -this.cam.y);
    this.drawApparatusFront(ctx);
    this.drawSparks(ctx, dt);
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
  R.drawPixels = function (ctx, alpha) {
    const sim = this.sim, n = sim.n;
    const w = Math.ceil(this.width * this.dpr), h = Math.ceil(this.height * this.dpr);
    if (this.pixelCanvas.width !== w || this.pixelCanvas.height !== h) {
      this.pixelCanvas.width = w; this.pixelCanvas.height = h;
      this.imageData = null;
    }
    if (!this.imageData) {
      this.imageData = this.pixelCtx.createImageData(w, h);
      this.pix32 = new Uint32Array(this.imageData.data.buffer);
    }
    const pix = this.pix32;
    pix.fill(0);
    const s = this.cam.scale * this.dpr;
    const ox = this.width / 2 * this.dpr - this.cam.x * s, oy = this.height / 2 * this.dpr - this.cam.y * s;
    const colors = this.colorMode === 'element' ? null : sim.molColor;
    const elColor = CE.Simulation.elColor;
    const x = sim.x, y = sim.y, el = sim.el;
    const hl = this.highlight;
    for (let i = 0; i < n; i++) {
      const px = Math.floor(x[i] * s + ox), py = Math.floor(y[i] * s + oy);
      const size = Math.max(1, Math.round(E[el[i]].radius * 1.2 * s));
      if (px < -size || py < -size || px >= w + size || py >= h + size) continue;
      let col = colors ? colors[i] : elColor[el[i]];
      if (hl && hl.has(i)) col = 0xff00ffff; // yellow
      const x0 = Math.max(0, px - (size >> 1)), y0 = Math.max(0, py - (size >> 1));
      const x1 = Math.min(w, x0 + size), y1 = Math.min(h, y0 + size);
      for (let yy = y0; yy < y1; yy++) {
        const row = yy * w;
        for (let xx = x0; xx < x1; xx++) pix[row + xx] = col;
      }
    }
    this.pixelCtx.putImageData(this.imageData, 0, 0);
    ctx.save();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.imageSmoothingEnabled = false;
    ctx.globalAlpha = alpha === undefined ? 1 : alpha;
    ctx.drawImage(this.pixelCanvas, 0, 0);
    ctx.restore();
  };

  // Real-world view: atoms are splatted onto a coarse grid, blurred into continuous material and
  // shaded by what the substance looks like (translucent liquid with a bright surface, opaque
  // solid, tinted gas). The simulation underneath is unchanged.
  R.drawReal = function (ctx, alpha) {
    const sim = this.sim, n = sim.n, W = P.W, H = P.H;
    const cw = Math.max(6, 2.6 / this.cam.scale); // grid cell in world units (~2.6 screen px)
    const gw = Math.ceil(W / cw), gh = Math.ceil(H / cw), N = gw * gh;
    if (!this._rg || this._rg.N !== N) {
      const ch = [];
      for (let k = 0; k < 11; k++) ch.push(new Float32Array(N));
      this._rg = { N: N, ch: ch, tmp: new Float32Array(N) };
      this.realCanvas = document.createElement('canvas');
      this.realCanvas.width = gw; this.realCanvas.height = gh;
      this.realCtx = this.realCanvas.getContext('2d');
      this.realImg = this.realCtx.createImageData(gw, gh);
    }
    // channels: 0 density, 1-3 colour, 4 opacity, 5 solid, 6 gas strength, 7-9 gas colour, 10 heat
    const ch = this._rg.ch;
    for (let k = 0; k < 11; k++) ch[k].fill(0);
    const [dens, cr, cg, cb, cop, csol, ga, gr, gg, gb] = ch;
    const species = sim.species || [], spIdx = sim.spIdx, x = sim.x, y = sim.y, el = sim.el;
    const area = 1 / (cw * cw);
    for (let i = 0; i < n; i++) {
      const sp = spIdx ? species[spIdx[i]] : null;
      const look = sp && sp.look;
      if (!look) continue;
      let gx = Math.floor(x[i] / cw), gy = Math.floor(y[i] / cw);
      if (gx < 0) gx = 0; else if (gx >= gw) gx = gw - 1;
      if (gy < 0) gy = 0; else if (gy >= gh) gy = gh - 1;
      const c = gy * gw + gx;
      const r = E[el[i]].radius;
      const w = Math.PI * r * r * area;
      // drawn as liquid/solid only when genuinely packed (several neighbours, persistently), so a
      // few gas molecules that happen to touch stay invisible vapour
      if ((sim.condF[i] > 0.6 && sim.contacts[i] >= 2) || (sim.nb[i] >= 2 && sim.contacts[i] >= 1) || sim.condF[i] > 0.9) {
        dens[c] += w;
        cr[c] += look.rgb[0] * w; cg[c] += look.rgb[1] * w; cb[c] += look.rgb[2] * w;
        cop[c] += look.opacity * w; csol[c] += (look.solid ? 1 : 0) * w;
      } else if (look.gas > 0) {
        const s = look.gas * w;
        ga[c] += s; gr[c] += look.rgb[0] * s; gg[c] += look.rgb[1] * s; gb[c] += look.rgb[2] * s;
      }
    }
    const R0 = Math.max(1, Math.round(20 / cw)), RG = Math.max(3, Math.round(160 / cw));
    for (let k = 0; k <= 5; k++) { boxBlur(ch[k], gw, gh, R0, this._rg.tmp); boxBlur(ch[k], gw, gh, R0, this._rg.tmp); }
    for (let k = 6; k <= 9; k++) { boxBlur(ch[k], gw, gh, RG, this._rg.tmp); boxBlur(ch[k], gw, gh, RG, this._rg.tmp); }
    const img = this.realImg.data;
    for (let gy = 0; gy < gh; gy++) {
      for (let gx = 0; gx < gw; gx++) {
        const c = gy * gw + gx, o = c * 4;
        // tinted gas
        let R = 0, G = 0, B = 0, A = 0;
        if (ga[c] > 1e-4) {
          A = Math.min(0.7, ga[c] * 1.6);
          R = gr[c] / ga[c]; G = gg[c] / ga[c]; B = gb[c] / ga[c];
        }
        // liquids and solids
        const d = dens[c];
        const cov = smooth(0.3, 0.55, d);
        if (cov > 0) {
          let r = cr[c] / d, g = cg[c] / d, b = cb[c] / d;
          const solid = csol[c] / d > 0.5;
          let a = solid ? cov : cov * Math.min(1, cop[c] / d);
          // light from above: brighten upward-facing surfaces, darken undersides
          const up = gy > 0 ? dens[c - gw] : 0, down = gy < gh - 1 ? dens[c + gw] : d;
          const edge = smooth(0.15, 0.5, d - up) - 0.6 * smooth(0.15, 0.5, d - down);
          const lift = solid ? 0.2 * edge : 0.7 * edge;
          if (lift > 0) { r += (255 - r) * lift; g += (255 - g) * lift; b += (255 - b) * lift; a = Math.min(1, a + 0.45 * lift); }
          else { r *= 1 + lift * 0.35; g *= 1 + lift * 0.35; b *= 1 + lift * 0.35; }
          // composite the material over the gas
          const outA = a + A * (1 - a);
          if (outA > 0) {
            R = (r * a + R * A * (1 - a)) / outA; G = (g * a + G * A * (1 - a)) / outA; B = (b * a + B * A * (1 - a)) / outA;
          }
          A = outA;
        }
        img[o] = R; img[o + 1] = G; img[o + 2] = B; img[o + 3] = A * 255;
      }
    }
    this.realCtx.putImageData(this.realImg, 0, 0);
    ctx.save();
    ctx.translate(this.width / 2, this.height / 2);
    ctx.scale(this.cam.scale, this.cam.scale);
    ctx.translate(-this.cam.x, -this.cam.y);
    ctx.globalAlpha = alpha;
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = 'high';
    ctx.drawImage(this.realCanvas, 0, 0, gw * cw, gh * cw);
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

  function smooth(a, b, v) {
    const t = Math.min(1, Math.max(0, (v - a) / (b - a)));
    return t * t * (3 - 2 * t);
  }

  // In-place separable box blur (running sums), radius r cells.
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
    for (let xx = 0; xx < w; xx++) {
      let sum = 0;
      for (let k = -r; k <= r; k++) sum += tmp[Math.min(h - 1, Math.max(0, k)) * w + xx];
      for (let yy = 0; yy < h; yy++) {
        a[yy * w + xx] = sum * inv;
        sum += tmp[Math.min(h - 1, yy + r + 1) * w + xx] - tmp[Math.max(0, yy - r) * w + xx];
      }
    }
  }

  // Zoomed in: ball-and-stick atoms with bonds.
  R.drawAtoms = function (ctx) {
    const sim = this.sim, n = sim.n, MB = P.MAXB;
    const s = this.cam.scale;
    const ox = this.width / 2 - this.cam.x * s, oy = this.height / 2 - this.cam.y * s;
    const x = sim.x, y = sim.y, el = sim.el;
    const margin = 60;
    const visible = this._visible && this._visible.length >= n ? this._visible : (this._visible = new Uint8Array(P.MAX));
    const sx = this._sx && this._sx.length >= n ? this._sx : (this._sx = new Float32Array(P.MAX));
    const sy = this._sy && this._sy.length >= n ? this._sy : (this._sy = new Float32Array(P.MAX));
    for (let i = 0; i < n; i++) {
      sx[i] = x[i] * s + ox; sy[i] = y[i] * s + oy;
      visible[i] = (sx[i] > -margin && sy[i] > -margin && sx[i] < this.width + margin && sy[i] < this.height + margin) ? 1 : 0;
    }
    ctx.save();
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    // bonds
    ctx.lineCap = 'round';
    const bondW = Math.max(1, 1.6 * s);
    for (let i = 0; i < n; i++) {
      if (!visible[i]) continue;
      for (let k = 0; k < sim.nb[i]; k++) {
        const j = sim.bp[i * MB + k];
        if (j < i && visible[j]) continue;
        const order = sim.bo[i * MB + k], ionic = sim.bc[i * MB + k] === 2;
        const dx = sx[j] - sx[i], dy = sy[j] - sy[i];
        const len = Math.sqrt(dx * dx + dy * dy) || 1;
        if (len > 40 * s) continue; // wrapped / stretched bond
        const nx = -dy / len, ny = dx / len;
        ctx.strokeStyle = ionic ? 'rgba(255,214,102,0.75)' : 'rgba(200,210,225,0.9)';
        ctx.lineWidth = order > 1 ? bondW * 0.7 : bondW;
        if (ionic) ctx.setLineDash([bondW * 1.5, bondW * 1.5]);
        const gap = bondW * 1.1;
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
      if (!visible[i]) continue;
      const z = el[i];
      const r = E[z].radius * 0.42 * s;
      ctx.fillStyle = this.colorMode === 'element' || r > 6 ? this.elColorCss[z] : cssFromPacked(sim.molColor[i]);
      ctx.beginPath(); ctx.arc(sx[i], sy[i], r, 0, Math.PI * 2); ctx.fill();
      if (r > 5) {
        ctx.strokeStyle = 'rgba(0,0,0,0.45)'; ctx.lineWidth = Math.max(0.5, r * 0.08); ctx.stroke();
      }
      if (hl && hl.has(i)) {
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
