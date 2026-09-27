// SVG charts for samaya Studio. Time runs on a calendar: the page passes a time axis (see
// calendar() in app.js) with a label per period, the weekends and the month boundaries. One
// colour per series name, the same in every chart of a case; colours are muted, as in
// industrial displays, so that status colours (green, amber, red) keep their meaning.
(function () {
  const PALETTE = ["#2F5F98", "#C27A2C", "#3F8A5A", "#7A5CA6", "#B0463D", "#2E8A8A", "#8C7250",
                   "#A2507A", "#5C6B7C", "#6C9BD2"];
  const colors = new Map();
  function color(name) {
    if (!colors.has(name)) colors.set(name, PALETTE[colors.size % PALETTE.length]);
    return colors.get(name);
  }
  function resetColors() { colors.clear(); }
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
  const fmt = (v) => Math.abs(v) >= 100 ? v.toLocaleString("en-IN", { maximumFractionDigits: 0 })
                                        : v.toLocaleString("en-IN", { maximumFractionDigits: 1 });
  const svg = (w, h, body) => `<svg viewBox="0 0 ${w} ${h}" class="chart" role="img">${body.join("")}</svg>`;
  // A readable nice maximum for a value axis.
  function niceMax(v) {
    if (!(v > 0)) return 1;
    const p = Math.pow(10, Math.floor(Math.log10(v))), m = v / p;
    return (m <= 1 ? 1 : m <= 2 ? 2 : m <= 2.5 ? 2.5 : m <= 5 ? 5 : 10) * p;
  }

  // The time axis under a plot from x0 to x1: weekend shading from y0 to y1, a label per period
  // (thinned when periods are narrow), and the month (or day, for hours) where it changes.
  function timeAxis(time, x0, x1, y0, y1) {
    const out = [], n = time.n, cw = (x1 - x0) / n;
    for (let t = 0; t < n; t++) {
      if (time.weekend(t)) out.push(`<rect x="${x0 + t * cw}" y="${y0}" width="${cw}" height="${y1 - y0}" class="wk"/>`);
    }
    const every = Math.max(1, Math.ceil(time.minWidth / cw));
    for (let t = 0; t < n; t++) {
      const x = x0 + (t + 0.5) * cw;
      const [a, b] = time.tick(t);
      if (t % every === 0) {
        out.push(`<text x="${x}" y="${y1 + 15}" class="ax" text-anchor="middle">${esc(a)}</text>`);
        if (b) out.push(`<text x="${x}" y="${y1 + 28}" class="ax2" text-anchor="middle">${esc(b)}</text>`);
      }
      const g = time.group(t);
      if (g) {
        out.push(`<line x1="${x0 + t * cw}" y1="${y0}" x2="${x0 + t * cw}" y2="${y1 + 44}" class="mline"/>`);
        out.push(`<text x="${x0 + t * cw + 4}" y="${y1 + 44}" class="axm">${esc(g)}</text>`);
      }
    }
    return out;
  }
  function valueAxis(x0, x1, y0, y1, vmax, unit) {
    const out = [];
    for (let k = 0; k <= 4; k++) {
      const y = y1 - (y1 - y0) * k / 4;
      out.push(`<line x1="${x0}" y1="${y}" x2="${x1}" y2="${y}" class="grid"/>`);
      out.push(`<text x="${x0 - 8}" y="${y + 4}" class="ax" text-anchor="end">${fmt(vmax * k / 4)}</text>`);
    }
    if (unit) out.push(`<text x="${x0 - 8}" y="${y0 - 12}" class="ax" text-anchor="end">${esc(unit)}</text>`);
    return out;
  }
  function legend(names, label, x = 64, y = 8) {
    const out = [];
    for (const nm of names) {
      out.push(`<rect x="${x}" y="${y}" width="11" height="11" fill="${color(nm)}"/>`);
      out.push(`<text x="${x + 17}" y="${y + 10}" class="lg">${esc(label(nm))}</text>`);
      x += 34 + 6.6 * label(nm).length;
    }
    return out;
  }

  // Stacked columns per period (for example CDU charge by crude, per day).
  function stacked(data, time, unit, label, w = 960, h = 330) {
    const x0 = 64, x1 = w - 12, y0 = 44, y1 = h - 52;
    const names = Object.keys(data).sort();
    const tot = [];
    for (let t = 0; t < time.n; t++) tot.push(names.reduce((s, k) => s + Math.max(0, data[k][t] || 0), 0));
    const vmax = niceMax(Math.max(1e-9, ...tot) * 1.05);
    const body = timeAxis(time, x0, x1, y0, y1).concat(valueAxis(x0, x1, y0, y1, vmax, unit), legend(names, label));
    const cw = (x1 - x0) / time.n, bw = Math.max(2, cw * 0.62);
    for (let t = 0; t < time.n; t++) {
      const x = x0 + (t + 0.5) * cw - bw / 2;
      let base = 0;
      for (const k of names) {
        const v = Math.max(0, data[k][t] || 0);
        if (v <= 0) continue;
        const yt = y1 - (y1 - y0) * (base + v) / vmax, hh = (y1 - y0) * v / vmax;
        body.push(`<rect x="${x}" y="${yt}" width="${bw}" height="${hh}" fill="${color(k)}"><title>${esc(label(k))}, ${esc(time.title(t))}: ${fmt(v)} ${esc(unit)}</title></rect>`);
        base += v;
      }
    }
    body.push(`<line x1="${x0}" y1="${y1}" x2="${x1}" y2="${y1}" class="base"/>`);
    return svg(w, h, body);
  }

  // Lines per series over time (for example tank stock per crude, end of day).
  function lines(data, time, unit, label, w = 960, h = 330) {
    const x0 = 64, x1 = w - 12, y0 = 44, y1 = h - 52;
    const names = Object.keys(data).sort();
    const vmax = niceMax(Math.max(1e-9, ...names.map((k) => Math.max(...data[k]))) * 1.05);
    const body = timeAxis(time, x0, x1, y0, y1).concat(valueAxis(x0, x1, y0, y1, vmax, unit), legend(names, label));
    const cw = (x1 - x0) / time.n;
    for (const k of names) {
      const pts = data[k].map((v, t) => `${(x0 + (t + 0.5) * cw).toFixed(1)},${(y1 - (y1 - y0) * Math.max(0, v) / vmax).toFixed(1)}`);
      body.push(`<polyline points="${pts.join(" ")}" fill="none" stroke="${color(k)}" stroke-width="2" stroke-linejoin="round"><title>${esc(label(k))}</title></polyline>`);
    }
    body.push(`<line x1="${x0}" y1="${y1}" x2="${x1}" y2="${y1}" class="base"/>`);
    return svg(w, h, body);
  }

  // Horizontal bars with names (for example marginal values of capacity).
  function hbars(items, label, unit, w = 960) {
    const h = 24 + 28 * items.length;
    const vmax = Math.max(1e-9, ...items.map((i) => Math.abs(i[1]))) * 1.08;
    const body = [];
    items.forEach(([k, v], i) => {
      const y = 10 + i * 28;
      const bw = (w - 360) * Math.abs(v) / vmax;
      body.push(`<text x="280" y="${y + 14}" class="ax" text-anchor="end">${esc(label(k))}</text>`);
      body.push(`<rect x="290" y="${y + 2}" width="${bw}" height="16" fill="${v >= 0 ? "var(--bar)" : "var(--bar-neg)"}"><title>${esc(label(k))}: ${fmt(v)} ${esc(unit || "")}</title></rect>`);
      body.push(`<text x="${296 + bw}" y="${y + 15}" class="ax">${fmt(v)}</text>`);
    });
    return svg(w, h, body);
  }

  // The berth schedule as a Gantt chart: a lane per kind of cargo (term, spot), a bar on the day
  // each cargo berths, coloured by crude and labelled with its cargo number.
  function gantt(arrivals, time, label, w = 960) {
    const lanes = [...new Set(arrivals.map((a) => a.lane))];
    if (!lanes.length) lanes.push("Term cargoes");
    const x0 = 150, x1 = w - 12, y0 = 38, lh = 40, y1 = y0 + lh * lanes.length;
    const h = y1 + 54;
    const crudes = [...new Set(arrivals.map((a) => a.crude))].sort();
    const body = legend(crudes, label, x0, 8).concat(timeAxis(time, x0, x1, y0, y1));
    const cw = (x1 - x0) / time.n;
    lanes.forEach((ln, i) => {
      const y = y0 + i * lh;
      body.push(`<line x1="${x0}" y1="${y + lh}" x2="${x1}" y2="${y + lh}" class="grid"/>`);
      body.push(`<text x="${x0 - 10}" y="${y + lh / 2 + 4}" class="ax" text-anchor="end">${esc(ln)}</text>`);
    });
    for (let t = 0; t <= time.n; t++) body.push(`<line x1="${x0 + t * cw}" y1="${y0}" x2="${x0 + t * cw}" y2="${y1}" class="vgrid"/>`);
    for (const a of arrivals) {
      const y = y0 + lanes.indexOf(a.lane) * lh;
      const x = x0 + a.t * cw;
      body.push(`<rect x="${x + 2}" y="${y + 7}" width="${cw - 4}" height="${lh - 14}" fill="${color(a.crude)}"><title>${esc(time.title(a.t))}: ${esc(a.cargo)}, ${esc(label(a.crude))}</title></rect>`);
      if (cw >= 26) body.push(`<text x="${x + cw / 2}" y="${y + lh / 2 + 4}" class="bar-lbl" text-anchor="middle">${esc(a.short)}</text>`);
    }
    body.push(`<rect x="${x0}" y="${y0}" width="${x1 - x0}" height="${y1 - y0}" class="frame"/>`);
    return svg(w, h, body);
  }

  // On/off by period, one row per unit (for example boilers by hour).
  function onoff(on, time, label, w = 960) {
    const names = Object.keys(on).sort();
    const x0 = 150, x1 = w - 12, y0 = 10, lh = 30, y1 = y0 + lh * names.length;
    const h = y1 + 54;
    const cw = (x1 - x0) / time.n;
    const body = timeAxis(time, x0, x1, y0, y1);
    names.forEach((k, i) => {
      const y = y0 + i * lh;
      body.push(`<text x="${x0 - 10}" y="${y + 19}" class="ax" text-anchor="end">${esc(label(k))}</text>`);
      for (let t = 0; t < time.n; t++) {
        const v = on[k][t] || 0;
        body.push(`<rect x="${x0 + t * cw + 1}" y="${y + 5}" width="${Math.max(1, cw - 2)}" height="${lh - 10}" class="${v > 0.5 ? "on" : "off"}"><title>${esc(label(k))}, ${esc(time.title(t))}: ${v > 0.5 ? "on" : "off"}</title></rect>`);
      }
    });
    return svg(w, h, body);
  }

  window.Charts = { stacked, lines, hbars, gantt, onoff, color, resetColors, esc };
})();
