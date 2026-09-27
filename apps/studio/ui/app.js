// samaya Studio: the page. Talks to the Windows host (apps/studio/main.cpp) through WebView2
// messages; outside the app (a plain browser) it runs in a mock mode on a recorded solve.
//
// Each solve is a case: a model file, its settings (including the start date that puts the
// model's periods on the calendar) and its result. Cases of the same model can be compared.
(function () {
  "use strict";
  const $ = (s) => document.querySelector(s);
  const esc = (s) => String(s ?? "").replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);
  const host = window.chrome && window.chrome.webview ? window.chrome.webview : null;
  const TOL = 1e-6;  // the verifier's default tolerances (include/samaya/verify.hpp)

  const state = { cases: [], selected: null, compare: null, rerun: null, nextId: 1, samples: [], gpu: false,
                  tab: "overview", filter: "" };

  // ---------------------------------------------------------------------------------------------
  // Formatting.
  function num(v, sig = 10) {
    if (v === null || v === undefined || typeof v !== "number" || !isFinite(v)) return "–";
    if (v !== 0 && (Math.abs(v) < 1e-3 || Math.abs(v) >= 1e12)) return v.toExponential(2);
    return v.toLocaleString("en-IN", { maximumSignificantDigits: sig });
  }
  function secs(v) { return !isFinite(v) ? "–" : v < 1 ? `${v.toFixed(3)} s` : v < 60 ? `${v.toFixed(2)} s` : `${Math.floor(v / 60)} min ${Math.round(v % 60)} s`; }
  function stem(path) { const n = String(path).split(/[\\/]/).pop(); return n.replace(/\.(mps|qps|txt)$/i, ""); }
  const fileName = (path) => String(path).split(/[\\/]/).pop();
  // Model names as a planner reads them: arab_light -> Arab Light, feed_cdu -> CDU feed parts, etc.
  const ACRONYMS = new Set(["cdu", "vdu", "fcc", "hcu", "dhds", "lpg", "atf", "lco", "vgo", "vr", "sr", "ln", "hn",
                            "gtg", "stg", "mp", "hp", "mw", "fo", "spm"]);
  function pretty(key) {
    return String(key).split("_").filter(Boolean).map((w) => ACRONYMS.has(w) ? w.toUpperCase()
      : /^\d+$/.test(w) ? w : w[0].toUpperCase() + w.slice(1)).join(" ");
  }
  const STATUS = {
    optimal: ["Optimal", "ok"], infeasible: ["Infeasible", "info"], unbounded: ["Unbounded", "info"],
    infeasible_or_unbounded: ["Infeasible or unbounded", "info"], time_limit: ["Time limit", "warn"],
    iteration_limit: ["Iteration limit", "warn"], node_limit: ["Node limit", "warn"], numerical_error: ["Numerical error", "bad"],
    invalid_model: ["Invalid model", "bad"], not_implemented: ["Not supported", "bad"], not_convex: ["Not convex", "bad"],
    not_solved: ["Not solved", "bad"], interrupted: ["Stopped", "warn"],
  };
  const statusOf = (s) => STATUS[s] || [String(s || "unknown").replace(/_/g, " "), "bad"];

  // ---------------------------------------------------------------------------------------------
  // Calendar. Model periods (day, week or hour) become dates from the case's start date.
  const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
  const DAYS = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
  const pad2 = (n) => String(n).padStart(2, "0");
  const dmy = (d) => `${pad2(d.getDate())}-${MONTHS[d.getMonth()]}-${d.getFullYear()}`;
  const dm = (d) => `${pad2(d.getDate())}-${MONTHS[d.getMonth()]}`;
  const hm = (d) => `${pad2(d.getHours())}:${pad2(d.getMinutes())}`;
  function isoDate(d) { return `${d.getFullYear()}-${pad2(d.getMonth() + 1)}-${pad2(d.getDate())}`; }
  function parseDate(s) {
    const m = /^(\d{4})-(\d{2})-(\d{2})$/.exec(s || "");
    return m ? new Date(Number(m[1]), Number(m[2]) - 1, Number(m[3])) : null;
  }
  function tomorrow() { const d = new Date(); return new Date(d.getFullYear(), d.getMonth(), d.getDate() + 1); }
  // The period of a model, from its name (the MRPL cases: crude schedule by day, plan by week,
  // utilities by hour).
  function periodUnit(name) {
    const n = String(name || "").toUpperCase();
    return n.includes("CRUDE") ? "day" : n.includes("PLAN") ? "week" : n.includes("UTILITY") ? "hour" : null;
  }
  function calendar(c, n) {
    const unit = periodUnit(c.result?.summary?.name) || "day";
    const start = parseDate(c.options.startDate) || tomorrow();
    const at = (t) => unit === "hour" ? new Date(start.getFullYear(), start.getMonth(), start.getDate(), t)
      : new Date(start.getFullYear(), start.getMonth(), start.getDate() + (unit === "week" ? 7 * t : t));
    const cell = (t) => unit === "hour" ? { dateTime: at(t) } : at(t);  // for tables and Excel
    return {
      n, unit, start, at, cell,
      minWidth: unit === "day" ? 26 : unit === "week" ? 46 : 34,
      tick: (t) => { const d = at(t); return unit === "day" ? [String(d.getDate()), DAYS[d.getDay()].slice(0, 2)] : unit === "week" ? [dm(d), ""] : [hm(d), ""]; },
      weekend: (t) => unit === "day" && (at(t).getDay() === 0 || at(t).getDay() === 6),
      group: (t) => {
        const d = at(t), p = t > 0 ? at(t - 1) : null;
        if (unit === "hour") return !p || d.getDate() !== p.getDate() ? `${DAYS[d.getDay()]} ${dm(d)}` : "";
        return !p || d.getMonth() !== p.getMonth() ? `${MONTHS[d.getMonth()]} ${d.getFullYear()}` : "";
      },
      title: (t) => { const d = at(t); return unit === "hour" ? `${DAYS[d.getDay()]} ${dmy(d)} ${hm(d)}` : unit === "week" ? `week from ${DAYS[d.getDay()]} ${dmy(d)}` : `${DAYS[d.getDay()]} ${dmy(d)}`; },
      header: unit === "hour" ? "Date and hour" : unit === "week" ? "Week from" : "Date",
      // The horizon, short enough for a summary box: the dates, then the length and the year.
      span: () => {
        const a = at(0), b = unit === "week" ? new Date(at(n - 1).getFullYear(), at(n - 1).getMonth(), at(n - 1).getDate() + 6) : at(n - 1);
        const len = unit === "hour" ? `${n} hours` : unit === "week" ? `${n} weeks` : `${n} days`;
        const year = a.getFullYear() === b.getFullYear() ? String(a.getFullYear()) : `${a.getFullYear()}–${b.getFullYear()}`;
        if (unit === "hour" && isoDate(a) === isoDate(b)) return [`${dm(a)}, ${hm(a)}–${hm(b)}`, `${len} · ${year}`];
        return [`${dm(a)} – ${dm(b)}`, `${len} · ${year}`];
      },
    };
  }
  // A table or spreadsheet cell as text on the page.
  function cellText(v) {
    if (v instanceof Date) return `${DAYS[v.getDay()]} ${dmy(v)}`;
    if (v && v.dateTime instanceof Date) return `${DAYS[v.dateTime.getDay()]} ${dmy(v.dateTime)} ${hm(v.dateTime)}`;
    return typeof v === "number" ? num(v) : String(v ?? "");
  }

  // ---------------------------------------------------------------------------------------------
  // Host bridge.
  function send(msg) { if (host) host.postMessage(msg); else mockHost(msg); }
  if (host) host.addEventListener("message", (e) => onHost(e.data));
  function onHost(m) {
    switch (m.type) {
      case "init":
        state.samples = m.samples || []; state.gpu = !!m.gpu;
        $("#version").textContent = `samaya ${m.version || ""} · ${m.cores || "?"} cores${m.gpu ? " · NVIDIA GPU" : ""}`;
        $("#gpu-note").textContent = m.gpu ? "An NVIDIA GPU was found: PDLP can run its iterations on it." : "";
        renderSamples();
        break;
      case "picked": addCases(m.files || []); break;
      case "log": { const c = byId(m.id); if (c) { c.log.push(m.line); if (isShown(c)) appendLog(m.line); } break; }
      case "result":
        // A file samaya could not read has a summary but no model: show its message as an error.
        if (m.summary && m.summary.status === "read_error") finish(m.id, (c) => { c.status = "error"; c.error = m.summary.message || "The file could not be read."; });
        else finish(m.id, (c) => { c.status = "done"; c.result = m; });
        break;
      case "error": finish(m.id, (c) => { c.status = "error"; c.error = m.message; }); break;
      case "cancelled": finish(m.id, (c) => { c.status = "cancelled"; }); break;
      case "saved": toast(`Saved <b>${esc(fileName(m.path))}</b>`, m.path); break;
      case "notice": toast(esc(m.message)); break;
    }
  }
  function finish(id, fn) {
    const c = byId(id);
    if (!c) return;
    fn(c);
    c.elapsed = (performance.now() - c.started) / 1000;
    render();
    pump();
  }

  // ---------------------------------------------------------------------------------------------
  // Cases and the queue (one solve at a time; samaya itself uses all cores).
  const byId = (id) => state.cases.find((c) => c.id === id);
  function addCases(files) {
    const added = [];
    for (const f of files) {
      if (!/\.(mps|qps)$/i.test(f.path || f.name)) { toast(`Skipped <b>${esc(f.name)}</b>: samaya reads .mps and .qps files.`); continue; }
      const base = stem(f.name || f.path), same = state.cases.filter((x) => x.path === f.path).length;
      const c = { id: state.nextId++, name: same ? `${base} (${same + 1})` : base, file: f.name || f.path, path: f.path,
                  status: "queued", log: [], options: options() };
      state.cases.unshift(c);
      added.push(c);
    }
    if (added.length) { state.selected = added[added.length - 1].id; state.compare = null; state.rerun = null; state.tab = "overview"; }
    render();
    pump();
  }
  function pump() {
    if (state.cases.some((c) => c.status === "solving")) return;
    const next = [...state.cases].reverse().find((c) => c.status === "queued");
    if (!next) return;
    next.status = "solving";
    next.started = performance.now();
    send({ type: "solve", id: next.id, path: next.path, options: next.options });
    render();
  }
  function options() {
    return {
      startDate: $("#opt-start").value || isoDate(tomorrow()),
      timeLimit: Math.max(1, Number($("#opt-time").value) || 60),
      threads: Number($("#opt-threads").value) || 0,
      mipGap: Math.max(0, Number($("#opt-gap").value) || 0) / 100,
      lpMethod: $("#opt-method").value,
      presolve: $("#opt-presolve").checked,
    };
  }
  function setOptions(o) {
    $("#opt-start").value = o.startDate || isoDate(tomorrow());
    $("#opt-time").value = o.timeLimit;
    $("#opt-threads").value = String(o.threads);
    $("#opt-gap").value = +(o.mipGap * 100).toPrecision(6);
    $("#opt-method").value = o.lpMethod;
    $("#opt-presolve").checked = o.presolve;
  }
  const METHOD = { auto: "automatic", dual: "dual simplex", primal: "primal simplex", barrier: "interior point", pdlp: "PDLP" };
  function settingsRows(o) {
    return [["Start date", o.startDate ? cellText(parseDate(o.startDate)) : "–"], ["Time limit", `${o.timeLimit} s`],
            ["Threads", o.threads ? String(o.threads) : "automatic"], ["MIP gap", `${+(o.mipGap * 100).toPrecision(6)} %`],
            ["LP method", METHOD[o.lpMethod] || o.lpMethod], ["Presolve", o.presolve ? "on" : "off"]];
  }
  const settingsText = (o) => settingsRows(o).map(([k, v]) => `${k.toLowerCase()} ${v}`).join(", ");

  // Back from a case: the settings page, filled with that case's settings, to solve it again.
  function back(c) {
    state.selected = null;
    state.compare = null;
    state.rerun = c ? c.id : null;
    if (c) setOptions(c.options);
    render();
    $("#main").scrollTop = 0;
  }
  function rerun() {
    const c = byId(state.rerun);
    state.rerun = null;
    if (c) addCases([{ name: c.file || c.name, path: c.path }]);
    else render();
  }
  // Finished cases of the same model file as c, oldest first.
  const siblings = (c) => state.cases.filter((x) => x.path === c.path && x.status === "done").reverse();
  function compare(c) {
    if (!c || siblings(c).length < 2) return;
    state.compare = c.path;
    state.selected = null;
    render();
    $("#main").scrollTop = 0;
  }

  // ---------------------------------------------------------------------------------------------
  // Toolbar, sidebar, status bar.
  function renderToolbar() {
    const c = byId(state.selected);
    const done = !!c && c.status === "done";
    const cmp = state.compare ? state.cases.find((x) => x.path === state.compare && x.status === "done") : null;
    $("#tb-rerun").disabled = !c || c.status === "queued" || c.status === "solving";
    $("#tb-stop").disabled = !(c && c.status === "solving");
    $("#tb-compare").disabled = !(c && siblings(c).length >= 2);
    $("#tb-excel").disabled = !done && !cmp;
    $("#tb-csv").disabled = !done && !cmp;
    $("#tb-sol").disabled = !done;
    const solving = state.cases.find((x) => x.status === "solving");
    const queued = state.cases.filter((x) => x.status === "queued").length;
    $("#sb-state").textContent = solving ? `Solving ${solving.name}…${queued ? ` (${queued} queued)` : ""}` : "Ready";
  }
  function caseMeta(c) {
    if (c.status === "solving") return ["run", `Solving… ${secs((performance.now() - c.started) / 1000)}`];
    if (c.status === "queued") return ["", "Queued"];
    if (c.status === "error") return ["bad", "Error"];
    if (c.status === "cancelled") return ["warn", "Stopped"];
    const s = c.result.summary || {}, [label, cls] = statusOf(s.status);
    return [s.verified ? (cls === "ok" ? "ok" : "warn") : "bad", `${label}${s.verified ? ", verified" : ", not verified"} · ${secs(s.solve_seconds ?? c.elapsed)}`];
  }
  function renderCases() {
    const box = $("#cases");
    $("#case-count").textContent = state.cases.length ? String(state.cases.length) : "";
    if (!state.cases.length) { box.innerHTML = `<div class="empty-note">No cases yet. Open a model or pick a sample.</div>`; return; }
    box.innerHTML = state.cases.map((c) => {
      const [dot, meta] = caseMeta(c);
      const active = c.id === state.selected || (state.compare && c.path === state.compare && c.status === "done");
      const when = c.options.startDate ? ` · from ${dm(parseDate(c.options.startDate))}` : "";
      return `<div class="case-item ${active ? "active" : ""}" data-id="${c.id}" title="${esc(c.path)}\n${esc(settingsText(c.options))}">
        <span class="dot ${dot}"></span><div class="nm">${esc(c.name)}</div><div class="meta">${esc(meta)}</div>
        <div class="meta">${esc((METHOD[c.options.lpMethod] || c.options.lpMethod) + when)}</div></div>`;
    }).join("");
    box.querySelectorAll(".case-item").forEach((el) => el.addEventListener("click", () => {
      state.selected = Number(el.dataset.id); state.compare = null; state.rerun = null; state.tab = "overview"; render();
    }));
  }
  const SAMPLE_INFO = {
    "mrpl_crude_small.mps": ["MILP", "14 days"], "mrpl_crude_medium.mps": ["MILP", "21 days"],
    "mrpl_plan_small.mps": ["LP", "4 weeks"], "mrpl_plan_medium.mps": ["LP", "13 weeks"],
    "mrpl_utility_small.mps": ["MILP", "24 hours"], "mrpl_utility_medium.mps": ["MILP", "48 hours"],
  };
  function renderSamples() {
    const box = $("#samples");
    box.innerHTML = state.samples.length ? state.samples.map((s, i) => {
      const [type, horizon] = SAMPLE_INFO[s.name] || ["", ""];
      return `<tr><td><b>${esc(s.label || s.name)}</b><div class="muted">${esc(s.name)}</div></td><td>${esc(type)}</td><td>${esc(horizon)}</td>
        <td><button class="btn small" data-i="${i}">Solve</button></td></tr>`;
    }).join("") : `<tr><td colspan="4" class="muted">No samples were installed.</td></tr>`;
    box.querySelectorAll("button[data-i]").forEach((el) => el.addEventListener("click", () => addCases([state.samples[Number(el.dataset.i)]])));
  }

  // ---------------------------------------------------------------------------------------------
  // Main area.
  const isShown = (c) => state.selected === c.id;
  function render() {
    renderCases();
    renderToolbar();
    const c = byId(state.selected);
    const page = c || state.compare;
    $("#home").classList.toggle("hidden", !!page);
    $("#case").classList.toggle("hidden", !page);
    const again = !page && byId(state.rerun);
    $("#rerun").classList.toggle("hidden", !again);
    if (again) $("#rerun-name").textContent = again.name;
    if (c) { $("#case").innerHTML = casePage(c); wireCase(c); }
    else if (state.compare) { $("#case").innerHTML = comparePage(state.compare); wireCase(null); }
  }

  const BACK = `<button class="back" data-act="back" title="Back to the settings (Alt+Left)">
      <svg viewBox="0 0 24 24"><path d="M15 5l-7 7 7 7"/></svg>Back · change settings and run again</button>`;
  function header(c, tags) {
    return BACK + `<div class="case-head"><div class="title"><h1>${esc(c.result?.summary?.name || c.name)}</h1>
      <div class="path">${esc(c.path)}</div></div><div class="tags">${tags}</div></div>`;
  }
  const tag = (text, cls) => `<span class="tag ${cls}">${esc(text)}</span>`;
  const panel = (title, body, aside) => `<div class="panel"><div class="panel-head">${esc(title)}${aside ? `<span class="aside">${aside}</span>` : ""}</div>${body}</div>`;

  function casePage(c) {
    if (c.status === "queued") return header(c, tag("Queued", "info")) +
      panel("Waiting", `<div class="panel-body">This case starts when the current solve finishes. <button class="btn small" data-act="remove">Remove</button></div>`);
    if (c.status === "solving") return header(c, tag("Solving", "info")) +
      panel("Solving", `<div class="panel-body"><div class="muted">Elapsed</div><div class="elapsed" id="elapsed">0.0 s</div>
        <div class="progress"><div></div></div><div class="hint">Time limit ${c.options.timeLimit} s. The answer is verified before it is shown. Stop (Esc) ends the solve.</div></div>`) +
      panel("Solver log", `<div class="log" id="live-log">${esc(c.log.join("\n"))}</div>`);
    if (c.status === "error") return header(c, tag("Error", "bad")) +
      panel("samaya could not solve this file", `<div class="panel-body">${esc(c.error)}</div>`) + panel("Solver log", `<div class="log">${colorLog(c.log.join("\n"))}</div>`);
    if (c.status === "cancelled") return header(c, tag("Stopped", "warn")) +
      panel("Stopped", `<div class="panel-body">The solve was stopped before it finished.</div>`) + panel("Solver log", `<div class="log">${colorLog(c.log.join("\n"))}</div>`);
    return donePage(c);
  }
  function colorLog(text) { return esc(text).replace(/^(.*\(verified\).*)$/m, '<span class="ok">$1</span>'); }

  function kpiCells(c) {
    const R = c.result, s = R.summary || {};
    const gap = isFinite(s.objective) && isFinite(s.dual_bound) ? Math.abs(s.objective - s.dual_bound) / Math.max(1e-9, Math.abs(s.objective)) : NaN;
    const plan = planContent(c), span = plan.time ? plan.time.span() : null;
    return [
      ["Objective" + (R.sense ? ` (${R.sense})` : ""), num(s.objective), ""],
      ["Best bound", num(s.dual_bound), isFinite(gap) ? `gap ${(gap * 100).toFixed(3)} %` : ""],
      ["Horizon", span ? span[0] : "–", span ? span[1] : "no calendar for this model"],
      ["Solve time", secs(s.solve_seconds), `read ${secs(s.read_seconds)}`],
      ["Model size", `${num(s.rows, 12)} × ${num(s.cols, 12)}`, `${num(s.nnz, 12)} nonzeros`],
      ["Integer variables", num(s.integers ?? 0, 12), `${s.class || ""} · ${num(s.nodes ?? 0, 12)} nodes`],
    ];
  }
  function donePage(c) {
    const s = c.result.summary || {}, [label, cls] = statusOf(s.status);
    const tags = tag(label, cls) + (s.verified ? tag("Verified", "ok") : tag("Not verified", "bad"));
    const kpis = kpiCells(c).map(([l, v, sub]) => `<div class="kpi"><div class="l">${esc(l)}</div><div class="v" title="${esc(v)}">${esc(v)}</div><div class="s">${esc(sub)}</div></div>`).join("");
    const unit = periodUnit(s.name);
    const TABS = [["overview", "Overview"], ["plan", unit === "week" ? "Plan" : unit ? "Schedule" : "Results"], ["verification", "Verification"], ["solution", "Solution"], ["log", "Log"]];
    const tabs = TABS.map(([t, n]) => `<button class="tab ${state.tab === t ? "active" : ""}" data-tab="${t}">${n}</button>`).join("");
    return header(c, tags) + `<div class="kpis">${kpis}</div><div class="tabs">${tabs}</div><div id="tab-body">${tabBody(c)}</div>`;
  }

  // ---------------------------------------------------------------------------------------------
  // Tabs.
  function tabBody(c) {
    switch (state.tab) {
      case "plan": return planTab(c);
      case "verification": return verificationTab(c);
      case "solution": return solutionTab(c);
      case "log": return panel("Solver log", `<div class="log">${colorLog(c.result.log || c.log.join("\n"))}</div>`);
      default: return overviewTab(c);
    }
  }

  function checks(c) {
    const R = c.result, s = R.summary || {}, k = R.recheck || null, mip = (s.integers || 0) > 0;
    const own = [["Primal feasibility (rows and bounds)", s.max_primal_violation, TOL]];
    if (s.status === "optimal") own.push([mip ? "Dual feasibility (final LP)" : "Dual feasibility (optimality certificate)", s.max_dual_violation, TOL]);
    const ind = k && k.available !== false ? [
      ["Rows: worst relative violation", k.rows, TOL], ["Bounds: worst relative violation", k.bounds, TOL],
      ["Integrality: worst distance to an integer", k.integrality, TOL], ["Objective recomputed from the model (relative difference)", k.objectiveDiff, TOL],
    ] : [];
    return { own, ind };
  }
  const passes = (rows) => rows.every(([, v, t]) => typeof v === "number" && v <= t);
  function checkTable(rows) {
    return `<table><thead><tr><th>Check</th><th class="num">Largest error</th><th class="num">Tolerance</th><th class="num">Result</th></tr></thead><tbody>` +
      rows.map(([n, v, t]) => { const ok = typeof v === "number" && v <= t; return `<tr><td>${esc(n)}</td><td class="num">${num(v)}</td><td class="num">${t.toExponential(0)}</td><td class="num ${ok ? "pass" : "fail"}">${ok ? "Pass" : "Fail"}</td></tr>`; }).join("") +
      `</tbody></table>`;
  }
  function verdict(c) {
    const s = c.result.summary || {}, { own, ind } = checks(c);
    const ok = s.verified && passes(own) && passes(ind);
    const icon = ok ? `<svg viewBox="0 0 24 24"><path d="M20 6 9 17l-5-5"/></svg>` : `<svg viewBox="0 0 24 24"><path d="M12 8v5M12 17h.01"/><circle cx="12" cy="12" r="10"/></svg>`;
    const what = s.status === "optimal" ? "Optimal solution" : s.status === "infeasible" ? "The model is infeasible (certificate checked)"
      : s.status === "unbounded" ? "The model is unbounded (ray checked)" : `Status: ${statusOf(s.status)[0].toLowerCase()}`;
    return `<div class="verdict ${ok ? "ok" : "bad"}">${icon}<div><b>${esc(what)}.</b> ${ok ? "Verified by samaya and re-checked independently from the model file: every check passed." : "See the checks below."}</div></div>`;
  }
  function overviewTab(c) {
    const R = c.result, s = R.summary || {}, plan = planContent(c), first = plan.charts[0];
    const model = `<table class="kv"><tbody>
        <tr><td>File</td><td>${esc(fileName(c.path))}</td></tr>
        <tr><td>Type</td><td>${esc(s.class || "")}, ${esc(R.sense || "")}</td></tr>
        <tr><td>Rows · columns · nonzeros</td><td>${num(s.rows, 12)} · ${num(s.cols, 12)} · ${num(s.nnz, 12)}</td></tr>
        <tr><td>Simplex iterations · nodes</td><td>${num(s.simplex_iterations ?? 0, 12)} · ${num(s.nodes ?? 0, 12)}</td></tr>
        ${settingsRows(c.options).map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(v)}</td></tr>`).join("")}</tbody></table>`;
    return verdict(c) + `<div class="grid2">${panel("Verification", checkTable(checks(c).own.concat(checks(c).ind)) +
        `<div class="note">samaya’s verifier, then an independent re-check from the model file. Details in the Verification tab.</div>`)}
      ${panel("Case", model)}</div>` +
      (first ? `<div style="margin-top:12px">${panel(first[0], `<div class="chart-wrap">${first[1]}</div>`, "more in the " + (periodUnit(s.name) === "week" ? "Plan" : "Schedule") + " tab")}</div>` : "");
  }
  function verificationTab(c) {
    const k = c.result.recheck, { own, ind } = checks(c);
    return verdict(c) + `<div class="grid2">${panel("samaya’s built-in verifier", checkTable(own) +
        `<div class="note">Recorded by samaya before it reported the answer (its --json run summary).</div>`)}
      ${panel("Independent re-check by samaya Studio", (ind.length ? checkTable(ind) : `<div class="panel-body muted">No solution to re-check for this status.</div>`) +
        `<div class="note">Recomputed here from the model file and the solution alone, with separate code.${k && k.rowName ? ` Worst row: ${esc(k.rowName)}.` : ""}</div>`)}</div>`;
  }
  function solutionTab(c) {
    const R = c.result, q = state.filter.toLowerCase(), LIMIT = 2000;
    const cols = (R.cols || []).filter(([n, v]) => Math.abs(v) > 1e-9 && (!q || n.toLowerCase().includes(q)));
    const rows = (R.rowsData || []).filter(([n]) => !q || n.toLowerCase().includes(q));
    const t1 = cols.slice(0, LIMIT).map(([n, v]) => `<tr><td>${esc(n)}</td><td class="num">${num(v)}</td></tr>`).join("");
    const t2 = rows.slice(0, LIMIT).map(([n, a, d]) => `<tr><td>${esc(n)}</td><td class="num">${num(a)}</td><td class="num">${num(d)}</td></tr>`).join("");
    const more = (n) => n > LIMIT ? `<div class="note">Showing ${LIMIT} of ${num(n, 12)}; the Excel file has all of them.</div>` : "";
    return `<input class="search" id="filter" placeholder="Filter variables and constraints by name" value="${esc(state.filter)}">
      <div class="grid2" style="margin-top:0">${panel(`Variables (${num(cols.length, 12)} nonzero${q ? ", matching" : ""})`,
        `<div class="scroll"><table><thead><tr><th>Column</th><th class="num">Value</th></tr></thead><tbody>${t1}</tbody></table></div>${more(cols.length)}`)}
      ${panel(`Constraints (${num(rows.length, 12)}${q ? ", matching" : ""})`,
        `<div class="scroll"><table><thead><tr><th>Row</th><th class="num">Activity</th><th class="num">Dual</th></tr></thead><tbody>${t2}</tbody></table></div>${more(rows.length)}`)}</div>
      ${R.truncated ? `<div class="note">This model is large: Studio shows part of the solution; the solution file has every value.</div>` : ""}`;
  }

  // ---------------------------------------------------------------------------------------------
  // The plan and schedule, in the MRPL cases' terms (a port of cases/report.py), else generic views.
  // Series keyed by name from columns prefix<name>_<period>; all-zero series are left out unless
  // keepZero (a unit that stays off is still a row of the commitment chart).
  function series(cols, prefix, keepZero) {
    const out = {};
    for (const [k, v] of cols) {
      if (!k.startsWith(prefix)) continue;
      const i = k.lastIndexOf("_"); const key = k.slice(prefix.length, i), t = k.slice(i + 1);
      if (!/^\d+$/.test(t)) continue;
      (out[key] = out[key] || [])[Number(t)] = v;
    }
    for (const k of Object.keys(out)) { const a = out[k]; for (let i = 0; i < a.length; i++) a[i] = a[i] || 0; if (!keepZero && !a.some((v) => Math.abs(v) > 1e-6)) delete out[k]; }
    return out;
  }
  const periods = (d) => Math.max(0, ...Object.values(d).map((a) => a.length));
  const round3 = (v) => Math.round(v * 1000) / 1000 + 0;
  // A table with a row per period (dates down, as schedulers read them) and a column per series.
  function periodTable(title, data, time, unit) {
    const names = Object.keys(data).sort();
    return { title, rows: [[time.header, ...names.map((k) => `${pretty(k)}${unit ? ` (${unit})` : ""}`)],
      ...Array.from({ length: time.n }, (_, t) => [time.cell(t), ...names.map((k) => round3(data[k][t] || 0))])] };
  }
  function planContent(c) {
    const R = c.result;
    if (c._plan && c._plan.key === c.options.startDate) return c._plan;
    Charts.resetColors();
    const name = String(R.summary?.name || "").toUpperCase(), cols = R.cols || [];
    const charts = [], tables = [];
    let time = null;
    if (name.includes("CRUDE")) {
      const tank = series(cols, "tank_"), charge = series(cols, "charge_");
      time = calendar(c, periods(tank) || periods(charge) || 14);
      const arrivals = cols.filter(([k, v]) => (k.startsWith("arrive_") || k.startsWith("spot_")) && v > 0.5).map(([k]) => {
        const p = k.split("_"), spot = p[0] === "spot";
        return { t: Number(p[p.length - 1]), lane: spot ? "Spot cargoes" : "Term cargoes", cargo: `${spot ? "spot" : "term"} cargo ${p[1]}`,
                 short: `${spot ? "S" : "C"}${p[1]}`, crude: p.slice(2, -1).join("_") };
      }).sort((a, b) => a.t - b.t);
      const skipped = cols.filter(([k, v]) => k.startsWith("skip_") && v > 0.5).map(([k]) => { const p = k.slice(5).split("_"); return `cargo ${p[0]} (${pretty(p.slice(1).join("_"))})`; });
      charts.push(["Berth schedule (single-point mooring, at most one cargo a day)", Charts.gantt(arrivals, time, pretty)]);
      charts.push(["Crude tank stock, end of day", Charts.lines(tank, time, "kt", pretty)]);
      charts.push(["CDU charge by crude", Charts.stacked(charge, time, "kt/day", pretty)]);
      tables.push({ title: "Berth schedule", rows: [["Date", "Cargo", "Crude"], ...arrivals.map((a) => [time.cell(a.t), a.cargo, pretty(a.crude)])],
                    note: skipped.length ? "Not lifted: " + skipped.join(", ") : "" });
      tables.push(periodTable("Tank stock (kt, end of day)", tank, time));
      tables.push(periodTable("CDU charge (kt per day)", charge, time));
    } else if (name.includes("PLAN")) {
      const run = series(cols, "run_");
      time = calendar(c, periods(run) || 4);
      const feed = series(cols, "feed_"), sell = series(cols, "sell_"), exp = series(cols, "export_");
      charts.push(["Crude processed by week", Charts.stacked(run, time, "kt/week", pretty)]);
      charts.push(["Unit feed by week", Charts.lines(feed, time, "kt/week", pretty)]);
      const caps = (R.rowsData || []).filter(([k, , d]) => k.includes("capacity") && Math.abs(d) > 1e-6).map(([k, , d]) => [k, Math.abs(d)]).sort((a, b) => b[1] - a[1]);
      const capLabel = (k) => { const m = /^(.*)_capacity_(\d+)$/.exec(k); return m ? `${pretty(m[1])} capacity, ${time.title(Number(m[2]))}` : pretty(k); };
      if (caps.length) charts.push(["Marginal value of capacity (binding limits only)", Charts.hbars(caps.slice(0, 12), capLabel, "lakh Rs per extra kt")]);
      tables.push(periodTable("Crude processed (kt per week)", run, time));
      tables.push(periodTable("Unit feed (kt per week)", feed, time));
      tables.push(periodTable("Domestic sales (kt per week)", sell, time));
      tables.push(periodTable("Exports (kt per week)", exp, time));
      if (caps.length) tables.push({ title: "Marginal value of capacity", rows: [["Constraint", "Week from", "Lakh Rs per extra kt"],
        ...caps.map(([k, v]) => { const m = /^(.*)_capacity_(\d+)$/.exec(k); return [m ? pretty(m[1]) + " capacity" : k, m ? time.cell(Number(m[2])) : "", Math.round(v * 100) / 100]; })] });
    } else if (name.includes("UTILITY")) {
      const all = series(cols, "", true), on = {}, flows = {};
      for (const [k, v] of Object.entries(all)) { if (k.endsWith("_on")) on[k.slice(0, -3)] = v; else if ((k.endsWith("_steam") || k.endsWith("_mw")) && !k.startsWith("stg") && v.some((x) => Math.abs(x) > 1e-6)) flows[k] = v; }
      time = calendar(c, periods(on) || 24);
      charts.push(["Unit commitment (on or off, by hour)", Charts.onoff(on, time, pretty)]);
      charts.push(["Steam (t/h) and power (MW) by hour", Charts.lines(flows, time, "", pretty)]);
      tables.push(periodTable("Unit commitment (1 = on)", on, time));
      tables.push(periodTable("Steam (t/h) and power (MW)", flows, time));
    } else {
      const top = cols.filter(([, v]) => Math.abs(v) > 1e-9).sort((a, b) => Math.abs(b[1]) - Math.abs(a[1])).slice(0, 15);
      if (top.length) charts.push(["Largest variable values", Charts.hbars(top, String, "")]);
      const duals = (R.rowsData || []).filter(([, , d]) => Math.abs(d) > 1e-9).sort((a, b) => Math.abs(b[2]) - Math.abs(a[2])).slice(0, 12).map(([k, , d]) => [k, d]);
      if (duals.length) charts.push(["Largest duals (shadow prices)", Charts.hbars(duals, String, "")]);
    }
    c._plan = { key: c.options.startDate, charts, tables, time };
    return c._plan;
  }
  function planTab(c) {
    const { charts, tables } = planContent(c);
    if (!charts.length && !tables.length) return panel("Results", `<div class="panel-body muted">No plan view for this model; see the Solution tab.</div>`);
    return `<div class="stack">` + charts.map(([t, s]) => panel(t, `<div class="chart-wrap">${s}</div>`)).join("") +
      tables.map((t) => panel(t.title, `<div class="scroll">${simpleTable(t.rows)}</div>${t.note ? `<div class="note">${esc(t.note)}</div>` : ""}`, "also in the Excel export")).join("") + `</div>`;
  }
  function simpleTable(rows) {
    const [h, ...b] = rows;
    const numeric = h.map((_, i) => b.length > 0 && b.every((r) => typeof r[i] === "number"));
    return `<table><thead><tr>${h.map((x, i) => `<th class="${numeric[i] ? "num" : ""}">${esc(x)}</th>`).join("")}</tr></thead><tbody>` +
      b.map((r) => `<tr>${r.map((x, i) => `<td class="${numeric[i] ? "num" : ""}">${esc(cellText(x))}</td>`).join("")}</tr>`).join("") + `</tbody></table>`;
  }

  // ---------------------------------------------------------------------------------------------
  // Comparing the cases of one model (as planners compare a base case with alternatives).
  function compareData(path) {
    const list = state.cases.filter((x) => x.path === path && x.status === "done").reverse();
    const base = list[0], bs = base.result.summary || {};
    const rows = [["", ...list.map((c) => c.name)]];
    const line = (label, f) => rows.push([label, ...list.map(f)]);
    line("Status", (c) => statusOf(c.result.summary.status)[0]);
    line("Verified", (c) => c.result.summary.verified ? "yes" : "no");
    line("Objective", (c) => c.result.summary.objective ?? "");
    line(`Change from ${base.name}`, (c) => c === base ? "base case" : isFinite(c.result.summary.objective) && isFinite(bs.objective) ? c.result.summary.objective - bs.objective : "");
    line("Best bound", (c) => c.result.summary.dual_bound ?? "");
    line("Solve time (s)", (c) => Math.round((c.result.summary.solve_seconds ?? 0) * 1000) / 1000);
    line("Simplex iterations", (c) => c.result.summary.simplex_iterations ?? 0);
    line("Branch-and-bound nodes", (c) => c.result.summary.nodes ?? 0);
    for (const [k] of settingsRows(base.options)) line(k, (c) => settingsRows(c.options).find((r) => r[0] === k)[1]);
    // The variables that differ most between the cases.
    const vals = list.map((c) => new Map(c.result.cols || []));
    const keys = new Set(); vals.forEach((m) => m.forEach((_, k) => keys.add(k)));
    const diffs = [];
    for (const k of keys) {
      const v = vals.map((m) => m.get(k) || 0), spread = Math.max(...v) - Math.min(...v);
      if (spread > 1e-6) diffs.push([k, spread, v]);
    }
    diffs.sort((a, b) => b[1] - a[1]);
    const drows = [["Variable", ...list.map((c) => c.name)], ...diffs.slice(0, 25).map(([k, , v]) => [k, ...v.map(round3)])];
    return { list, base, rows, drows, ndiff: diffs.length };
  }
  function comparePage(path) {
    const { list, base, rows, drows, ndiff } = compareData(path);
    // The best objective, marked only when it is better than another case's (by more than 1e-9 relative).
    const sense = base.result.sense === "maximize" ? 1 : -1, objs = list.map((c) => c.result.summary.objective);
    const top = Math.max(...objs.filter(isFinite).map((v) => sense * v)), low = Math.min(...objs.filter(isFinite).map((v) => sense * v));
    const best = top - low > 1e-9 * Math.max(1, Math.abs(top)) ? list.find((c) => sense * c.result.summary.objective === top) : null;
    const head = `<tr>${rows[0].map((h, i) => `<th class="${i ? "num" : ""}">${esc(h)}${best && list[i - 1] === best ? " (best)" : ""}</th>`).join("")}</tr>`;
    // Rows whose values differ between the cases are in bold.
    const differs = (r) => r.slice(1).some((x) => cellText(x) !== cellText(r[1]));
    const body = rows.slice(1).map((r) => `<tr class="${!r[0].startsWith("Change") && !r[0].startsWith("Solve") && differs(r) ? "changed" : ""}">${r.map((x, i) => i === 0 ? `<td>${esc(x)}</td>` :
      `<td class="num ${list[i - 1] === best && r[0] === "Objective" ? "best" : ""} ${r[0].startsWith("Change") && typeof x === "number" ? (x > 1e-9 ? "diff-pos" : x < -1e-9 ? "diff-neg" : "") : ""}">${esc(typeof x === "number" && r[0].startsWith("Change") && x > 0 ? "+" + num(x) : cellText(x))}</td>`).join("")}</tr>`).join("");
    return BACK + `<div class="case-head"><div class="title"><h1>Compare cases: ${esc(base.result.summary.name || base.name)}</h1>
        <div class="path">${esc(path)}</div></div><div class="tags">${tag(`${list.length} cases`, "info")}</div></div>` +
      panel("Results and settings", `<table>${`<thead>${head}</thead>`}<tbody>${body}</tbody></table>`, `base case: ${esc(base.name)}`) +
      `<div style="margin-top:12px">${panel(ndiff ? `Largest differences in the solution (${num(Math.min(25, ndiff), 12)} of ${num(ndiff, 12)} variables that differ)` : "Differences in the solution",
        ndiff ? `<div class="scroll">${simpleTable(drows)}</div>` : `<div class="panel-body muted">The cases have the same solution.</div>`, "also in the Excel export")}</div>`;
  }

  // ---------------------------------------------------------------------------------------------
  // Excel and CSV.
  function sheets(c) {
    const R = c.result, s = R.summary || {}, { own, ind } = checks(c), plan = planContent(c);
    const span = plan.time ? plan.time.span() : null;
    const out = [["Summary", [["Item", "Value"], ["Model", s.name || c.name], ["Case", c.name], ["File", c.path], ["Sense", R.sense || ""],
      ["Status", s.status], ["Verified by samaya", s.verified ? "yes" : "no"], ["Objective", s.objective ?? ""], ["Best bound", s.dual_bound ?? ""],
      ["Horizon", span ? `${span[0]} (${span[1]})` : ""],
      ["Solve time (s)", s.solve_seconds ?? ""], ["Rows", s.rows], ["Columns", s.cols], ["Nonzeros", s.nnz], ["Integer columns", s.integers ?? 0],
      ["Simplex iterations", s.simplex_iterations ?? 0], ["Branch-and-bound nodes", s.nodes ?? 0],
      ...settingsRows(c.options), ["Exported", { dateTime: new Date() }]]]];
    out.push(["Verification", [["Check", "Checked by", "Largest error", "Tolerance", "Passed"],
      ...own.map(([n, v, t]) => [n, "samaya verifier", v, t, v <= t ? "yes" : "no"]),
      ...ind.map(([n, v, t]) => [n, "samaya Studio re-check (from the model file)", v, t, v <= t ? "yes" : "no"])]]);
    for (const t of plan.tables) out.push([t.title, t.note ? t.rows.concat([[t.note]]) : t.rows]);
    out.push(["Variables", [["Column", "Value"], ...(R.cols || [])]]);
    out.push(["Constraints", [["Row", "Activity", "Dual"], ...(R.rowsData || [])]]);
    out.push(["Solver log", [["samaya output"], ...String(R.log || "").split("\n").map((l) => [l])]]);
    return out;
  }
  function compareSheets(path) {
    const { rows, drows } = compareData(path);
    return [["Compare cases", rows], ["Largest differences", drows]];
  }
  function exportTarget() {
    const c = byId(state.selected);
    if (c && c.status === "done") return { name: c.name, sheets: () => sheets(c) };
    if (state.compare) { const b = state.cases.find((x) => x.path === state.compare); return { name: `${stem(b.path)}_compare`, sheets: () => compareSheets(state.compare) }; }
    return null;
  }
  function exportExcel() {
    const t = exportTarget();
    if (t) save(`${t.name}_results.xlsx`, Xlsx.workbook(t.sheets()), "excel");
  }
  function exportCsv() {
    const t = exportTarget();
    if (!t) return;
    const enc = new TextEncoder();
    const files = t.sheets().map(([n, rows]) => ({ name: `${n.replace(/[^\w.-]+/g, "_")}.csv`, data: enc.encode(Xlsx.csv(rows)) }));
    save(`${t.name}_csv.zip`, Xlsx.zip(files), "zip");
  }
  function save(name, bytes, kind) {
    if (host) { send({ type: "save", name, kind, base64: Xlsx.base64(bytes) }); return; }
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([bytes])); a.download = name; a.click();
    toast(`Downloaded <b>${esc(name)}</b>`);
  }

  // ---------------------------------------------------------------------------------------------
  // Wiring.
  function wireCase(c) {
    const page = $("#case");
    page.querySelectorAll("[data-tab]").forEach((b) => b.addEventListener("click", () => { state.tab = b.dataset.tab; render(); }));
    page.querySelectorAll("[data-act]").forEach((b) => b.addEventListener("click", () => {
      const a = b.dataset.act;
      if (a === "back") back(c || state.cases.find((x) => x.path === state.compare && x.status === "done"));
      if (a === "remove" && c) { state.cases = state.cases.filter((x) => x !== c); state.selected = null; render(); }
    }));
    const f = page.querySelector("#filter");
    if (f && c) { f.addEventListener("input", () => { state.filter = f.value; const pos = f.selectionStart; $("#tab-body").innerHTML = tabBody(c); wireCase(c); const g = $("#filter"); g.focus(); g.setSelectionRange(pos, pos); }); }
    const live = page.querySelector("#live-log"); if (live) live.scrollTop = live.scrollHeight;
  }
  function appendLog(line) {
    const live = $("#live-log");
    if (!live) return;
    live.textContent += (live.textContent ? "\n" : "") + line;
    live.scrollTop = live.scrollHeight;
  }
  setInterval(() => {
    const c = state.cases.find((x) => x.status === "solving");
    if (!c) return;
    const t = secs((performance.now() - c.started) / 1000);
    const e = $("#elapsed"); if (e && isShown(c)) e.textContent = t;
    const item = document.querySelector(`.case-item[data-id="${c.id}"] .meta`);
    if (item) item.textContent = `Solving… ${t}`;
  }, 100);

  function toast(html, path) {
    const t = document.createElement("div");
    t.className = "toast";
    t.innerHTML = `<div class="msg">${html}</div>` + (path ? `<a data-p="open">Open</a><a data-p="reveal">Show in folder</a>` : "");
    t.querySelectorAll("a").forEach((a) => a.addEventListener("click", () => send({ type: a.dataset.p, path })));
    $("#toasts").appendChild(t);
    setTimeout(() => t.remove(), 7000);
  }

  // Theme: light by default (office software); the choice is kept for the next start.
  function setTheme(dark) {
    document.documentElement.dataset.theme = dark ? "dark" : "light";
    $("#theme-btn").textContent = dark ? "Light theme" : "Dark theme";
    try { localStorage.setItem("studio.theme", dark ? "dark" : "light"); } catch (e) { /* storage may be off */ }
    if (host) send({ type: "theme", dark });
  }
  let stored = null;
  try { stored = localStorage.getItem("studio.theme"); } catch (e) { /* storage may be off */ }

  const pick = () => send({ type: "pick" });
  $("#tb-open").addEventListener("click", pick);
  $("#browse-link").addEventListener("click", (e) => { e.preventDefault(); pick(); });
  $("#tb-new").addEventListener("click", () => back(null));
  $("#home-btn").addEventListener("click", () => back(null));
  $("#tb-rerun").addEventListener("click", () => back(byId(state.selected)));
  $("#tb-stop").addEventListener("click", () => { const c = byId(state.selected); if (c) send({ type: "cancel", id: c.id }); });
  $("#tb-compare").addEventListener("click", () => compare(byId(state.selected)));
  $("#tb-excel").addEventListener("click", exportExcel);
  $("#tb-csv").addEventListener("click", exportCsv);
  $("#tb-sol").addEventListener("click", () => { const c = byId(state.selected); if (c) send({ type: "saveCopy", id: c.id, name: `${c.name}.sol` }); });
  $("#rerun-btn").addEventListener("click", rerun);
  $("#rerun-close").addEventListener("click", () => { state.rerun = null; render(); });
  $("#theme-btn").addEventListener("click", () => setTheme(document.documentElement.dataset.theme !== "dark"));
  document.addEventListener("keydown", (e) => {
    const c = byId(state.selected);
    if (e.ctrlKey && e.key.toLowerCase() === "o") { e.preventDefault(); pick(); }
    if (e.ctrlKey && e.key.toLowerCase() === "s") { e.preventDefault(); exportExcel(); }
    if (e.key === "Escape" && c && c.status === "solving") send({ type: "cancel", id: c.id });
    if (e.altKey && e.key === "ArrowLeft" && (c || state.compare)) { e.preventDefault(); back(c || state.cases.find((x) => x.path === state.compare && x.status === "done")); }
  });

  // Drag and drop anywhere in the window. The files go to the host as objects, which gives it
  // their real paths (WebView2 postMessageWithAdditionalObjects).
  let depth = 0;
  document.addEventListener("dragenter", (e) => { e.preventDefault(); depth++; document.body.classList.add("dragging"); $("#drop").classList.add("over"); });
  document.addEventListener("dragleave", () => { if (--depth <= 0) { depth = 0; document.body.classList.remove("dragging"); $("#drop").classList.remove("over"); } });
  document.addEventListener("dragover", (e) => e.preventDefault());
  document.addEventListener("drop", (e) => {
    e.preventDefault(); depth = 0; document.body.classList.remove("dragging"); $("#drop").classList.remove("over");
    const files = e.dataTransfer?.files;
    if (!files || !files.length) return;
    if (host && host.postMessageWithAdditionalObjects) host.postMessageWithAdditionalObjects({ type: "dropped" }, files);
    else addCases([...files].map((f) => ({ name: f.name, path: f.name })));
  });

  // ---------------------------------------------------------------------------------------------
  // Mock host (plain browser): replays a recorded solve of the crude scheduling case.
  function mockHost(msg) {
    const M = window.SAMAYA_MOCK;
    if (msg.type === "ready") setTimeout(() => onHost({ type: "init", version: "0.1.0", cores: 16, gpu: false,
      samples: [{ name: "mrpl_crude_small.mps", label: "Crude receipt schedule", path: "samples\\mrpl_crude_small.mps" },
                { name: "mrpl_plan_small.mps", label: "Refinery plan", path: "samples\\mrpl_plan_small.mps" },
                { name: "mrpl_utility_small.mps", label: "Steam and power", path: "samples\\mrpl_utility_small.mps" }] }), 0);
    if (msg.type === "pick") onHost({ type: "picked", files: [{ name: "mrpl_crude_small.mps", path: M.path }] });
    if (msg.type === "solve") {
      const lines = String(M.log).split("\n");
      if (demo.has("instant")) {  // screenshots: no timers
        lines.forEach((l) => onHost({ type: "log", id: msg.id, line: l }));
        if (!demo.has("solving")) onHost(Object.assign({}, M, { id: msg.id }));
        return;
      }
      lines.forEach((l, i) => setTimeout(() => onHost({ type: "log", id: msg.id, line: l }), 40 * i));
      setTimeout(() => onHost(Object.assign({}, M, { id: msg.id })), 40 * lines.length + 200);
    }
    if (msg.type === "cancel") onHost({ type: "cancelled", id: msg.id });
  }

  // Test hooks for the host's screenshot mode (--capture, --eval): open a tab or a page directly.
  window.studio = {
    select(tab) { state.tab = tab; render(); },
    home() { back(null); },
    back() { back(byId(state.selected)); },
    rerun,
    compare() { compare(byId(state.selected) || state.cases.find((x) => x.status === "done")); },
    add(files) { addCases(files); },
    options(o) { setOptions(Object.assign(options(), o)); },
    scroll(y) { $("#main").scrollTop = y; },
    theme(dark) { setTheme(!!dark); },
    // The workbook that Export to Excel would save, as base64 (for automated checks).
    workbookBase64() { const t = exportTarget(); return t ? Xlsx.base64(Xlsx.workbook(t.sheets())) : ""; },
    state,
  };

  // Mock-mode demo hook for screenshots: index.html#demo&instant&tab=plan
  const demo = new URLSearchParams(location.hash.slice(1));
  $("#opt-start").value = isoDate(tomorrow());
  setTheme(stored === "dark");
  render();
  send({ type: "ready" });
  if (!host && demo.has("demo")) {
    setTimeout(() => {
      if (demo.get("start")) $("#opt-start").value = demo.get("start");
      addCases([{ name: "mrpl_crude_small.mps", path: window.SAMAYA_MOCK.path }]);
      if (demo.get("tab")) window.studio.select(demo.get("tab"));
      if (demo.has("back")) window.studio.back();
    }, 0);
  }
})();
