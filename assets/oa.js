// The supplemental appendix's interactive figures: Observable Plot over the datasets of python/paper/siteData.py, in
// the palette of the paper's figures (_generated/style.js, written from python/paper/figures.py). Plot, d3 and
// Inputs are pinned here, so the site draws the same in a year as today.
import * as Plot from "https://cdn.jsdelivr.net/npm/@observablehq/plot@0.6.16/+esm";
import * as d3 from "https://cdn.jsdelivr.net/npm/d3@7.9.0/+esm";
import * as Inputs from "https://cdn.jsdelivr.net/npm/@observablehq/inputs@0.12.0/+esm";
import {STYLE} from "../_generated/style.js";

export {Plot, d3, Inputs, STYLE};
const INK = STYLE.ink, SERIES = STYLE.series;
const near = (a, b) => Math.abs(a - b) < 1e-9;

// ---------------------------------------------------------------------------------------------------------
// The selection in the address: #key=value&..., so a link opens the same view.

export function readState(key, fallback) {
  const v = new URLSearchParams(location.hash.slice(1)).get(key);
  return v === null ? fallback : v;
}

export function writeState(obj) {
  const p = new URLSearchParams(location.hash.slice(1));
  for (const [k, v] of Object.entries(obj)) p.set(k, v);
  history.replaceState(null, "", "#" + p.toString());
}

// A control whose value is kept in the address under `key`. options: [[label, value], ...]; the value given in
// the address wins over `value` when it is one of the options.
export function control(kind, options, {key, label, value} = {}) {
  const values = options.map(o => o[1]);
  const fromHash = readState(key, null);
  const start = values.find(v => String(v) === fromHash) ?? value ?? values[0];
  const make = kind === "radio" ? Inputs.radio : Inputs.select;
  const input = make(new Map(options), {label, value: start});
  input.classList.add("oa-control");
  input.addEventListener("input", () => writeState({[key]: input.value}));
  return input;
}

export function checkboxes(options, {key, label, value} = {}) {
  const values = options.map(o => o[1]);
  const fromHash = readState(key, null);
  const start = fromHash !== null ? fromHash.split(",").filter(v => values.includes(v)) : (value ?? values);
  const input = Inputs.checkbox(new Map(options), {label, value: start});
  input.classList.add("oa-control");
  input.addEventListener("input", () => writeState({[key]: input.value.join(",")}));
  return input;
}

// ---------------------------------------------------------------------------------------------------------
// Units and formats, as the paper prints them.

export const OUTCOMES = {
  tau:      {label: "Equilibrium tax rate", axis: "Percentage points", level: "Tax rate (%)",
             change: v => signed(v, 2) + " p.p.", fmt: v => d3.format(".2f")(v) + "%"},
  sr:       {label: "Savings over GDP", axis: "Percentage points of GDP", level: "Savings over GDP (%)",
             change: v => signed(v, 2) + " p.p.", fmt: v => d3.format(".2f")(v) + "%"},
  workweek: {label: "Average workweek", axis: "Hours per week", level: "Average workweek (hours)",
             change: v => signed(v, 2) + " hours", fmt: v => d3.format(".2f")(v) + " hours"},
  theta:    {label: "Pension design θ", axis: "θ in 2020", level: "Pension design θ",
             change: v => signed(v, 3), fmt: v => d3.format(".3f")(v)},
  nu:       {label: "Workers per retiree", axis: null, level: "Workers per retiree ν",
             change: v => signed(v, 3), fmt: v => d3.format(".3f")(v)},
  iota:     {label: "Informal savings ratio ι", axis: null, level: "Informal savings ratio ι",
             change: v => signed(100*v, 1) + " p.p.", fmt: v => d3.format(".3f")(v)}
};

export function signed(v, digits) {
  const r = +v.toFixed(digits);
  return r === 0 ? (0).toFixed(digits) : (r > 0 ? "+" : "−") + Math.abs(r).toFixed(digits);
}

// The light tone of a colour: 60% of it over white, as figures.py's TINT.
function tint(c) {
  const a = d3.rgb(c), w = 0.6;
  return d3.rgb(w*a.r + (1-w)*255, w*a.g + (1-w)*255, w*a.b + (1-w)*255).formatHex();
}

// Three values of a parameter in the paper's three colours, the middle one being the paper's.
export function seriesColours(values) {
  return new Map(values.map((v, i) => [v, STYLE.rhoColours[i % 3]]));
}

export function legend(items) {
  const div = document.createElement("div");
  div.className = "oa-legend";
  for (const {label, colour, shape = "square", dashed = false} of items) {
    const s = document.createElement("span");
    const sw = shape === "line"
      ? `<svg width="22" height="10"><line x1="1" y1="5" x2="21" y2="5" stroke="${colour}" stroke-width="2.2" ${dashed ? 'stroke-dasharray="4,3"' : ""}/></svg>`
      : shape === "open"
        ? `<svg width="12" height="12"><circle cx="6" cy="6" r="4.2" fill="white" stroke="${colour}" stroke-width="1.6"/></svg>`
        : shape === "dot"
          ? `<svg width="12" height="12"><circle cx="6" cy="6" r="4.6" fill="${colour}"/></svg>`
          : `<svg width="14" height="10"><rect width="14" height="10" fill="${colour}"/></svg>`;
    s.innerHTML = sw + `<span>${label}</span>`;
    div.append(s);
  }
  return div;
}

// A box drawn at the width of its container and redrawn when that changes. OJS's width is main's, which spans
// the page once the page holds a .column-page block, so it serves only as the first guess. draw(w) returns a node
// or an array of nodes.
function fitted(draw, width) {
  const box = document.createElement("div");
  let drawn = Math.floor(width);
  box.replaceChildren(...[draw(drawn)].flat());
  const ro = new ResizeObserver(() => {
    if (!box.isConnected) return ro.disconnect();
    const w = Math.floor(box.clientWidth);
    if (w > 0 && Math.abs(w - drawn) > 1) { drawn = w; box.replaceChildren(...[draw(w)].flat()); }
  });
  ro.observe(box);
  return box;
}

// ---------------------------------------------------------------------------------------------------------
// The counterfactual figure: rows of counterfactuals, two facets side by side on one axis, colour the series
// (rho or xi). mode "given": figure 6.1's bars of the change, the design row in two tones (theta = 0 full,
// theta = 1 light) and ageing as the acute bar with a tick for mild. mode "chosen": figure 7.1's dumbbells, the
// open marker pinned and the filled one chosen; the design as a level with the corners and each facet's
// observed design drawn.
//   rows:    [{label, scenario, light?, tick?}]   scenario is the dataset key; light/tick name the second key
//   facets:  [{label, filter}]                    filter(d) selects the facet's records
//   series:  {field, values, labels}              the colour dimension
//   observed: {facetLabel: theta}                 mode "chosen", outcome theta
export function counterfactuals(data, {rows, facets, series, outcome, mode = "given", observed = {}, width = 640,
                                       note = null}) {
  const o = OUTCOMES[outcome];
  const colours = seriesColours(series.values);
  const label = new Map(series.values.map((v, i) => [v, series.labels[i]]));
  const marks = [], bars = [], ticks = [], pins = [], unchanged = [];
  const level = mode === "chosen" && outcome === "theta";
  const value = d => level ? d.level : d.change;
  for (const f of facets) {
    const sub = data.filter(d => f.filter(d) && d.outcome === outcome);
    for (const r of rows) {
      let vals = [];
      for (const s of series.values) {
        const at = sub.filter(d => near(d[series.field], s));
        if (mode === "given") {
          const main = at.find(d => d.scenario === r.scenario);
          const second = r.light ?? r.tick;
          const other = second ? at.find(d => d.scenario === second) : null;
          if (main) { bars.push({facet: f.label, row: r.label, series: s, value: main.change, fill: colours.get(s),
                                 text: `${r.label}${r.light ? " (θ = 0)" : r.tick ? " (acute)" : ""}, ${label.get(s)}: ${o.change(main.change)}`});
                      vals.push(main.change); }
          if (other && r.light) { bars.push({facet: f.label, row: r.label, series: s, value: other.change,
                                             fill: tint(colours.get(s)),
                                             text: `${r.label} (θ = 1), ${label.get(s)}: ${o.change(other.change)}`});
                                  vals.push(other.change); }
          if (other && r.tick) { ticks.push({facet: f.label, row: r.label, series: s, value: other.change,
                                             stroke: colours.get(s),
                                             text: `${r.label} (mild), ${label.get(s)}: ${o.change(other.change)}`});
                                 vals.push(other.change); }
        } else {
          const pin = at.find(d => d.scenario === r.scenario && d.reading === "pinned");
          const cho = at.find(d => d.scenario === r.scenario && d.reading === "chosen");
          if (pin && cho) {
            pins.push({facet: f.label, row: r.label, series: s, pinned: value(pin), chosen: value(cho),
                       stroke: colours.get(s),
                       text: `${r.label}, ${label.get(s)}: ${level ? "θ " : ""}${(level ? o.fmt : o.change)(value(pin))} pinned, ${(level ? o.fmt : o.change)(value(cho))} chosen`});
            vals.push(value(pin) - (level ? (observed[f.label] ?? 0) : 0), value(cho) - (level ? (observed[f.label] ?? 0) : 0));
          }
        }
      }
      const zero = level ? 5e-4 : 0.05;
      if (vals.length && vals.every(v => Math.abs(v) < zero)) unchanged.push({facet: f.label, row: r.label});
      if (!vals.length) unchanged.push({facet: f.label, row: r.label, text: "not run"});
    }
  }
  bars.sort((a, b) => Math.abs(b.value) - Math.abs(a.value));     // the shorter bar drawn over the longer
  const xs = [0, ...bars.map(d => d.value), ...ticks.map(d => d.value),
              ...pins.flatMap(d => [d.pinned, d.chosen]), ...(level ? [0, 1] : [])];
  const nS = series.values.length;
  if (level) {
    marks.push(Plot.ruleX([0, 1], {stroke: INK.grid, strokeWidth: 1.5}));
    marks.push(Plot.ruleX(facets.filter(f => observed[f.label] !== undefined),
                          {x: f => observed[f.label], fx: f => f.label, stroke: INK.muted, strokeDasharray: "3,3"}));
  } else {
    marks.push(Plot.ruleX([0], {stroke: INK.primary, strokeWidth: 1}));
  }
  marks.push(Plot.barX(bars, {x1: 0, x2: "value", y: "series", fx: "facet", fy: "row", fill: "fill",
                              title: "text", tip: true}));
  marks.push(Plot.tickX(ticks, {x: "value", y: "series", fx: "facet", fy: "row", stroke: "stroke", strokeWidth: 4}));
  marks.push(Plot.tickX(ticks, {x: "value", y: "series", fx: "facet", fy: "row", stroke: "white", strokeWidth: 2,
                                title: "text", tip: true}));
  marks.push(Plot.link(pins, {x1: "pinned", x2: "chosen", y1: "series", y2: "series", fx: "facet", fy: "row",
                              stroke: "stroke", strokeWidth: 2}));
  marks.push(Plot.dot(pins, {x: "pinned", y: "series", fx: "facet", fy: "row", r: 4, fill: "white",
                             stroke: "stroke", strokeWidth: 1.6}));
  marks.push(Plot.dot(pins, {x: "chosen", y: "series", fx: "facet", fy: "row", r: 4.4, fill: "stroke",
                             title: "text", tip: true}));
  marks.push(Plot.text(unchanged, {fx: "facet", fy: "row", frameAnchor: "left", textAnchor: "start", dx: 6,
                                   text: d => d.text ?? "unchanged", fill: INK.muted, fontStyle: "italic", fontSize: 12}));
  // One chart per facet, stacked, on one x scale, so the facets compare by position: the grid runs through every
  // chart and the axis is drawn once, under the last.
  const charts = w => facets.flatMap((f, k) => {
    const last = k === facets.length - 1;
    const title = document.createElement("div");
    title.className = "oa-facet-title";
    title.textContent = f.label;
    return [title, Plot.plot({
      width: w, height: rows.length * (nS * 10 + 14) + (last ? 50 : 8), marginLeft: 150, marginRight: 14, marginTop: 4,
      marginBottom: last ? 46 : 4, style: {fontSize: "12px"},
      x: {domain: d3.extent(xs), nice: true, insetLeft: 4, insetRight: 4, axis: null},
      y: {domain: series.values, axis: null, padding: 0.12},
      fx: {domain: [f.label], axis: null},
      fy: {domain: rows.map(r => r.label), label: null, padding: 0.18},
      color: {type: "identity"},
      marks: [Plot.gridX(), ...(last ? [Plot.axisX({label: o.axis, labelAnchor: "center", labelOffset: 36})] : []),
              ...marks]
    })];
  });
  const wrap = document.createElement("div");
  wrap.className = "oa-chart";
  const items = series.values.map((v, i) => ({label: series.labels[i], colour: colours.get(v)}));
  if (mode === "chosen") items.push({label: "pinned", colour: INK.secondary, shape: "open"},
                                    {label: "chosen", colour: INK.secondary, shape: "dot"});
  wrap.append(legend(items), fitted(charts, width));
  if (note) { const p = document.createElement("div"); p.className = "oa-chartnote"; p.textContent = note; wrap.append(p); }
  return wrap;
}

// ---------------------------------------------------------------------------------------------------------
// Paths over calendar years: one panel per variable, model periods as markers joined by lines. lines:
// [{label, colour, dashed?, symbol?, records}] with records [{year, value}]. rule: a year to mark.

export function pathPanel(variable, lines, {width = 320, height = 250, rule = null, ruleLabel = null,
                                             observed = []} = {}) {
  const o = OUTCOMES[variable] ?? {fmt: d3.format(".3f"), level: variable};
  const pts = lines.flatMap(l => l.records.filter(r => Number.isFinite(r.value)).map(r =>
    ({...r, key: l.label, colour: l.colour, symbol: l.symbol ?? "circle", dashed: !!l.dashed,
      text: `${l.label}, ${r.year}: ${o.fmt(r.value)}`})));
  const ys = pts.map(d => d.value).concat(observed.map(d => d.value));
  const marks = [];
  if (rule !== null) {
    marks.push(Plot.ruleX([rule], {stroke: INK.grid, strokeWidth: 1.5}));
  }
  for (const l of lines) {
    const recs = l.records.filter(r => Number.isFinite(r.value));
    marks.push(Plot.line(recs, {x: "year", y: "value", stroke: l.colour, strokeWidth: 2,
                                strokeDasharray: l.dashed ? "5,4" : null}));
  }
  marks.push(Plot.dot(pts, {x: "year", y: "value", fill: d => d.dashed ? "white" : d.colour, stroke: "colour",
                            symbol: "symbol", r: 3.6, strokeWidth: 1.5, title: "text", tip: true}));
  if (observed.length) marks.push(Plot.dot(observed, {x: "year", y: "value", r: 6, stroke: d => d.colour,
                                                      fill: "none", strokeWidth: 1.4, title: d => d.text, tip: true}));
  return Plot.plot({
    width, height, marginLeft: 48, marginTop: 18, style: {fontSize: "12px"},
    x: {label: null, tickFormat: "d", ticks: [...new Set(pts.map(d => d.year))]},
    y: {label: null, grid: true, nice: true, domain: d3.extent(ys)},
    color: {type: "identity"}, symbol: {type: "identity"},
    marks
  });
}

// Two or more path panels in a row, each titled, with one legend above; one column in a box narrower than 560px.
export function pathGrid(panels, items, {width = 640, columns = 2, height = 250} = {}) {
  const wrap = document.createElement("div");
  wrap.className = "oa-chart";
  wrap.append(legend(items), fitted(width => {
    const n = width < 560 ? 1 : columns;
    const grid = document.createElement("div");
    grid.className = "oa-grid";
    grid.style.gridTemplateColumns = `repeat(${n}, minmax(0, 1fr))`;
    const w = Math.max(240, Math.floor((width - 24*(n - 1)) / n));
    for (const p of panels) {
      const box = document.createElement("div");
      const h = document.createElement("div");
      h.className = "oa-panel-title";
      h.textContent = p.title;
      box.append(h, pathPanel(p.variable, p.lines, {width: w, height, rule: p.rule, ruleLabel: p.ruleLabel,
                                                    observed: p.observed ?? []}));
      if (p.readout) { const r = document.createElement("div"); r.className = "oa-readout"; r.textContent = p.readout; box.append(r); }
      grid.append(box);
    }
    return grid;
  }, width));
  return wrap;
}

// ---------------------------------------------------------------------------------------------------------
// A2: the Argentine reform's paths at one rho, the paper's rho = 1 faint behind it, and under each panel the
// changes in 2010 and 2040 -- the two points of figure 5.2 at that rho.

const ARGPANELS = [
  {variable: "tau", title: "Tax rate τ (%)", change: v => signed(v, 2) + " p.p."},
  {variable: "sr", title: "Savings over GDP (%)", change: v => signed(v, 2) + " p.p."},
  {variable: "workweek", title: "Average workweek (hours)", change: v => signed(v, 2) + " hours"},
  {variable: "iota", title: "Informal savings ratio ι", change: v => signed(100*v, 1) + " p.p."}
];

export function reformPaths(data, rho, {width = 640} = {}) {
  const pre = INK.muted, post = SERIES[0];
  const panels = ARGPANELS.map(p => {
    const at = r => data.filter(d => d.variable === p.variable && near(d.rho, r)).sort((a, b) => a.year - b.year);
    const sel = at(rho), ref = near(rho, 1) ? [] : at(1);
    const lines = [{label: `Pre-reform, ρ = ${rho.toFixed(1)}`, colour: pre, records: sel.map(d => ({year: d.year, value: d.pre}))},
                   {label: `Reform, ρ = ${rho.toFixed(1)}`, colour: post, records: sel.map(d => ({year: d.year, value: d.reform}))}];
    if (ref.length) lines.unshift({label: "Pre-reform, ρ = 1.0", colour: d3.color(pre).copy({opacity: 0.35}).formatRgb(), dashed: true,
                                   records: ref.map(d => ({year: d.year, value: d.pre}))},
                                  {label: "Reform, ρ = 1.0", colour: d3.color(post).copy({opacity: 0.35}).formatRgb(), dashed: true,
                                   records: ref.map(d => ({year: d.year, value: d.reform}))});
    const g = y => sel.find(d => d.year === y);
    const readout = `Change: ${p.change(g(2010).reform - g(2010).pre)} in 2010, ${p.change(g(2040).reform - g(2040).pre)} in 2040`;
    return {title: p.title, variable: p.variable, lines, rule: 2010, ruleLabel: "reform", readout};
  });
  const items = [{label: "pre-reform", colour: pre, shape: "line"}, {label: "reform", colour: post, shape: "line"}];
  if (!near(rho, 1)) items.push({label: "ρ = 1, the paper's", colour: "#b8b7b2", shape: "line", dashed: true});
  return pathGrid(panels, items, {width, columns: 2, height: 230});
}

// ---------------------------------------------------------------------------------------------------------
// The rows and series the pages share: section 6's counterfactuals in figure 6.1's order, section 7's in figure
// 7.1's, the two parameter grids and section 7's two costs.

export const ROWS_GIVEN = [
  {label: "Pension design", scenario: "theta0", light: "theta1"},
  {label: "Ageing", scenario: "acute", tick: "mild"},
  {label: "French income distr.", scenario: "income"},
  {label: "French voting", scenario: "voting"},
  {label: "French leisure", scenario: "leisure"}
];
export const ROWS_CHOSEN = [
  {label: "Acute ageing", scenario: "acute"},
  {label: "French income distr.", scenario: "income"},
  {label: "French voting", scenario: "voting"},
  {label: "French leisure", scenario: "leisure"}
];
export const RHO = {field: "rho", values: [0.5, 1, 2], labels: ["ρ = 0.5", "ρ = 1, the paper's", "ρ = 2"]};
export const XI = {field: "xi", values: [0.2, 0.3, 0.4], labels: ["ξ = 0.2", "ξ = 0.3, the paper's", "ξ = 0.4"]};
export const GIVEN_OUTCOMES = [["Tax rate", "tau"], ["Savings over GDP", "sr"], ["Workweek", "workweek"]];
export const CHOSEN_OUTCOMES = [["Pension design θ", "theta"], ["Tax rate", "tau"], ["Savings over GDP", "sr"],
                                ["Workweek", "workweek"]];
// Section 7's two costs of redistribution, the `cost` column of the endogenous-design datasets: linear in the tax
// rate (sections 7.1-7.2, the default) and with the tax rate squared (section 7.3).
export const COSTFORMS = [["linear", "linear"], ["squared", "squared"]];
export const COSTNAMES = {linear: "The linear cost", squared: "The squared cost"};
export const COUNTRIES = [["U.S.", "US"], ["UK", "UK"], ["France", "FR"]];
const COUNTRYSTYLE = {US: {colour: SERIES[0], symbol: "circle"}, UK: {colour: SERIES[1], symbol: "square"},
                      FR: {colour: INK.secondary, symbol: "triangle"}};
export function countryStyle(c) { return COUNTRYSTYLE[c]; }
