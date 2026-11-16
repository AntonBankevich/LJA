#include "viewer_page.hpp"

namespace dot_viewer {
    const char *const VIEWER_PAGE = R"PAGE(<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>Dot viewer</title>
<style>
  * { box-sizing: border-box; }
  body { margin: 0; font: 14px system-ui, sans-serif; display: flex; flex-direction: column; height: 100vh; background: #fafafa; }
  #controls { padding: 8px 12px; border-bottom: 1px solid #ccc; background: #fff; display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }
  #controls textarea { width: 380px; height: 34px; font: 13px monospace; resize: vertical; }
  #controls input[type=number] { width: 70px; }
  #controls label { display: flex; align-items: center; gap: 4px; }
  .sep { width: 1px; align-self: stretch; background: #ddd; }
  #nav { padding: 4px 12px; background: #fff; border-bottom: 1px solid #eee; display: flex; gap: 16px; align-items: center; flex-wrap: wrap; font-size: 13px; }
  #crumbs a { color: #1f4e79; cursor: pointer; text-decoration: underline; }
  #legend { display: flex; gap: 12px; align-items: center; color: #444; margin-left: auto; }
  .sw { display: inline-block; width: 22px; height: 14px; vertical-align: middle; margin-right: 4px; border: 1px solid #333; }
  .sw.l2 { background: #9dc3e6; border-radius: 4px; outline: 1px solid #333; outline-offset: 1px; }
  .sw.l1 { background: #f4b183; clip-path: polygon(30% 0, 70% 0, 100% 30%, 100% 70%, 70% 100%, 30% 100%, 0 70%, 0 30%); border: none; }
  .sw.l0 { background: #fff; border-radius: 50%; }
  .sw.ctx { background: #eee; border-style: dashed; }
  .sw.fr { background: yellow; border-radius: 50%; }
  #status { padding: 4px 12px; font-size: 13px; color: #333; border-bottom: 1px solid #eee; background: #fff; min-height: 24px; }
  #status.error { color: #b00020; }
  #view { flex: 1; overflow: auto; cursor: grab; position: relative; }
  #view.dragging { cursor: grabbing; }
  #view svg { display: block; }
  #view g.node.clickable { cursor: pointer; }
  #view g.node.clickable:hover polygon, #view g.node.clickable:hover ellipse, #view g.node.clickable:hover path { stroke: #e65100; stroke-width: 3; }
  #busy { position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%); background: #000b; color: #fff; padding: 12px 20px; border-radius: 6px; display: none; text-align: center; }
  button { padding: 4px 10px; }
</style>
</head>
<body>
<div id="controls">
  <button id="top" title="Show the fully folded graph">Folded graph</button>
  <label title="Leave out top-level vertices that have no edges"><input id="hideIsolated" type="checkbox" checked> hide isolated</label>
  <span class="sep"></span>
  <select id="mode">
    <option value="vertex">Around vertex</option>
    <option value="path">Along path</option>
  </select>
  <textarea id="query" placeholder="Vertex ids / labels / tooltips. Path: tokens separated by spaces, commas or ->"></textarea>
  <label>radius <input id="radius" type="number" min="0" value="1"></label>
  <label>limit <input id="limit" type="number" min="1" value="300"></label>
  <label><input id="highlight" type="checkbox" checked> highlight path</label>
  <button id="draw">Draw</button>
  <span class="sep"></span>
  <label>engine <select id="engine">
    <option>dot</option><option>neato</option><option>fdp</option><option>sfdp</option><option>twopi</option><option>circo</option>
  </select></label>
  <button id="back" disabled>Back</button>
  <button id="zoomIn">+</button>
  <button id="zoomOut">&minus;</button>
  <button id="fit">Fit</button>
  <button id="dlSvg" disabled>SVG</button>
  <button id="dlDot" disabled>DOT</button>
</div>
<div id="nav">
  <span id="crumbs"></span>
  <span id="legend">
    <span><span class="sw l2"></span>level 2: unbranching path</span>
    <span><span class="sw l1"></span>level 1: local component</span>
    <span><span class="sw l0"></span>original vertex</span>
    <span><span class="sw ctx"></span>context</span>
    <span><span class="sw fr"></span>has hidden neighbours (original view)</span>
  </span>
</div>
<div id="status"></div>
<div id="view"></div>
<div id="busy">Rendering&hellip;<br><small>the full folded graph can take a minute the first time</small></div>
<script>
(function () {
  const $ = id => document.getElementById(id);
  const view = $('view');
  let state = null;          // last successful response, plus .trail (folded views leading to it)
  const history = [];        // previous states for Back
  let scale = 1;
  let baseSize = null;       // natural svg size in px

  function setStatus(text, isError) {
    $('status').textContent = text;
    $('status').className = isError ? 'error' : '';
  }

  async function post(url, params) {
    $('busy').style.display = 'block';
    try {
      const resp = await fetch(url, { method: 'POST', body: new URLSearchParams(params) });
      const text = await resp.text();
      try { return JSON.parse(text); } catch (e) { return { ok: false, message: text }; }
    } catch (e) {
      return { ok: false, message: String(e) };
    } finally {
      $('busy').style.display = 'none';
    }
  }

  function applyScale() {
    const svg = view.querySelector('svg');
    if (!svg || !baseSize) return;
    svg.setAttribute('width', baseSize.w * scale);
    svg.setAttribute('height', baseSize.h * scale);
  }

  function viewName(v) { return v === 'top' ? 'Folded graph' : v; }

  function renderCrumbs() {
    const crumbs = $('crumbs');
    crumbs.innerHTML = '';
    const trail = state.trail || [];
    trail.forEach((v, i) => {
      if (i > 0) crumbs.append(' › ');
      if (i === trail.length - 1 && state.kind === 'fold') {
        const b = document.createElement('b'); b.textContent = viewName(v); crumbs.append(b);
      } else {
        const a = document.createElement('a'); a.textContent = viewName(v);
        a.onclick = () => openFold(v, trail.slice(0, i), null);
        crumbs.append(a);
      }
    });
    if (state.kind === 'orig') {
      if (trail.length) crumbs.append(' › ');
      const b = document.createElement('b');
      b.textContent = state.origLabel ? state.origLabel + ' (original graph)' : 'original graph';
      crumbs.append(b);
    }
  }

  function show(newState, anchor) {
    // anchor: {id, x, y} screen position of an svg node to keep fixed across the re-render
    state = newState;
    view.innerHTML = state.svg.substring(state.svg.indexOf('<svg'));
    const svg = view.querySelector('svg');
    const vb = svg.viewBox.baseVal;
    baseSize = vb && vb.width ? { w: vb.width * 4 / 3, h: vb.height * 4 / 3 } : { w: svg.width.baseVal.value, h: svg.height.baseVal.value };
    applyScale();
    const frontier = new Set(state.frontier);
    for (const g of view.querySelectorAll('g.node')) {
      const kind = g.id.charAt(0), idx = parseInt(g.id.substring(1));
      // Folded views: everything opens something. Original view: only yellow vertices expand.
      if (state.kind === 'fold' || (kind === 'v' && frontier.has(idx))) g.classList.add('clickable');
    }
    if (anchor) {
      const el = document.getElementById(anchor.id);
      if (el) {
        const r = el.getBoundingClientRect();
        view.scrollLeft += r.left + r.width / 2 - anchor.x;
        view.scrollTop += r.top + r.height / 2 - anchor.y;
      }
    }
    let text;
    if (state.kind === 'fold') {
      text = state.message + ' Click a folded vertex to open it, an original vertex to explore it in the original graph.';
    } else {
      text = state.nodes.length + ' vertices drawn, ' + state.frontier.length + ' with hidden neighbours (yellow).';
      if (state.truncated) text += ' Truncated by the vertex limit.';
      if (state.message) text += ' ' + state.message;
    }
    setStatus(text, false);
    renderCrumbs();
    $('back').disabled = history.length === 0;
    $('dlSvg').disabled = $('dlDot').disabled = false;
  }

  function common() {
    return { limit: $('limit').value, engine: $('engine').value };
  }

  function push(res, trail, anchor, resetView) {
    if (state) history.push(state);
    res.trail = trail;
    if (resetView) scale = 1;
    show(res, anchor);
    if (resetView) { view.scrollLeft = 0; view.scrollTop = 0; }
  }

  async function openFold(v, parentTrail, anchor) {
    const res = await post('/api/fold', { view: v, engine: $('engine').value, hide_isolated: $('hideIsolated').checked ? '1' : '0' });
    if (!res.ok) { setStatus(res.message, true); return; }
    push(res, parentTrail.concat([v]), anchor, !anchor);
  }

  async function draw() {
    const res = await post('/api/query', Object.assign(common(), {
      mode: $('mode').value, query: $('query').value, radius: $('radius').value,
      highlight: $('highlight').checked ? '1' : '0'
    }));
    if (!res.ok) { setStatus(res.message, true); return; }
    push(res, [], null, true);
  }

  async function openOriginal(idx) {
    const res = await post('/api/query', Object.assign(common(), { mode: 'index', query: String(idx) }));
    if (!res.ok) { setStatus(res.message, true); return; }
    push(res, state.trail || [], null, true);
  }

  async function openComponent(idx) {
    const res = await post('/api/query', Object.assign(common(), { mode: 'component', query: String(idx) }));
    if (!res.ok) { setStatus(res.message, true); return; }
    res.origLabel = 'C' + idx;
    push(res, state.trail || [], null, true);
  }

  async function expand(idx, anchor) {
    const res = await post('/api/expand', Object.assign(common(), {
      nodes: state.nodes.join(','), highlight: state.highlight.join(','), vertex: idx
    }));
    if (!res.ok) { setStatus(res.message, true); return; }
    res.origLabel = state.origLabel;
    push(res, state.trail, anchor, false);
  }

  $('draw').onclick = draw;
  $('top').onclick = () => openFold('top', [], null);
  $('hideIsolated').onchange = () => { if (state && state.kind === 'fold' && state.view === 'top') openFold('top', [], null); };
  $('query').addEventListener('keydown', e => {
    if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); draw(); }
  });
  $('back').onclick = () => { if (history.length) show(history.pop(), null); };

  function zoomAt(factor, cx, cy) {
    const rect = view.getBoundingClientRect();
    const px = (cx === undefined ? rect.width / 2 : cx - rect.left);
    const py = (cy === undefined ? rect.height / 2 : cy - rect.top);
    const x = view.scrollLeft + px, y = view.scrollTop + py;
    scale = Math.min(20, Math.max(0.01, scale * factor));
    const old = { w: view.scrollWidth, h: view.scrollHeight };
    applyScale();
    view.scrollLeft = x * view.scrollWidth / old.w - px;
    view.scrollTop = y * view.scrollHeight / old.h - py;
  }
  $('zoomIn').onclick = () => zoomAt(1.25);
  $('zoomOut').onclick = () => zoomAt(0.8);
  $('fit').onclick = () => {
    if (!baseSize) return;
    scale = Math.min(view.clientWidth / baseSize.w, view.clientHeight / baseSize.h);
    applyScale();
  };
  view.addEventListener('wheel', e => {
    if (!e.ctrlKey) return;
    e.preventDefault();
    zoomAt(e.deltaY < 0 ? 1.15 : 1 / 1.15, e.clientX, e.clientY);
  }, { passive: false });

  // Drag to pan; a press without movement counts as a click.
  let drag = null;
  view.addEventListener('mousedown', e => {
    if (e.button !== 0) return;
    drag = { x: e.clientX, y: e.clientY, sl: view.scrollLeft, st: view.scrollTop, moved: false };
  });
  window.addEventListener('mousemove', e => {
    if (!drag) return;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    if (Math.abs(dx) + Math.abs(dy) > 4) { drag.moved = true; view.classList.add('dragging'); }
    if (drag.moved) { view.scrollLeft = drag.sl - dx; view.scrollTop = drag.st - dy; }
  });
  window.addEventListener('mouseup', () => { view.classList.remove('dragging'); });
  view.addEventListener('click', e => {
    const moved = drag && drag.moved;
    drag = null;
    if (moved || !state) return;
    const g = e.target.closest('g.node');
    if (!g || !g.classList.contains('clickable')) return;
    e.preventDefault();
    const kind = g.id.charAt(0), idx = parseInt(g.id.substring(1));
    if (kind === 'p') {
      openFold('P' + idx, state.trail || [], null);
    } else if (kind === 'c') {
      openComponent(idx);
    } else if (kind === 'v' && state.kind === 'fold') {
      openOriginal(idx);
    } else if (kind === 'v') {
      const r = g.getBoundingClientRect();
      expand(idx, { id: g.id, x: r.left + r.width / 2, y: r.top + r.height / 2 });
    }
  });

  function download(name, text, type) {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([text], { type }));
    a.download = name;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }
  $('dlSvg').onclick = () => state && download('subgraph.svg', state.svg, 'image/svg+xml');
  $('dlDot').onclick = () => state && download('subgraph.dot', state.dot, 'text/vnd.graphviz');

  openFold('top', [], null);
})();
</script>
</body>
</html>
)PAGE";
}
