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
  #controls textarea { width: 420px; height: 34px; font: 13px monospace; resize: vertical; }
  #controls input[type=number] { width: 70px; }
  #controls label { display: flex; align-items: center; gap: 4px; }
  #status { padding: 4px 12px; font-size: 13px; color: #333; border-bottom: 1px solid #eee; background: #fff; min-height: 24px; }
  #status.error { color: #b00020; }
  #view { flex: 1; overflow: auto; cursor: grab; position: relative; }
  #view.dragging { cursor: grabbing; }
  #view svg { display: block; }
  #view g.node.frontier { cursor: pointer; }
  #view g.node.frontier:hover polygon, #view g.node.frontier:hover ellipse, #view g.node.frontier:hover path { stroke: #e65100; stroke-width: 3; }
  #busy { position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%); background: #000a; color: #fff; padding: 10px 18px; border-radius: 6px; display: none; }
  button { padding: 4px 10px; }
</style>
</head>
<body>
<div id="controls">
  <select id="mode">
    <option value="vertex">Around vertex</option>
    <option value="path">Along path</option>
  </select>
  <textarea id="query" placeholder="Vertex ids / labels / tooltips. Path: tokens separated by spaces, commas or ->"></textarea>
  <label>radius <input id="radius" type="number" min="0" value="1"></label>
  <label>limit <input id="limit" type="number" min="1" value="300"></label>
  <label>engine <select id="engine">
    <option>dot</option><option>neato</option><option>fdp</option><option>sfdp</option><option>twopi</option><option>circo</option>
  </select></label>
  <label><input id="highlight" type="checkbox" checked> highlight path</label>
  <button id="draw">Draw</button>
  <button id="back" disabled>Back</button>
  <button id="zoomIn">+</button>
  <button id="zoomOut">&minus;</button>
  <button id="fit">Fit</button>
  <button id="dlSvg" disabled>SVG</button>
  <button id="dlDot" disabled>DOT</button>
</div>
<div id="status"></div>
<div id="view"></div>
<div id="busy">Rendering&hellip;</div>
<script>
(function () {
  const $ = id => document.getElementById(id);
  const view = $('view');
  let state = null;          // last successful response
  const history = [];        // previous states for Back
  let scale = 1;
  let baseSize = null;       // natural svg size in px

  function setStatus(text, isError) {
    $('status').textContent = text;
    $('status').className = isError ? 'error' : '';
  }

  fetch('/api/info').then(r => r.json()).then(info => {
    if (!state) setStatus('Graph: ' + info.nodes + ' vertices, ' + info.edges + ' edges. Yellow vertices are clickable: they have neighbours not yet drawn.');
  });

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

  function show(newState, anchor) {
    // anchor: {idx, x, y} screen position of a node to keep fixed across the re-render
    state = newState;
    view.innerHTML = state.svg.substring(state.svg.indexOf('<svg'));
    const svg = view.querySelector('svg');
    const vb = svg.viewBox.baseVal;
    baseSize = vb && vb.width ? { w: vb.width * 4 / 3, h: vb.height * 4 / 3 } : { w: svg.width.baseVal.value, h: svg.height.baseVal.value };
    applyScale();
    const frontier = new Set(state.frontier);
    for (const g of view.querySelectorAll('g.node')) {
      const idx = parseInt(g.id.substring(1));
      if (frontier.has(idx)) g.classList.add('frontier');
    }
    if (anchor) {
      const el = document.getElementById('v' + anchor.idx);
      if (el) {
        const r = el.getBoundingClientRect();
        view.scrollLeft += r.left + r.width / 2 - anchor.x;
        view.scrollTop += r.top + r.height / 2 - anchor.y;
      }
    }
    let text = state.nodes.length + ' vertices drawn, ' + state.frontier.length + ' with hidden neighbours (yellow).';
    if (state.truncated) text += ' Truncated by the vertex limit.';
    if (state.message) text += ' ' + state.message;
    setStatus(text, false);
    $('back').disabled = history.length === 0;
    $('dlSvg').disabled = $('dlDot').disabled = false;
  }

  function common() {
    return { limit: $('limit').value, engine: $('engine').value };
  }

  async function draw() {
    const res = await post('/api/query', Object.assign(common(), {
      mode: $('mode').value, query: $('query').value, radius: $('radius').value,
      highlight: $('highlight').checked ? '1' : '0'
    }));
    if (!res.ok) { setStatus(res.message, true); return; }
    if (state) history.push(state);
    scale = 1;
    show(res, null);
    view.scrollLeft = 0; view.scrollTop = 0;
  }

  async function expand(idx, anchor) {
    const res = await post('/api/expand', Object.assign(common(), {
      nodes: state.nodes.join(','), highlight: state.highlight.join(','), vertex: idx
    }));
    if (!res.ok) { setStatus(res.message, true); return; }
    history.push(state);
    show(res, anchor);
  }

  $('draw').onclick = draw;
  $('query').addEventListener('keydown', e => {
    if (e.key === 'Enter' && !e.shiftKey) { e.preventDefault(); draw(); }
  });
  $('back').onclick = () => { if (history.length) show(history.pop(), null); };

  function zoomAt(factor, cx, cy) {
    const rect = view.getBoundingClientRect();
    const px = (cx === undefined ? rect.width / 2 : cx - rect.left);
    const py = (cy === undefined ? rect.height / 2 : cy - rect.top);
    const x = view.scrollLeft + px, y = view.scrollTop + py;
    scale = Math.min(20, Math.max(0.02, scale * factor));
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
    if (!g || !g.classList.contains('frontier')) return;
    e.preventDefault();
    const r = g.getBoundingClientRect();
    expand(parseInt(g.id.substring(1)), { idx: parseInt(g.id.substring(1)), x: r.left + r.width / 2, y: r.top + r.height / 2 });
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
})();
</script>
</body>
</html>
)PAGE";
}
