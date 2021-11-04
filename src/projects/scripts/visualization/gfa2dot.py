#Copyright by tdvorkina
import networkx as nx

import sys
from pathlib import Path

def get_vertex(e, edges, edges_rev, used, vertex_cnt):
    is_loop = False
    if e in edges:
        for ee in edges[e]:
            if e == ee[0]:
                is_loop = True
                break
        if is_loop:
            in_vertex, in_len = -1, -1
            for ee in edges[e]:
                in_len = ee[1]
                if ee[0] in used:
                    in_vertex = used[ee[0]]["out"]
            out_vertex, out_len = -1, -1
            if e in edges_rev:
                for ee in edges_rev[e]:
                    out_len = ee[1]
                    if ee[0] in used:
                        out_vertex = used[ee[0]]["in"]
            if in_vertex == -1:
                in_vertex, in_len = out_vertex, out_len
            if out_vertex == -1:
                out_vertex, out_len == in_vertex, in_len
            if in_vertex == -1 and out_vertex == -1:
                vertex_cnt += 1
                in_vertex, out_vertex = vertex_cnt, vertex_cnt
                in_len, out_len = edges[e][0][1], edges[e][0][1]
            return in_vertex, in_len, out_vertex, out_len, vertex_cnt

    in_vertex, in_len = -1, -1
    if e in edges:
       for ee in edges[e]:
            in_len = ee[1]
            if ee[0] in used:
                in_vertex = used[ee[0]]["out"]
       if len(edges[e]) > 0 and in_vertex == -1:
            ee = edges[e][0][0]
            for eee in edges_rev[ee]:
                if eee[0] in used:
                    in_vertex = used[eee[0]]["in"]
    if in_vertex == -1:
        vertex_cnt += 1
        in_vertex = vertex_cnt

    out_vertex, out_len = -1, -1
    if e in edges_rev:
        for ee in edges_rev[e]:
            out_len = ee[1]
            if ee[0] in used:
                out_vertex = used[ee[0]]["in"]
        if len(edges_rev[e]) > 0 and out_vertex == -1:
            ee = edges_rev[e][0][0]
            for eee in edges[ee]:
                if eee[0] in used:
                    out_vertex = used[eee[0]]["out"]
    if out_vertex == -1:
         vertex_cnt += 1
         out_vertex = vertex_cnt
    return in_vertex, in_len, out_vertex, out_len, vertex_cnt

gfa = sys.argv[1]
outfilename = Path(sys.argv[2])
edges_color = {}
if len(sys.argv) > 3:
    color_file = Path(sys.argv[3])
    with color_file.open() as f:
        for ln in f.readlines():
            e_id, color = ln.strip().split("\t")
            if ("-" not in e_id) and ("+" not in e_id):
                e_id = e_id + "+"
            elif e_id.startswith("-"):
                e_id = e_id[1:] + "-"
            elif e_id.startswith("+"):
                e_id = e_id[1:] + "+"
            edges_color[e_id] = color

edges = {}
edges_rev = {}
edges_len = {}
with open(gfa, "r") as fin:
    for ln in fin.readlines():
        if ln.startswith("S"):
            _, e_id, e_str = ln.strip().split("\t")[:3]
            edges_len[e_id + "+"] = len(e_str)
            edges_len[e_id + "-"] = len(e_str)
        elif ln.startswith("L"):
            _, e1, or1, e2, or2, v_len = ln.strip().split("\t")
            if e1 + or1 not in edges:
                edges[e1 + or1] = []
            if e2 + or2 not in edges_rev:
                edges_rev[e2 + or2] = []
            edges[e1 + or1].append([e2 + or2, v_len[:-1]])
            edges_rev[e2 + or2].append([e1 + or1, v_len[:-1]])

            or1 = "+" if or1 == "-" else "-"
            or2 = "+" if or2 == "-" else "-"
            if e2 + or2 not in edges:
                edges[e2 + or2] = []
            if e1 + or1 not in edges_rev:
                edges_rev[e1 + or1] = []
            edges[e2 + or2].append([e1 + or1, v_len[:-1]])
            edges_rev[e1 + or1].append([e2 + or2, v_len[:-1]])

graph = nx.MultiDiGraph()
links = []
used = {}
vertex_cnt = 0
vertex_len = [(i, {"label": "-1"}) for i in range(2*len(edges_len))]
for e in edges_len:
    in_vertex, in_len, out_vertex, out_len, vertex_cnt = get_vertex(e, edges, edges_rev, used, vertex_cnt)
    vertex_len[in_vertex] = (in_vertex, {"label": in_len})
    vertex_len[out_vertex] = (out_vertex, {"label": out_len})
    color = "black" if e not in edges_color else edges_color[e]
    links.append([in_vertex, in_len, out_vertex, out_len, e, str(edges_len[e]), color])
    used[e] = {"in": in_vertex, "out": out_vertex}


vertex_len = vertex_len[1:vertex_cnt + 1]
graph.add_nodes_from(vertex_len)

with open(outfilename.with_suffix(".tsv"), "w") as fout:
    for l in links:
        fout.write("\t".join([str(it) for it in l]) + "\n")
        graph.add_edge(l[0], l[2], label = l[5], color = l[6])

pos = nx.nx_agraph.graphviz_layout(graph)
nx.draw(graph, pos=pos)
nx.drawing.nx_agraph.write_dot(graph, outfilename)
