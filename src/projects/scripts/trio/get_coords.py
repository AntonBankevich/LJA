#!/usr/bin/python3
import sys
import os
import random
import logging
from Bio.Seq import Seq
#import networkx as nx
inf = 1e8


# TODO: add logging
# https://github.com/ablab/IsoQuant/blob/master/isoquant.py#L248
# https://docs.python.org/3/howto/logging.html


# for each edge 4 vertices

class Vertex:
    def __init__(self, ver_id):
        # vertices are segments starts & ends, edges are either links(empty seq) or segments (stored in segments ends)
        self.seq = ""
        self.outgoing = []
        self.incoming = []
        self.getInnerId() = ver_id
        self.k = 0
        self.rc_id = (self.getInnerId() // 4) * 4 + 3 - (self.getInnerId() % 4)

    def __str__(self):
        return (f'getInnerId() {self.getInnerId()}, out: {self.outgoing}, inc: {self.incoming}')

    # starts, end, rc_start, rc_end   0,1,2,3, rc_id = 3 - id

    # loops and rc loops
    def is_compressible(self):
        if (len(self.incoming) == 1 and len(self.outgoing) == 1):
            if self.incoming[0] // 2 != self.outgoing[0] // 2:
                return True
        return False

    def indegree(self):
        return len(self.incoming)

    def outdegree(self):
        return len(self.outgoing)

class Edge:
    def __init__(self, edge_id, start_vertex, end_vertex, seq):
        self.start_vertex = start_vertex
        self.edge_id = edge_id
        self.end_vertex = end_vertex
        #temporary for local runs
        self.seq = ""
        self.seq = seq
        self.label = str(self.get_external_id()) + self.get_orientation()

    def get_orientation(self):
        if (self.edge_id % 2 == 0):
            return "+"
        else:
            return "-"

    def get_external_id(self):
        return self.edge_id // 2

    # we really use it, is it bad for RC?
    def get_rc_id(self):
        return (self.edge_id // 2) * 2 + (self.edge_id % 2 + 1) % 2

    def __str__(self):
        return (f'getInnerId() {self.edge_id}, getStart: {self.start_vertex}, end: {self.end_vertex}')

    def length(self):
        return len(self.seq)


class node_stat:
    def __init__(self, length, cov, seq):
        self.length = int(length)
        self.cov = float(cov)
        self.seq = seq

class HaplotypeStats:
    def __init__(self, triobin_str):
        #32      m       0       273     28      390     22      24      112906  17
        #s->seq[i].name, type, s->cnt[i].sc[0], s->cnt[i].sc[1],c[0<<2|2], c[2<<2|0], c[0<<2|1], c[1<<2|0], s->cnt[i].nk, c[0])
        arr = triobin_str.split()
        self.haplotype  = arr[1]
        self.getInnerId() = int(arr[0])
        self.decisive_strips = [int(arr[2]), int(arr[3])]
        self.decisive_counts = [int(arr[4]), int(arr[5])]
        self.total_kmers = int(arr[8])
        #TODO: arr 6, 7, 9
    def is_undefined(self):
        return self.haplotype != "m" and self.haplotype != "p"

# vertices with variable k, edges: gfa_id-> gfa_id*2, gfa_id*2 +1
class Graph:
    def __init__(self, vertices, edges):
        self.vertices = vertices
        self.edges = edges
        self.edge_next_id = max(edges.keys()) + 1

    def indegree(self, vid):
        return len(self.vertices[vid].incoming)

    def outdegree(self, vid):
        return len(self.vertices[vid].outgoing)

    def IsTip(self, eid):
        start_vid = self.edges[eid].start_vertex
        end_vid = self.edges[eid].end_vertex
        if self.indegree(start_vid) == 0 or self.outdegree(end_vid) == 0:
            return True
        else:
            return False

    #simplified check, works only for trivial cases
    def IsBridge(self, eid):
        start_vid = self.edges[eid].start_vertex
        end_vid = self.edges[eid].end_vertex
        bridge = not self.IsTip(eid)
        for alt in self.vertices[start_vid].outgoing:
            if alt != eid and not self.IsTip(alt):
                bridge = False

        for alt in self.vertices[end_vid].incoming:
            if alt != eid and not self.IsTip(alt):
                bridge = False

        return bridge

    def AddEdge(self, start_vertex, end_vertex, new_seq, new_label):
        eid = self.edge_next_id
        self.edges[eid] = Edge(eid, start_vertex, end_vertex, new_seq)
        self.edges[eid].label = new_label
        self.vertices[start_vertex].outgoing.append(eid)
        self.vertices[end_vertex].incoming.append(eid)
        self.edge_next_id += 1

    def InternalRemoveEdge (self, eid):
        incident = [self.edges[eid].start_vertex, self.edges[eid].end_vertex]
        for vid in incident:
            if eid in self.vertices[vid].incoming:
                self.vertices[vid].incoming.remove(eid)
            if eid in self.vertices[vid].outgoing:
                self.vertices[vid].outgoing.remove(eid)
        del self.edges[eid]

    def compress_vertex_if_needed(self, vid):
        if vid in self.vertices.keys():
            if self.indegree(vid) == 1 and self.outdegree(vid) == 1:
                # Need to fix for reverse complement
                if vid != self.vertices[vid].rc_id:
                    work_vs = [self.vertices[vid], self.vertices[self.vertices[vid].rc_id]]
                    to_remove = set()
                    for work_v in work_vs:
                        # dirty -first edge is not changed so we use correct edge even after added RC edge
                        start_e = self.edges[work_v.incoming[0]]
                        end_e = self.edges[work_v.outgoing[0]]
                        to_remove.add(work_v.incoming[0])
                        to_remove.add(work_v.outgoing[0])
                        k = work_v.k
                        if k >= len(start_e.seq):
                            print(f'wrong getEdge length {work_v.incoming[0]}')
                            return
                        new_seq = start_e.seq[:-k] + end_e.seq
                        new_label = start_e.label + end_e.label
#                        print (f'{new_label} {start_e.label} {end_e.label}')
                      
                        self.AddEdge(start_e.start_vertex, end_e.end_vertex, new_seq, new_label)
                    for eid in to_remove:
                        self.InternalRemoveEdge(eid)

    def remove_edge_gfa_id(self, gfa_id):
        to_compress = set()
        for eid in range(gfa_id * 2, gfa_id * 2 + 2):
            if eid in self.edges.keys():
                to_compress.add(self.edges[eid].start_vertex)
                to_compress.add(self.edges[eid].end_vertex)
                # special function, to get used in vertices
                self.InternalRemoveEdge(eid)

        for vid in to_compress:
            self.compress_vertex_if_needed(vid)
    
    def get_internal_length(self, eid):
        getEdge = self.edges[eid]
        res = len(getEdge.seq) - self.vertices[getEdge.start_vertex].k -self.vertices[getEdge.end_vertex].k
        return res

    def print_to_fasta(self, outfile):
        used = set()
        out_f = open(outfile, 'w')
        for e in self.edges.keys():
            if not(self.edges[e].get_external_id() in used):
                out_f.write(">" + self.edges[e].label + '\n')
                out_f.write(self.edges[e].seq + '\n')
                used.add(self.edges[e].get_external_id())

    def get_coords(self, line):
        arr = line.split('_')
        getStart = 0
        for eid in arr:
            getEdge =  self.edges[int(eid) * 2]
            next = getStart +  len(getEdge.seq) - self.vertices[getEdge.end_vertex].k
            print(f'{eid} : {getStart} {next}')
            getStart = next

    def print_to_gfa(self, outfile):
        labels = {}
        canonic = {}
        out_f = open(outfile, 'w')

        for e in self.edges.keys():
            if not(e in labels):
                dirty_label = self.edges[e].label.replace("+", "P").replace("-", "M")
                out_f.write(f'S\t{dirty_label}\t{self.edges[e].seq}\n')
                labels[e] = dirty_label
                canonic[e] = True
                labels[self.edges[e].get_rc_id()] = dirty_label
                canonic[self.edges[e].get_rc_id()] = False
        for vid in self.vertices.keys():
            if self.vertices[vid].rc_id < vid:
                continue
            for eid in self.vertices[vid].incoming:
                inc_label = "+"
                if canonic[eid] == False:
                    inc_label = "-"
                for out_eid in self.vertices[vid].outgoing:
                    out_label = "+"
                    if canonic[out_eid] == False:
                        out_label = "-"
                    out_f.write(f'L\t{labels[eid]}\t{inc_label}\t{labels[out_eid]}\t{out_label}\t{self.vertices[vid].k}M\n')

    def primitive_clean(self):
        changed = True
        max_deletion_length = 1000000
        tips = 0
        bulges = 0
        while changed:
            changed = False
            eids = list(self.edges.keys())
            for eid in eids:
                if not(eid in self.edges.keys()):
                    continue
                v_start_id = self.edges[eid].start_vertex
                v_start = self.vertices[v_start_id]
                if v_start.outdegree() == 2:
                    tid = v_start.outgoing[0]
                    bid = v_start.outgoing[1]
                    if self.edges[tid].end_vertex == self.edges[bid].end_vertex and self.edges[tid].get_external_id()!= self.edges[bid].get_external_id():
                        self.remove_edge_gfa_id(self.edges[tid].get_external_id())
                        changed = True
                        bulges += 1
                elif v_start.indegree() == 0 and v_start.outdegree() == 1:
                    v_end_id = self.edges[eid].end_vertex
                    v_end = self.vertices[v_end_id]
                    eid = v_start.outgoing[0]
                    if v_end.indegree() >= 2 and v_end.outdegree() >= 1 and self.edges[eid].length() < max_deletion_length:
                        self.remove_edge_gfa_id(self.edges[eid].get_external_id())
                        changed = True
                        tips += 1
        print (f'Removed {tips} tips and {bulges} bulges in primitive simplification')

    def print_to_dot(self, outfile, colors):
        dot_graph = nx.MultiDiGraph()
        vertices_to_dot = {}
        count = 0
        for vid in self.vertices.keys():
            dot_graph.add_node(vid, label=self.vertices[vid].k)
        for eid in self.edges.keys():
            e_color = "black"
            ext_id = self.edges[eid].get_external_id()
            if ext_id in colors.keys():
                e_color = colors[ext_id]
            dot_graph.add_edge(self.edges[eid].start_vertex, self.edges[eid].end_vertex, label=self.edges[eid].label, color=e_color)
        pos = nx.nx_agraph.graphviz_layout(dot_graph)
        nx.draw(dot_graph, pos=pos)
        nx.draw(dot_graph)
        nx.drawing.nx_agraph.write_dot(dot_graph, outfile)


def rc(seq):
    seqS = Seq(seq)

    return seqS.reverse_complement().__str__()


'''
    res = seq
    l = len(seq)
    table = {}
    table['A']= 'G'
    table['G'] = 'A'
    table['C'] = 'T'
    table['T'] = 'C'
    table['N'] = 'N'
    for i in range (0, l):
        if not seq[i] in table:
            print (f'{seq} {i} strange symbol {seq[i]}')
            exit(0)
        res[l - 1 - i] = table[seq[i]]
    return res
'''


def get_ids(link_name):
    arr = link_name.split()
    res = [arr[1], arr[3]]
    return res


def get_small_component(getInnerId(), min_size, max_size, min_cov, neighbours, global_used, segments):
    used = set()
    if getInnerId() in global_used:
        return used
    in_job = [getInnerId()]
    while in_job:
        cur_v = in_job.pop()
        if cur_v not in used:
            used.add(cur_v)
            for v in neighbours[cur_v]:
                if v not in used:
                    in_job.append(v)
    #    print (id + " " + str(len(used)))
    global_used.update(used)
    total_len = 0
    total_cov = 0.0
    for f in used:
        total_len += segments[f].length
        total_cov += segments[f].length * segments[f].cov
    total_cov /= total_len
    #    if total_len >= low_cutoff and total_len <= high_cutoff and total_cov > min_coverage:
    #    if total_len >= ids[id].length and total_len > 1000 and total_cov > 10 and len(used) < 2:
    if total_len < max_size and total_len > min_size and total_cov > min_cov and total_len > segments[getInnerId()].length:
        return used
    else:
        return set()


# params: fastg file, ids file, max component size
# looks for ids that contains in small connected components of size 1< SIZE <= max_cutoff


def get_header_id(header):
    pos = header.find("edge_")
    return header[:pos]

#awful but working, just as me
def construct_graph(edge_component, segments, links):
    vertices = {}
    equivalents = {}
    canonic_ids = {}
    id_count = 0
    edge_count = 0
    edges_to_id = {}
    edges = {}
    for e in edge_component:
        e_str = e
        e = int(e)
        # starts, end, rc_start, rc_end   0,1,2,3, rc_id = 3 - id
        for i in range(0, 4):
            vertices[id_count] = Vertex(id_count)
            equivalents[id_count] = set({id_count})
            id_count += 1
        edges_to_id[e] = edge_count
        vertices[edge_count * 4].outgoing.append(e * 2)
        vertices[edge_count * 4 + 2].outgoing.append(e * 2 + 1)
        vertices[edge_count * 4 + 1].incoming.append(e * 2)
        vertices[edge_count * 4 + 3].incoming.append(e * 2 + 1)
        edges[e * 2] = Edge(e * 2, edge_count * 4, edge_count * 4 + 1, segments[e_str].seq)
        edges[e * 2 + 1] = Edge(e * 2 + 1, edge_count * 4 + 2, edge_count * 4 + 3, rc(segments[e_str].seq))
        edge_count += 1

    # adding edges from lines
    for e in edge_component:
        if e not in links:
            continue
        for l in links[e]:
            arr = l.split()
            # L	edge_43116	-	edge_6653	+	0M

            edge_start = e
            edge_end = arr[3]
            first_edge = int(edge_start) * 2
            second_edge = int(edge_end) * 2
            overlap = int(arr[5][:-1])
            link_start_shift = 1
            if arr[2] == '-':
                first_edge += 1
            link_end_shift = 0
            if arr[4] == '-':
                second_edge += 1

            start_id = edges[first_edge].end_vertex
            end_id = edges[second_edge].start_vertex
            #            print (f'from link {l} adding edge {start_id} {end_id} ')
            #            vertices[start_id].next.append(end_id)
            #            vertices[end_id].prev.append(start_id)
            equivalents[start_id].update(equivalents[end_id])
            for tmp in equivalents[start_id]:
                if tmp != start_id:
                    equivalents[tmp].update(equivalents[start_id])
            vertices[start_id].k = overlap
            vertices[end_id].k = overlap
            # rc_Link
            start_id = vertices[start_id].rc_id
            end_id = vertices[end_id].rc_id
            #            print(f'from link {l} adding edge {end_id} {start_id} ')
            #            vertices[end_id].next.append(start_id)
            #            vertices[start_id].prev.append(end_id)
            equivalents[start_id].update(equivalents[end_id])
            for tmp in equivalents[start_id]:
                if tmp != start_id:
                    equivalents[tmp].update(equivalents[start_id])

            vertices[start_id].k = overlap
            vertices[end_id].k = overlap
    for v in vertices.keys():
        canonic_ids[v] = v
        for w in equivalents[v]:
            if w < canonic_ids[v]:
                canonic_ids[v] = w
    canonic_vertices = {}
    for v in vertices.keys():
        if canonic_ids[v] == v:
            canonic_vertices[v] = Vertex(v)
            canonic_vertices[v].rc_id = canonic_ids[vertices[v].rc_id]
            canonic_vertices[v].k = vertices[v].k
    # canonic ids is smallest vetween all equivalents, so it exists
    for v in vertices.keys():
        #        print(f'{v}  {len(vertices[v].outgoing)} {len(vertices[v].incoming)}')
        for i in vertices[v].outgoing:
            #            print(f'  {i}')
            canonic_vertices[canonic_ids[v]].outgoing.append(i)
        for i in vertices[v].incoming:
            #            print(f'  {i}')
            canonic_vertices[canonic_ids[v]].incoming.append(i)
    #    for v in canonic_vertices.keys():
    #        print(canonic_vertices[v])
    for getInnerId() in edges.keys():
        edges[getInnerId()].start_vertex = canonic_ids[edges[getInnerId()].start_vertex]
        edges[getInnerId()].end_vertex = canonic_ids[edges[getInnerId()].end_vertex]
    #        print(edges[id])
    G = Graph(canonic_vertices, edges)
    return G

#deprecated
def get_unbranching_paths(graph):
    used_edges = set()
    for e_id in graph.edges.keys():
        if graph.edges[e_id].get_external_id() not in used_edges and (
        not graph.vertices[graph.edges[e_id].start_vertex].is_compressible()):
            label = ""
            current = e_id
            while True:
                label += str(graph.edges[current].get_external_id()) + "_" + graph.edges[current].get_orientation() + "_"
                cur_v = graph.edges[current].end_vertex
                used_edges.add(graph.edges[current].get_external_id())
                if len(graph.vertices[cur_v].outgoing) != 1:
                    break
                if len(graph.vertices[cur_v].incoming) != 1:
                    break
                current = graph.vertices[cur_v].outgoing[0]
                if graph.edges[current].get_external_id():
                    break
            print(f'PATH {label}')


def get_bulges(graph):
    used = set()
    bulges = {}
    b_file = open("bulges.txt", "w")
    for vid in graph.vertices.keys():
        v = graph.vertices[vid]
        if len(v.outgoing) == 2:
            e1 = graph.edges[v.outgoing[0]]
            e2 = graph.edges[v.outgoing[1]]
            if e1.end_vertex == e2.end_vertex:
                if not (e1.get_external_id() in used) and (not e2.get_external_id() in used):
                    used.add(e1.get_external_id())
                    used.add(e2.get_external_id())
                    l1 = e1.length()
                    l2 = e2.length()
                    v_len = v.k + graph.vertices[e2.end_vertex].k
                    l1 -= v_len
                    l2 -= v_len
                    b_file.write(f'{e1.get_external_id()} {e2.get_external_id()} {l1} {l2} \n')
                    bulges[e1] = e2
                    bulges[e2] = e1
    return bulges

#update only fixable bulges
def update_fixable_haplotypes(bulges, haplotypes):
    count = 0
    for e in bulges.keys():
        h = e.get_external_id()
        if haplotypes[h].is_undefined():
            if haplotypes[bulges[e].get_external_id()].haplotype == "m":
                haplotypes[h].haplotype = "p"
                count += 1
            elif haplotypes[bulges[e].get_external_id()].haplotype == "p":
                haplotypes[h].haplotype = "m"
                count += 1
    print (f'Updated {count} fixable bulges')

    return count


def assign_ambiguous_haplotypes(bulges, haplotypes, graph):
    count = 0
    short_length = 100
    for e in bulges.keys():
        h = e.get_external_id()
        #should be any difference here?
        if haplotypes[h].is_undefined():
            if haplotypes[bulges[e].get_external_id()].is_undefined():
                l1 = graph.get_internal_length(e.edge_id)
                l2 = graph.get_internal_length(bulges[e].edge_id)
                top = haplotypes[h]
                bottom = haplotypes[bulges[e].get_external_id()]
                paternal_count = top.decisive_counts[0] * bottom.decisive_counts[1]
                maternal_count = top.decisive_counts[1] * bottom.decisive_counts[0]
                top_total = top.decisive_counts[0] + top.decisive_counts[1]
                bottom_total = bottom.decisive_counts[0] + bottom.decisive_counts[1]
                decision = 'a'

                if maternal_count > paternal_count or (top_total == 0 and bottom.decisive_counts[0] > bottom.decisive_counts[1]) or (bottom_total == 0 and top.decisive_counts[1] > top.decisive_counts[0]):
                    decision = "m"
                elif paternal_count > maternal_count or (top_total == 0 and bottom.decisive_counts[1] > bottom.decisive_counts[0]) or (bottom_total == 0 and top.decisive_counts[0] > top.decisive_counts[1]):
                    decision = "p"
                else:
                    if l1 < short_length and l2 < short_length:
                        if random.randint(0, 1) == 0:
                            decision = "m"
                        else:
                            decision = "p"
                if decision == "m":
                    haplotypes[h].haplotype = "m"
                    haplotypes[bulges[e].get_external_id()].haplotype = "p"
                    count +=1
                elif decision == "p":
                    haplotypes[h].haplotype = "p"
                    haplotypes[bulges[e].get_external_id()].haplotype = "m"
                    count +=1
                    #others left unfixed. Looks like this is case

    print (f'Updated {count} unfixable short bulges')

def remove_haplotype(haplotypes, graph, to_remove):
    removed = 0
    bridges = 0
    changed = True

    while changed:
        changed = False
        for eid in list(graph.edges):
            if eid not in graph.edges:
                continue
            f = graph.edges[eid].get_external_id()

            if f in haplotypes:
                if haplotypes[f].haplotype == to_remove:
                    if graph.IsBridge(eid):
                        bridges += 1
                        print(f'Bridge {eid}')
                        continue
                    graph.remove_edge_gfa_id(f)
                    removed += 1
                    changed = True
    print (f'Skipped {bridges} bridges')
    print (f'Removed {removed} paternal edges')

    return removed

def get_start_end_vertex(edge_component, segments, edges_to_id):
    max_l = 0
    max_e = ""
    for e in edge_component:
        print(f'getEdge {e} length {segments[e].length}')
        if segments[e].length > max_l:
            max_l = segments[e].length
            max_e = e
    print(max_e)
    return [edges_to_id[max_e] * 4 + 1, edges_to_id[max_e] * 4]

def construct_hashes(seq):
    k = 31
    res = set()
    base = 7
    mod = 123492373
    cur_base = 1
    hash = ord(seq[0])
    for i in range(1, k):
        hash = (hash * base + ord(seq[i])) % mod
        cur_base = (base * cur_base) %mod
    shift_pow = (base * cur_base) %mod
    for i in range(k, len(seq)):
        res.add(hash)
        hash *= base
        hash %= mod
        hash = (hash + mod - ((ord(seq[i - k]) * shift_pow) % mod)) % mod
        hash = (hash + ord(seq[i])) % mod
    res.add(hash)
    return res

def get_distinctive_counts(seq_a, seq_b):
    set_a = construct_hashes(seq_a)
    set_b = construct_hashes(seq_b)
    dist_a = 0
    dist_b = 0
    for s in set_a:
        if not (s in set_b):
            dist_a +=1
    for s in set_b:
        if not (s in set_a):
            dist_b += 1

    return [dist_a, dist_b]

def print_table(bulges, haplotypes, graph):
    used = set()
    print(f'topid\tbottomid\tstart_v.k\tend_v.k\ttop_len\tbottom_len\ttop_class\tbottom_class\ttop_p_strip\ttop_m_strip\tbottom_p_strip\tbottom_m_strip\ttop_p_count\ttop_m_count\tbottom_p_count\tbottom_m_count\tdistinct_top\tdistinct_bottom\t')
    for top in bulges.keys():
        topid = top.get_external_id()
        if topid in used:
            continue
        bottom = bulges[top]
        bulge_dist = get_distinctive_counts(top.seq, bottom.seq)
        bid =bottom.get_external_id()
        used.add(topid)
        used.add(bid)
        table_str = (f'{topid}\t{bid}\t' +
                     f'{graph.vertices[top.start_vertex].k}\t{graph.vertices[top.end_vertex].k}\t{top.length() -graph.vertices[top.start_vertex].k - graph.vertices[top.end_vertex].k }\t{bottom.length()  -graph.vertices[top.start_vertex].k - graph.vertices[top.end_vertex].k}\t' +
                     f'{haplotypes[topid].haplotype}\t{haplotypes[bid].haplotype}\t' +
                     f'{haplotypes[topid].decisive_strips[0]}\t{haplotypes[topid].decisive_strips[1]}\t' +
                     f'{haplotypes[bid].decisive_strips[0]}\t{haplotypes[bid].decisive_strips[1]}\t' +
                     f'{haplotypes[topid].decisive_counts[0]}\t{haplotypes[topid].decisive_counts[1]}\t' +
                     f'{haplotypes[bid].decisive_counts[0]}\t{haplotypes[bid].decisive_counts[1]}\t' +
                     f'{bulge_dist[0]}\t{bulge_dist[1]}')
        print (table_str)

def run_extraction(graph_f, haplotypes_f, upd_line):
    neighbours = {}
    segments = {}
    links = {}

    pred = {}
    good = set()
    for line in open(graph_f, 'r'):
        if line[0] == "L":
            arr = get_ids(line)
            if len(arr) <= 1:
                print(line)
                exit()
            if not (arr[0] in neighbours.keys()):
                neighbours[arr[0]] = set()
            if not (arr[1] in neighbours.keys()):
                neighbours[arr[1]] = set()
            neighbours[arr[0]].add(arr[1])
            neighbours[arr[1]].add(arr[0])
            if not arr[0] in links:
                links[arr[0]] = []
            if not arr[1] in links:
                links[arr[1]] = []
            links[arr[0]].append(line)

        elif line[0] == "S":
            arr = line.split()
            length = len(arr[2])
            #            cov = arr[3].split(':')[2]
            cov = 1
            segments[arr[1]] = node_stat(length, cov, arr[2])
            neighbours[arr[1]] = set()

    unique = set()
    total = 0
    print("Constructing graph...")

    graph = construct_graph(segments.keys(), segments, links)
    graph.print_to_gfa("check.gfa")
    graph.get_coords(upd_line)
    exit(0)
#    graph.print_to_dot("tst.dot", {})
#    graph.print_to_dot("tst.dot", {})
    bulges = get_bulges(graph)
#    print_table(bulges, haplotypes, graph)
    update_fixable_haplotypes(bulges, haplotypes)
    assign_ambiguous_haplotypes(bulges, haplotypes, graph)
    remove_haplotype(haplotypes,graph,"m")
    '''
    removed = 0
    for f in haplotypes.keys():
        #TODO parameter
        if haplotypes[f].haplotype == "m":
            graph.remove_edge_gfa_id(int(f))
            removed += 1
    print (f'Removed {removed} paternal edges')
    '''
    graph.print_to_fasta("maternal.fasta")
    graph.print_to_gfa("maternal.gfa")

    graph.primitive_clean()
    graph.print_to_fasta("additional_cleaning.fasta")
    graph.print_to_gfa("additional_cleaning.gfa")
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f'Script for gfa compression after removal of one haplotype')
        print(f'Usage: {sys.argv[0]} <graph.gfa> <haplotype.bin>')
        exit(0)
    graph = sys.argv[1]
 

    haplotypes_list = sys.argv[2]
    line = sys.argv[3].strip()
    random.seed(239)
    run_extraction(graph, haplotypes_list, line)

# print (total)
# for f in unique:
#    print (f)
