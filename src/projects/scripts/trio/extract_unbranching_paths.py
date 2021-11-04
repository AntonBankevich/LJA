#!/usr/bin/python3
import sys
import os
import random
import logging
from Bio.Seq import Seq
inf = 1e8
#TODO: add logging
# https://github.com/ablab/IsoQuant/blob/master/isoquant.py#L248
# https://docs.python.org/3/howto/logging.html


#for each edge 4 vertices

class Vertex:
    def __init__ (self, ver_id):
#vertices are segments starts & ends, edges are either links(empty seq) or segments (stored in segments ends)
        self.seq = ""
        self.outgoing = []
        self.incoming = []
        self.id = ver_id
        self.size = 0
        self.rc_id = (self.id // 4) * 4 + 3 - (self.id % 4)
    def __str__(self):
        return (f'id {self.id}, out: {self.outgoing}, inc: {self.incoming}')
#starts, end, rc_start, rc_end   0,1,2,3, rc_id = 3 - id

#loops and rc loops
    def is_compressible(self):
        if (len(self.incoming) == 1 and len(self.outgoing) == 1):
            if self.incoming[0]//2 != self.outgoing[0] //2:
                return True
        return False

class Edge:
    def __init__ (self, edge_id, start_vertex, end_vertex, seq):
        self.start_vertex = start_vertex
        self.edge_id = edge_id
        self.end_vertex = end_vertex
        self.seq = seq

    def get_orientation(self):
        if (self.edge_id %2 == 0):
            return "+"
        else:
            return "-"
    def get_external_id (self):
        return self.edge_id // 2
    def __str__(self):
        return (f'id {self.edge_id}, start: {self.start_vertex}, end: {self.end_vertex}')

    def length(self):
        return len(self.seq)
class node_stat:
    def __init__ (self, length, cov, seq):   
        self.length = int(length)
        self.cov = float (cov)
        self.seq = seq

#class Graph:
    
def rc (seq):

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


def get_small_component(id, min_size, max_size, min_cov, neighbours, global_used, segments):
    used = set()
    if id in global_used:
        return used
    in_job = [id]
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
    if total_len < max_size and total_len > min_size and total_cov > min_cov and total_len > segments[id].length:
        return used
    else:
        return set()
#params: fastg file, ids file, max component size
#looks for ids that contains in small connected components of size 1< SIZE <= max_cutoff 


def get_header_id(header):
    pos = header.find("edge_")
    return header[:pos]

def construct_graph(edge_component, segments, links, forbidden):
    vertices = {}
    equivalents = {}
    canonic_ids = {}
    id_count = 0
    edge_count = 0
    edges_to_id = {}
    edges = {}
    for e in edge_component:
        if e in forbidden:
            continue
        e_str = e
        e = int(e)
        #starts, end, rc_start, rc_end   0,1,2,3, rc_id = 3 - id
        for i in range(0, 4):
            vertices[id_count] = Vertex(id_count)
            equivalents[id_count] = set({id_count})
            id_count += 1
        edges_to_id[e] = edge_count
        vertices[edge_count * 4].outgoing.append(e * 2)
        vertices[edge_count * 4 + 2].outgoing.append(e * 2+1)
        vertices[edge_count * 4 +1 ].incoming.append(e * 2)
        vertices[edge_count * 4 +3 ].incoming.append(e * 2+1)
        edges[e*2] = Edge(e*2, edge_count * 4, edge_count * 4 + 1, segments[e_str].seq)
        edges[e*2+1] = Edge(e*2+1, edge_count * 4+2, edge_count * 4 + 3, rc(segments[e_str].seq))
        edge_count += 1
        
#adding edges from lines
    for e in edge_component:
        if e not in links:
            continue
        if e in forbidden:
            continue
        for l in links[e]:
            arr = l.split()
    #L	edge_43116	-	edge_6653	+	0M
           
            edge_start = e
            edge_end = arr[3]
            if edge_end in forbidden:
                continue
            first_edge = int(edge_start) * 2
            second_edge = int(edge_end) * 2
            overlap = int(arr[5][:-1])
            link_start_shift = 1
            if arr[2] == '-':
                first_edge+=1
            link_end_shift = 0
            if arr[4] == '-':
                second_edge +=1


            start_id = edges[first_edge].end_vertex
            end_id = edges[second_edge].start_vertex
#            print (f'from link {l} adding edge {start_id} {end_id} ')
#            vertices[start_id].next.append(end_id)
#            vertices[end_id].prev.append(start_id)
            equivalents[start_id].update(equivalents[end_id])
            for tmp in equivalents[start_id]:
                if tmp != start_id:
                    equivalents[tmp].update(equivalents[start_id])
            vertices[start_id].size = overlap
            vertices[end_id].size = overlap         
#rc_Link
            start_id = vertices[start_id].rc_id
            end_id = vertices[end_id].rc_id
#            print(f'from link {l} adding edge {end_id} {start_id} ')
#            vertices[end_id].next.append(start_id)
#            vertices[start_id].prev.append(end_id)
            equivalents[start_id].update(equivalents[end_id])
            for tmp in equivalents[start_id]:
                if tmp != start_id:
                    equivalents[tmp].update(equivalents[start_id])

            vertices[start_id].size = overlap
            vertices[end_id].size = overlap         
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
            canonic_vertices[v].size = vertices[v].size
#canonic ids is smallest vetween all equivalents, so it exists
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
    for id in edges.keys():
        
        edges[id].start_vertex = canonic_ids[edges[id].start_vertex]
        edges[id].end_vertex = canonic_ids[edges[id].end_vertex]
#        print(edges[id])
    return [canonic_vertices, edges]

def get_unbranching_paths(vertices, edges):
    used_edges = set()
    for e_id in edges.keys():
        if edges[e_id].get_external_id() not in used_edges and (not vertices[edges[e_id].start_vertex].is_compressible()):
            label = ""
            current = e_id
            while True:
                label += str(edges[current].get_external_id()) + "_" + edges[current].get_orientation() + "_"
                cur_v = edges[current].end_vertex
                used_edges.add(edges[current].get_external_id())
                if len(vertices[cur_v].outgoing) != 1:
                    break
                if len(vertices[cur_v].incoming) != 1:
                    break
                current = vertices[cur_v].outgoing[0]
                if edges[current].get_external_id():
                    break
            print (f'PATH {label}')

def print_bulges(vertices, edges):
    used = set()
    for vid in vertices.keys():
        v = vertices[vid]
        if len(v.outgoing) == 2:

            e1 = edges[v.outgoing[0]]
            e2 = edges[v.outgoing[1]]
            if e1.end_vertex == e2.end_vertex:
                if not (e1.get_external_id() in used) and (not e2.get_external_id() in used):
                    used.add(e1.get_external_id())
                    used.add(e2.get_external_id())
                    l1 = e1.length()
                    l2 = e2.length()
                    v_len = v.size + vertices[e2.end_vertex].size
                    l1 -= v_len
                    l2 -= v_len
                    print (f'{e1.get_external_id()} {e2.get_external_id()} {l1} {l2} ')

def get_bugles():
    return 0



def get_start_end_vertex(edge_component, segments, edges_to_id):
    max_l = 0
    max_e = ""
    for e in edge_component:
        print (f'edge {e} length {segments[e].length}')
        if segments[e].length > max_l:
            max_l = segments[e].length
            max_e = e
    print (max_e)
    return [edges_to_id[max_e] * 4 + 1, edges_to_id[max_e] * 4]



def run_extraction(graph_f, forbidden_f):
    neighbours = {}
    segments = {}
    links = {}

    pred = {}
    good = set()
    for line in open(graph_f, 'r'):
        if line[0] == "L":
            arr = get_ids(line)
            if len(arr) == 0:
                print (line)
                exit()
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
    forbidden = set()
    for line in open (forbidden_f, 'r'):
        forbidden.add (line.strip())
    print("Constructing graph...")

    [vertices, edges] = construct_graph(segments.keys(), segments, links, forbidden)
    get_unbranching_paths(vertices, edges)
    print_bulges(vertices, edges)
    print ("Constructed")
    
if __name__ == "__main__":
    if len (sys.argv) != 3:
        print(f'Script for gfa compression after removal of forbidden edges')
        print(f'Usage: {sys.argv[0]} <graph.gfa> <forbidden.list>')
        exit(0)
    graph = sys.argv[1]
    forbidden_list = sys.argv[2]
    run_extraction(graph, forbidden_list)
#print (total)
#for f in unique:
#    print (f)
