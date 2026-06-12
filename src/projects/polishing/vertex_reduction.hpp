#pragma once
#include <unordered_map>
#include "assembly_graph/assembly_graph_base.hpp"

std::unordered_map<ag::VertexId, Segment<ag::Vertex>> ConstructReduction(ag::AssemblyGraph &graph, size_t min_overlap, size_t max_repeat);
