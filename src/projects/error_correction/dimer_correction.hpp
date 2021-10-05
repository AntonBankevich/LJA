#pragma once
#include "dbg/graph_modification.hpp"
#include "dbg/compact_path.hpp"

GraphAlignment correctFromStart(const GraphAlignment &al, double reliable_coverage);
size_t CorrectDimers(logging::Logger &logger, RecordStorage &reads_storage, size_t k, size_t threads, double reliable_coverage);