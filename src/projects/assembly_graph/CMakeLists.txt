project(assembly_graph)

add_library(lja_ag STATIC assembly_graph_base.cpp assembly_graph.cpp graph_paths.cpp vertex_resolution.cpp edge_code_listener.cpp data_structures/read_alignment_storage.cpp data_structures/suffix_tracker.cpp random_access_paths.cpp graph_listeners.cpp alibned_read.cpp data_structures/read_logger.cpp data_structures/aligned_read_listeners.cpp data_structures/component.cpp ag_algorithms.cpp data_structures/splitters.cpp
        data_structures/aligned_read_statistics_tracker.cpp
        data_structures/coverage_samples.cpp
        dll_path.cpp
        dll_path_storage.cpp)
target_link_libraries (lja_ag m ${OpenMP_CXX_FLAGS} lja_sequence lja_common lja_alignment stdc++fs Threads::Threads)
