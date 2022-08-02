project(nano)

add_executable(nano nano.cpp ReadsAligner.cpp SGraphBuilder.cpp GraphSimplificator.cpp VertexExtender.cpp VertexExtender.h TipResolver.cpp TipResolver.h)
target_link_libraries(nano lja_ec lja_dbg lja_common lja_sequence edlib m)

add_library(ont_graph_simplification nano.cpp ReadsAligner.cpp SGraphBuilder.cpp GraphSimplificator.cpp VertexExtender.cpp VertexExtender.h TipResolver.cpp TipResolver.h)
target_link_libraries(ont_graph_simplification lja_ec lja_dbg lja_common lja_sequence edlib m)