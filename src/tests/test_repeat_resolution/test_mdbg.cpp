//
// Created by Andrey Bzikadze on 1/13/21.
//

#include "repeat_resolution/mdbg.hpp"
#include "repeat_resolution/mdbg_inc.hpp"
#include "repeat_resolution/paths.hpp"
#include "gtest/gtest.h"

using namespace repeat_resolution;

using RawEdgeInfo =
std::vector<std::tuple<RRVertexType, RRVertexType, std::string>>;

std::vector<SuccinctEdgeInfo>
GetEdgeInfo(std::map<RRVertexType, dbg::Vertex> &vertexes,
            std::vector<dbg::Edge> &edges, const RawEdgeInfo &raw_edge_info,
            int k, bool unique) {
    for (const auto &[st, en, str] : raw_edge_info) {
        Sequence seq(str);
        vertexes.emplace(st, st);
        dbg::Vertex &st_v = vertexes.at(st);
        st_v.seq = seq.Prefix(k);

        edges.emplace_back(&st_v, nullptr, seq.Subseq(k));
    }

    std::vector<SuccinctEdgeInfo> edge_info;
    for (auto it = raw_edge_info.begin(); it!=raw_edge_info.end(); ++it) {
        edge_info.push_back({std::get<0>(*it), std::get<1>(*it),
                             &(edges[it - raw_edge_info.begin()]), unique});
    }
    return edge_info;
}

using RawVertexInfo =
std::unordered_map<RRVertexType, std::pair<std::string, bool>>;

std::pair<bool, bool> CompareVertexes(const MultiplexDBG &graph,
                                      const RawVertexInfo &vertex_info) {
    bool IndexSetsEqual = [&graph, &vertex_info]() {
      std::unordered_set<RRVertexType> obs_vertex_set;
      for (const auto &vertex : graph) {
          obs_vertex_set.emplace(vertex);
      }

      std::unordered_set<RRVertexType> true_vertex_set;
      for (const auto &[index, vertex_prop] : vertex_info) {
          true_vertex_set.emplace(index);
      }
      return obs_vertex_set==true_vertex_set;
    }();

    bool PropsEquals = [&graph, &vertex_info]() {
      RawVertexInfo obs_props;
      for (const auto &vertex : graph) {
          const RRVertexProperty &vertex_prop = graph.node_prop(vertex);
          obs_props.emplace(vertex,
                            std::make_pair(vertex_prop.Seq().ToSequence().str(),
                                           vertex_prop.IsFrozen()));
      }

      return obs_props==vertex_info;
    }();
    return {IndexSetsEqual, PropsEquals};
}

bool CompareEdges(const MultiplexDBG &graph, const RawEdgeInfo &edge_info) {
    int cnt = 0;
    for (const auto &vertex : graph) {
        auto[nbr_begin, nbr_end] = graph.out_neighbors(vertex);
        for (auto nbr_it = nbr_begin; nbr_it!=nbr_end; ++nbr_it) {
            MDBGSeq seq =
                graph.GetEdgeSequence(graph.find(vertex), nbr_it, false, false);

            std::tuple<RRVertexType, RRVertexType, std::string> edge{
                vertex, nbr_it->first, seq.ToSequence().str()};
            if (std::find(edge_info.begin(), edge_info.end(), edge)==
                edge_info.end()) {
                std::cout << std::get<0>(edge) << " " << std::get<1>(edge)
                          << " "
                          << std::get<2>(edge) << "\n";
                std::cout
                    << "Found an edge that is not present among true edges";
                return false;
            }
            ++cnt;
        }
    }
    return cnt==edge_info.size();
}

TEST(DB1, Basic) {
    const size_t k = 2;

    const bool frozen = false;
    const bool unique = false;
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    RawEdgeInfo raw_edge_info{{0, 2, "CCT"},  // 0
                              {1, 2, "GACT"}, // 1
                              {2, 3, "CTAG"}, // 2
                              {3, 4, "AGTT"}, // 3
                              {3, 5, "AGC"},  // 4
                              {2, 4, "CTT"}}; // 5
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, unique);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 2, 3}});
      _path_vector.emplace_back(RRPath{"1", std::list<size_t>{1, 5}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);

    RawVertexInfo vertex_info = {{0, {"CC", frozen}}, {1, {"GA", frozen}},
                                 {2, {"CT", frozen}}, {3, {"AG", frozen}},
                                 {4, {"TT", frozen}}, {5, {"GC", frozen}}};
    {
        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, raw_edge_info));
    }
}

TEST(DBSingleEdge1, Basic) {
    const size_t k = 2;

    const bool frozen = false;
    const bool unique = false;
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<std::tuple<RRVertexType, RRVertexType, std::string>>
        raw_edge_info{{0, 1, "ACGTTGCA"}}; // 0
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, unique);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);

    {
        RawVertexInfo
            vertex_info = {{0, {"ACG", frozen}}, {1, {"GCA", frozen}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
            {0, 1, "ACGTTGCA"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, raw_edge_info));
    }
}

TEST(DBSingleEdge2, Basic) {
    const size_t k = 2;

    const bool frozen = false;
    const bool unique = false;
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<std::tuple<RRVertexType, RRVertexType, std::string>>
        raw_edge_info{{0, 1, "ACGCA"}}; // 0
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, unique);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo
            vertex_info = {{0, {"ACG", frozen}}, {1, {"GCA", frozen}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
            {0, 1, "ACGCA"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, raw_edge_info));
    }
}

TEST(DBSingleEdge3, Basic) {
    const size_t k = 2;

    const bool frozen = false;
    const bool unique = false;
    std::vector<std::tuple<RRVertexType, RRVertexType, std::string>>
        raw_edge_info{{0, 1, "ACGTGCA"}}; // 0
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, unique);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    int N = 5;
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + N, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info = {{0, {"ACGTGCA", true}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>>
            post_raw_edge{};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

TEST(DBStVertex, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "AAAAA"}, {0, 2, "AAACA"}, {0, 3, "AAA"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{4, {"AAA", false}},
                                  {1, {"AAA", false}},
                                  {5, {"AAA", false}},
                                  {2, {"ACA", false}},
                                  {6, {"AAA", true}}};

        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {4, 1, "AAAAA"}, {5, 2, "AAACA"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

TEST(DBEvVertex, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 3, "AAAAA"}, {1, 3, "AACAA"}, {2, 3, "AAA"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);

    {
        RawVertexInfo vertex_info{{0, {"AAA", false}},
                                  {1, {"AAC", false}},
                                  {4, {"AAA", false}},
                                  {5, {"CAA", false}},
                                  {2, {"AAA", true}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 4, "AAAAA"}, {1, 5, "AACAA"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph 1-in >1-out
TEST(DB1inVertex, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "AACAG"}, {1, 2, "AGACC"}, {1, 3, "AGATT"}, {1, 4, "AGAGG"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"AAC", false}},
                                  {1, {"CAG", false}},
                                  {2, {"ACC", false}},
                                  {3, {"ATT", false}},
                                  {4, {"AGG", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 1, "AACAG"}, {1, 2, "CAGACC"}, {1, 3, "CAGATT"},
            {1, 4, "CAGAGG"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph 1-in >1-out with 1-in transforming into a vertex
TEST(DB1inVertex, WithShortEdge) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "CAG"}, {1, 2, "AGACC"}, {1, 3, "AGATT"}, {1, 4, "AGAGG"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"CAG", false}},
                                  {2, {"ACC", false}},
                                  {3, {"ATT", false}},
                                  {4, {"AGG", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 2, "CAGACC"}, {0, 3, "CAGATT"}, {0, 4, "CAGAGG"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph >1-in 1-out
TEST(DB1outVertex, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 3, "CCAGA"}, {1, 3, "TTAGA"}, {2, 3, "GGAGA"}, {3, 4, "GAAAA"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"CCA", false}},
                                  {1, {"TTA", false}},
                                  {2, {"GGA", false}},
                                  {3, {"GAA", false}},
                                  {4, {"AAA", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 3, "CCAGAA"}, {1, 3, "TTAGAA"}, {2, 3, "GGAGAA"},
            {3, 4, "GAAAA"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph >1-in 1-out with 1-in transforming into a vertex
TEST(DB1outVertex, WithShortEdge) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 3, "CCAGA"}, {1, 3, "TTAGA"}, {2, 3, "GGAGA"}, {3, 4, "GAA"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"CCA", false}},
                                  {1, {"TTA", false}},
                                  {2, {"GGA", false}},
                                  {3, {"GAA", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 3, "CCAGAA"}, {1, 3, "TTAGAA"}, {2, 3, "GGAGAA"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with a complex vertex (2in-2out)
TEST(DBComplexVertex, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 2, "ACAAA"}, {1, 2, "GGAAA"}, {2, 3, "AATGC"}, {2, 4, "AATT"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 2}});
      _path_vector.emplace_back(RRPath{"1", std::list<size_t>{1, 3}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}},
                                  {1, {"GGA", false}},
                                  {4, {"ATT", false}},
                                  {3, {"TGC", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 3, "ACAAATGC"}, {1, 4, "GGAAATT"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with a complex vertex (loop)
TEST(DBComplexVertexLoop1, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 2, "ACAAA"}, {2, 2, "AAGAA"}, {2, 3, "AATGC"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 1}});
      _path_vector.emplace_back(RRPath{"1", std::list<size_t>{1, 2}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}}, {3, {"TGC", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 3, "ACAAAGAATGC"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with a complex vertex (loop + another traversal)
TEST(DBComplexVertexLoop2, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 2, "ACAAA"},
        {2, 2, "AAGAA"},
        {2, 3, "AATGC"},
        {4, 2, "GGAA"},
        {2, 5, "AATG"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 1}});
      _path_vector.emplace_back(RRPath{"1", std::list<size_t>{1, 2}});
      _path_vector.emplace_back(RRPath{"2", std::list<size_t>{3, 4}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}},
                                  {3, {"TGC", false}},
                                  {4, {"GGA", false}},
                                  {5, {"ATG", false}}};

        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 3, "ACAAAGAATGC"}, {4, 5, "GGAATG"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with a complex vertex (two loops)
TEST(DBComplexVertexLoop3, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 2, "ACAAA"}, {2, 2, "AAGAA"}, {2, 3, "AATGC"},
        {4, 2, "GGAA"}, {2, 2, "AAA"}, {2, 5, "AATG"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 1, 2}});
      _path_vector.emplace_back(RRPath{"1", std::list<size_t>{3, 4, 5}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}},
                                  {3, {"TGC", false}},
                                  {4, {"GGA", false}},
                                  {5, {"ATG", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 3, "ACAAAGAATGC"}, {4, 5, "GGAAATG"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with a complex vertex (multiple loops)
TEST(DBComplexVertexLoop4, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "ACAAA"}, {1, 1, "AAGAA"}, {1, 1, "AACAA"},
        {1, 1, "AATAA"}, {1, 1, "AAAAA"}, {1, 2, "AATGC"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector
          .emplace_back(RRPath{"0", std::list<size_t>{0, 1, 2, 3, 4, 5}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}}, {2, {"TGC", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 2, "ACAAAGAACAATAAAAATGC"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with a complex vertex (multiple loops, several traversals)
TEST(DBComplexVertexLoop5, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "ACAAA"},  // 0
        {1, 1, "AAGAA"},  // 1
        {1, 1, "AACAA"},  // 2
        {1, 1, "AATAA"},  // 3
        {1, 1, "AAAAA"},  // 4
        {1, 2, "AATGC"},  // 5
        {3, 1, "ACAAA"},  // 6
        {1, 4, "AATGC"},  // 7
        {5, 1, "ACAAA"},  // 8
        {1, 6, "AATGC"}}; // 9
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 1, 2, 5}});
      _path_vector.emplace_back(RRPath{"1", std::list<size_t>{6, 3, 4, 7}});
      _path_vector.emplace_back(RRPath{"2", std::list<size_t>{8, 9}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {

        RawVertexInfo vertex_info{{0, {"ACA", false}}, {2, {"TGC", false}},
                                  {3, {"ACA", false}}, {4, {"TGC", false}},
                                  {5, {"ACA", false}}, {6, {"TGC", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 2, "ACAAAGAACAATGC"}, {3, 4, "ACAAATAAAAATGC"},
            {5, 6, "ACAAATGC"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

/*
// graph with two buldges and loops inside
TEST(DBBuldges1, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "ACAAA"},  // 0
        {1, 1, "AAGAA"},  // 1
        {1, 2, "AACGC"},  // 2
        {0, 1, "ACTAA"},  // 3
        {1, 1, "AAAAA"},  // 4
        {1, 2, "AATGC"},  // 5
        {0, 1, "ACAAA"},  // 6
        {1, 2, "AATGC"}}; // 7
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 1, 2}});
      _path_vector.emplace_back(RRPath{"1", std::list<size_t>{3, 4, 5}});
      _path_vector.emplace_back(RRPath{"2", std::list<size_t>{6, 7}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{16, {"ACA", false}}, {3, {"CGC", false}},
                                  {17, {"ACT", false}}, {4, {"TGC", false}},
                                  {18, {"ACA", false}}, {5, {"TGC", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {16, 3, "ACAAAGAACGC"}, {17, 4, "ACTAAAAATGC"},
            {18, 5, "ACAAATGC"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}
*/

// graph with complex vertex and 4 connections
TEST(DBComplexVertexConn4, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 2, "ACAAA"}, {1, 2, "GGAAA"}, {2, 3, "AATGC"}, {2, 4, "AATT"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 2}});
      _path_vector.emplace_back(RRPath{"1", std::list<size_t>{0, 3}});
      _path_vector.emplace_back(RRPath{"2", std::list<size_t>{1, 2}});
      _path_vector.emplace_back(RRPath{"3", std::list<size_t>{1, 3}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}}, {1, {"GGA", false}},
                                  {3, {"TGC", false}}, {4, {"ATT", false}},
                                  {5, {"AAA", false}}, {6, {"AAA", false}},
                                  {7, {"AAT", false}}, {8, {"AAT", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 5, "ACAAA"}, {5, 7, "AAAT"}, {7, 3, "AATGC"}, {5, 8, "AAAT"},
            {1, 6, "GGAAA"}, {6, 7, "AAAT"}, {6, 8, "AAAT"}, {8, 4, "AATT"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with complex vertex and 3 connections
TEST(DBComplexVertexConn3, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 2, "ACAAA"}, {1, 2, "GGAAA"}, {2, 3, "AATGC"}, {2, 4, "AATT"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 2}});
      _path_vector.emplace_back(RRPath{"1", std::list<size_t>{0, 3}});
      // _path_vector.emplace_back(RRPath{"2", std::list<size_t>{1, 2}});
      _path_vector.emplace_back(RRPath{"3", std::list<size_t>{1, 3}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}}, {1, {"GGA", false}},
                                  {3, {"TGC", false}}, {4, {"ATT", false}},
                                  {5, {"AAA", false}}, {8, {"AAT", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 5, "ACAAA"},
            {5, 3, "AAATGC"},
            {5, 8, "AAAT"},
            {1, 8, "GGAAAT"},
            {8, 4, "AATT"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with complex vertex and 3 connections
TEST(DBComplexVertexConn3_2, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 2, "ACAAA"}, {1, 2, "GGAAA"}, {2, 3, "AATGC"}, {2, 4, "AATT"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 2}});
      // _path_vector.emplace_back(RRPath{"1", std::list<size_t>{0, 3}});
      _path_vector.emplace_back(RRPath{"2", std::list<size_t>{1, 2}});
      _path_vector.emplace_back(RRPath{"3", std::list<size_t>{1, 3}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}}, {1, {"GGA", false}},
                                  {3, {"TGC", false}}, {4, {"ATT", false}},
                                  {6, {"AAA", false}}, {7, {"AAT", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 7, "ACAAAT"},
            {7, 3, "AATGC"},
            {6, 7, "AAAT"},
            {1, 6, "GGAAA"},
            {6, 4, "AAATT"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with a complex vertex (loop)
TEST(DBComplexVertexLoop6, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "ACAAA"}, {1, 1, "AAGAA"}, {1, 2, "AATGC"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 1, 1, 2}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    int N = 4;
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + N, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo
            vertex_info{{0, {"ACAAAG", false}}, {2, {"GAATGC", false}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 2, "ACAAAGAAGAATGC"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with a complex vertex (loop)
TEST(DBComplexVertexLoop7, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        // {0, 1, "ACAAA"},
        {1, 1, "AAGAA"}};
    // {1, 2, "AATGC"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      // _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 2}});
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 0}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{1, {"AA", true}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {1, 1, "AAGAA"}};

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

// graph with a complex vertex (loop)
TEST(DBComplexVertexLoop8, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "ACAAA"}, {1, 1, "AAGAA"}, {1, 2, "AATGC"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 2}});
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{1, 1}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}},
                                  {2, {"TGC", false}},
                                  {6, {"AAG", true}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 2, "ACAAATGC"}, {6, 6, "AAGAAG"}};

//        {
//            for (auto vertex : mdbg) {
//                auto [begin, end] = mdbg.out_neighbors(vertex);
//                for (auto it = begin; it != end; ++it) {
//                    std::cout << vertex << " " << it->first << " " << it->second.prop().Size() << " " <<
//                        mdbg.GetEdgeSequence(mdbg.find(vertex), it, false, false).ToSequence().str() << "\n";
//                }
//            }
//        }
        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

/*
// graph with a complex vertex (loop) + rc
TEST(DBComplexVertexLoop9RC, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "ACAAA"}, {1, 1, "AACGTTGCAA"}, {1, 2, "AATGC"},
        {3, 4, "GCATT"}, {4, 4, "TTGCAACGTT"}, {4, 5, "TTTGT"},};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 2}});
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{1, 1}});
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{3, 5}});
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{4, 4}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, true);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", false}},
                                  {2, {"TGC", false}},
                                  {9, {"AAC", true}},
                                  {3, {"GCA", false}},
                                  {5, {"TGT", false}},
                                  {11, {"GTT", true}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {0, 2, "ACAAATGC"}, {9, 9,   "AACGTTGCAAC"},
            {3, 5, "GCATTTGT"}, {11, 11, "GTTGCAACGTT"}};

//        {
//            for (auto vertex : mdbg) {
//                auto [begin, end] = mdbg.out_neighbors(vertex);
//                for (auto it = begin; it != end; ++it) {
//                    std::cout << vertex << " " << it->first << " " <<
//                        mdbg.GetEdgeSequence(mdbg.find(vertex), it, false, false).ToSequence().str() << "\n";
//                }
//            }
//        }
        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
        mdbg.AssertValidity();
    }
}
 */

// graph with a two loops + rc
TEST(DBDoubleLoopRC, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 0, "AACGTCGCAA"}, {1, 1, "TTGCGACGTT"},
        {0, 0, "AAA"}, {1, 1, "TTT"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{0, 2, 0}});
      _path_vector.emplace_back(RRPath{"0", std::list<size_t>{1, 3, 1}});

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, true);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{9, {"AAA", true}},
                                  {5, {"TTT", true}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge{
            {9, 9, "AAACGTCGCAAA"}, {5, 5, "TTTGCGACGTTT"}};

//        {
//            for (auto vertex : mdbg) {
//                auto [begin, end] = mdbg.out_neighbors(vertex);
//                for (auto it = begin; it != end; ++it) {
//                    std::cout << vertex << " " << it->first << " " <<
//                              mdbg.GetEdgeSequence(mdbg.find(vertex), it, false, false).ToSequence().str() << "\n";
//                }
//            }
//        }
        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
        mdbg.AssertValidity();
    }
}

// graph with a single edge that will be isolated
TEST(DBIsolate, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info{
        {0, 1, "ACA"}};
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info{{0, {"ACA", true}}};
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge;

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
    }
}

TEST(DBEmptyGraph, Basic) {
    const size_t k = 2;

    std::vector<std::tuple<uint64_t, uint64_t, std::string>> raw_edge_info;
    std::map<RRVertexType, dbg::Vertex> vertexes;
    std::vector<dbg::Edge> edges;
    std::vector<SuccinctEdgeInfo> edge_info =
        GetEdgeInfo(vertexes, edges, raw_edge_info, k, false);

    RRPaths paths = []() {
      std::vector<RRPath> _path_vector;

      return PathsBuilder::FromPathVector(_path_vector);
    }();

    MultiplexDBG mdbg(edge_info, k, &paths, false);
    logging::Logger logger;

    MultiplexDBGIncreaser k_increaser{k, k + 1, logger, true};
    k_increaser.IncreaseUntilSaturation(mdbg);
    {
        RawVertexInfo vertex_info;
        std::vector<std::tuple<uint64_t, uint64_t, std::string>> post_raw_edge;

        auto[VertexIndexSetsEqual, VertexPropsEquals] =
        CompareVertexes(mdbg, vertex_info);
        ASSERT_TRUE(VertexIndexSetsEqual);
        ASSERT_TRUE(VertexPropsEquals);
        ASSERT_TRUE(CompareEdges(mdbg, post_raw_edge));
        ASSERT_TRUE(mdbg.IsFrozen());
    }
}