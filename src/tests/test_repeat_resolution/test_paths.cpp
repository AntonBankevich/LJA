//
// Created by Andrey Bzikadze on 12/10/21.
//

#include "gtest/gtest.h"
#include "repeat_resolution/paths.hpp"

using namespace repeat_resolution;

TEST(RRPathsTest, Basic) {
    std::vector<RRPath> _path_vector;
    _path_vector.emplace_back(
        RRPath{"0", std::list<size_t>{1, 2, 3, 4, 5, 2, 6, 7, 8, 9, 10}});
    _path_vector.emplace_back(
        RRPath{"1", std::list<size_t>{11, 12, 2, 13, 14, 15, 2, 17, 18}});
    _path_vector.emplace_back(RRPath{"2", std::list<size_t>{2}});
    _path_vector.emplace_back(RRPath{"3", std::list<size_t>{2, 19}});
    _path_vector.emplace_back(RRPath{"4", std::list<size_t>{5, 2}});

    RRPaths paths = PathsBuilder::FromPathVector(_path_vector);
    const auto &path_vector = paths.GetPaths();
    const auto &ei2p = paths.GetEdge2Pos();
    const auto &eip2p = paths.GetEdgepair2Pos();
    {
        std::vector<RRPath> path_vector_ref;
        path_vector_ref.emplace_back(
            RRPath{"0", std::list<size_t>{1, 2, 3, 4, 5, 2, 6, 7, 8, 9, 10}});
        path_vector_ref.emplace_back(
            RRPath{"1", std::list<size_t>{11, 12, 2, 13, 14, 15, 2, 17, 18}});
        path_vector_ref.emplace_back(RRPath{"2", std::list<size_t>{2}});
        path_vector_ref.emplace_back(RRPath{"3", std::list<size_t>{2, 19}});
        path_vector_ref.emplace_back(RRPath{"4", std::list<size_t>{5, 2}});
        ASSERT_EQ(_path_vector, path_vector_ref);
    }

    {
        std::unordered_map<size_t, size_t> index_cnt_ref{
            {1, 1}, {2, 7}, {3, 1}, {4, 1}, {5, 2}, {6, 1},
            {7, 1}, {8, 1}, {9, 1}, {10, 1}, {11, 1}, {12, 1},
            {13, 1}, {14, 1}, {15, 1}, {17, 1}, {18, 1}, {19, 1}};
        for (const auto &pair : ei2p) {
            ASSERT_NE(index_cnt_ref.find(pair.first), index_cnt_ref.end());
            ASSERT_EQ(pair.second.size(), index_cnt_ref.at(pair.first));
        }
    }
    {
        std::unordered_map<PairEdgeIndexType, size_t, PairEdgeIndexHash>
            index_cnt_ref{{{1, 2}, 1}, {{2, 3}, 1}, {{2, 19}, 1}, {{20, 2}, 1},
                          {{3, 4}, 1}, {{4, 5}, 1}, {{5, 2}, 2}, {{2, 6}, 1},
                          {{6, 7}, 1}, {{7, 8}, 1}, {{8, 9}, 1}, {{9, 10}, 1},
                          {{11, 12}, 1}, {{12, 2}, 1}, {{2, 13}, 1},
                          {{13, 14}, 1},
                          {{14, 15}, 1}, {{15, 2}, 1}, {{2, 17}, 1},
                          {{17, 18}, 1}};
        for (const auto &pair : eip2p) {
            ASSERT_NE(index_cnt_ref.find(pair.first), index_cnt_ref.end());
            ASSERT_EQ(pair.second.size(), index_cnt_ref.at(pair.first));
        }
    }

    paths.Remove(2);
    paths.assert_validity();
    {
        std::vector<RRPath> path_vector_ref;
        path_vector_ref.emplace_back(
            RRPath{"0", std::list<size_t>{1, 3, 4, 5, 6, 7, 8, 9, 10}});
        path_vector_ref.emplace_back(
            RRPath{"1", std::list<size_t>{11, 12, 13, 14, 15, 17, 18}});
        path_vector_ref.emplace_back(RRPath{"2", std::list<size_t>{}});
        path_vector_ref.emplace_back(RRPath{"3", std::list<size_t>{19}});
        path_vector_ref.emplace_back(RRPath{"4", std::list<size_t>{5}});
        ASSERT_EQ(path_vector, path_vector_ref);
    }
    {
        std::unordered_map<size_t, size_t> index_cnt_ref{
            {1, 1}, {3, 1}, {4, 1}, {5, 2}, {6, 1}, {7, 1},
            {8, 1}, {9, 1}, {10, 1}, {11, 1}, {12, 1}, {13, 1},
            {14, 1}, {15, 1}, {17, 1}, {18, 1}, {19, 1}};
        for (const auto &pair : ei2p) {
            ASSERT_NE(index_cnt_ref.find(pair.first), index_cnt_ref.end());
            ASSERT_EQ(pair.second.size(), index_cnt_ref.at(pair.first));
        }
    }
    {
        std::unordered_map<PairEdgeIndexType, size_t, PairEdgeIndexHash>
            index_cnt_ref{{{1, 3}, 1}, {{3, 4}, 1}, {{4, 5}, 1},
                          {{5, 6}, 1}, {{6, 7}, 1}, {{7, 8}, 1},
                          {{8, 9}, 1}, {{9, 10}, 1}, {{11, 12}, 1},
                          {{12, 13}, 1}, {{13, 14}, 1}, {{14, 15}, 1},
                          {{15, 17}, 1}, {{17, 18}, 1}};
        for (const auto &pair : eip2p) {
            ASSERT_NE(index_cnt_ref.find(pair.first), index_cnt_ref.end());
            ASSERT_EQ(pair.second.size(), index_cnt_ref.at(pair.first));
        }
    }
    paths.Add(1, 3, 2);
    paths.assert_validity();
    {
        std::vector<RRPath> path_vector_ref;
        path_vector_ref.emplace_back(
            RRPath{"0", std::list<size_t>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}});
        path_vector_ref.emplace_back(
            RRPath{"1", std::list<size_t>{11, 12, 13, 14, 15, 17, 18}});
        path_vector_ref.emplace_back(RRPath{"2", std::list<size_t>{}});
        path_vector_ref.emplace_back(RRPath{"3", std::list<size_t>{19}});
        path_vector_ref.emplace_back(RRPath{"4", std::list<size_t>{5}});
        ASSERT_EQ(path_vector, path_vector_ref);
    }
    {
        std::unordered_map<size_t, size_t> index_cnt_ref{
            {1, 1}, {2, 1}, {3, 1}, {4, 1}, {5, 2}, {6, 1}, {7, 1},
            {8, 1}, {9, 1}, {10, 1}, {11, 1}, {12, 1}, {13, 1}, {14, 1},
            {15, 1}, {17, 1}, {18, 1}, {19, 1}, {20, 1}};
        for (const auto &pair : ei2p) {
            ASSERT_NE(index_cnt_ref.find(pair.first), index_cnt_ref.end());
            ASSERT_EQ(pair.second.size(), index_cnt_ref.at(pair.first));
        }
    }
    {
        std::unordered_map<PairEdgeIndexType, size_t, PairEdgeIndexHash>
            index_cnt_ref{{{1, 2}, 1}, {{2, 3}, 1}, {{3, 4}, 1},
                          {{4, 5}, 1}, {{5, 6}, 1}, {{6, 7}, 1},
                          {{7, 8}, 1}, {{8, 9}, 1}, {{9, 10}, 1},
                          {{11, 12}, 1}, {{12, 13}, 1}, {{13, 14}, 1},
                          {{14, 15}, 1}, {{15, 17}, 1}, {{17, 18}, 1}};
        for (const auto &pair : eip2p) {
            ASSERT_NE(index_cnt_ref.find(pair.first), index_cnt_ref.end());
            ASSERT_EQ(pair.second.size(), index_cnt_ref.at(pair.first));
        }
    }

    paths.Merge(4, 5);
    paths.assert_validity();
    {
        std::vector<RRPath> path_vector_ref;
        path_vector_ref.emplace_back(
            RRPath{"0", std::list<size_t>{1, 2, 3, 4, 6, 7, 8, 9, 10}});
        path_vector_ref.emplace_back(
            RRPath{"1", std::list<size_t>{11, 12, 13, 14, 15, 17, 18}});
        path_vector_ref.emplace_back(RRPath{"2", std::list<size_t>{}});
        path_vector_ref.emplace_back(RRPath{"3", std::list<size_t>{19}});
        path_vector_ref.emplace_back(RRPath{"4", std::list<size_t>{4}});
        ASSERT_EQ(path_vector, path_vector_ref);
    }
    {
        std::unordered_map<size_t, size_t> index_cnt_ref{
            {1, 1}, {2, 1}, {3, 1}, {4, 2}, {6, 1}, {7, 1},
            {8, 1}, {9, 1}, {10, 1}, {11, 1}, {12, 1}, {13, 1},
            {14, 1}, {15, 1}, {17, 1}, {18, 1}, {19, 1}, {20, 1}};
        for (const auto &pair : ei2p) {
            ASSERT_NE(index_cnt_ref.find(pair.first), index_cnt_ref.end());
            ASSERT_EQ(pair.second.size(), index_cnt_ref.at(pair.first));
        }
    }
    {
        std::unordered_map<PairEdgeIndexType, size_t, PairEdgeIndexHash>
            index_cnt_ref{{{1, 2}, 1}, {{2, 3}, 1}, {{3, 4}, 1},
                          {{4, 6}, 1}, {{6, 7}, 1}, {{7, 8}, 1},
                          {{8, 9}, 1}, {{9, 10}, 1}, {{11, 12}, 1},
                          {{12, 13}, 1}, {{13, 14}, 1}, {{14, 15}, 1},
                          {{15, 17}, 1}, {{17, 18}, 1}};
        for (const auto &pair : eip2p) {
            ASSERT_NE(index_cnt_ref.find(pair.first), index_cnt_ref.end());
            ASSERT_EQ(pair.second.size(), index_cnt_ref.at(pair.first));
        }
    }
}

TEST(RRPathsTest, MergeIterDereference) {
    std::vector<RRPath> _path_vector;
    _path_vector.emplace_back(RRPath{"0", std::list<size_t>{1, 2}});
    _path_vector.emplace_back(RRPath{"1", std::list<size_t>{2, 3}});

    RRPaths paths = PathsBuilder::FromPathVector(_path_vector);
    paths.Merge(1, 2);
}
