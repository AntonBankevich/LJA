// NuclDeck against a plain string model: random single and bulk pushes (deck, Sequence) and pops at both ends, so the
// circular buffer wraps around and grows in every state, with bulk words crossing the buffer end.

#include "gtest/gtest.h"
#include "sequences/nucl_buffer.hpp"
#include "sequences/sequence.hpp"
#include <random>
#include <string>

namespace {
    std::string Str(const NuclDeck &deck) {
        std::string res;
        for (size_t i = 0; i < deck.size(); i++)
            res += "ACGT"[deck[i]];
        return res;
    }

    std::string RandomNucls(std::mt19937 &rng, size_t len) {
        std::string res;
        for (size_t i = 0; i < len; i++)
            res += "ACGT"[rng() % 4];
        return res;
    }

    NuclDeck DeckOf(const std::string &s) {
        NuclDeck deck;
        for (char c : s)
            deck.push_back((unsigned char) dignucl(c));
        return deck;
    }
}

TEST(NuclDeck, RandomOperationsMatchStringModel) {
    std::mt19937 rng(17);
    for (size_t round = 0; round < 50; round++) {
        NuclDeck deck(1 + rng() % 40);
        std::string model;
        for (size_t step = 0; step < 300; step++) {
            size_t op = rng() % 9;
            size_t len = rng() % 100;
            std::string add = RandomNucls(rng, len);
            if (op == 0) {
                deck.push_back((unsigned char) dignucl(add.empty() ? 'A' : add[0]));
                model += add.empty() ? "A" : add.substr(0, 1);
            } else if (op == 1) {
                deck.push_front((unsigned char) dignucl(add.empty() ? 'C' : add[0]));
                model = (add.empty() ? "C" : add.substr(0, 1)) + model;
            } else if (op == 2) {
                deck.push_back(DeckOf(add));
                model += add;
            } else if (op == 3) {
                deck.push_front(DeckOf(add));
                model = add + model;
            } else if (op == 4) {
                deck.push_back(Sequence(add));
                model += add;
            } else if (op == 5) {
                deck.push_front(Sequence(add));
                model = add + model;
            } else if (op == 6 && !model.empty()) {
                size_t cnt = rng() % (model.size() + 1);
                deck.pop_back(cnt);
                model.resize(model.size() - cnt);
            } else if (op == 7 && !model.empty()) {
                size_t cnt = rng() % (model.size() + 1);
                deck.pop_front(cnt);
                model = model.substr(cnt);
            } else if (op == 8) {
                deck.push_back(deck);
                model += model;
                if (model.size() > 5000) {
                    deck.pop_front(model.size() - 100);
                    model = model.substr(model.size() - 100);
                }
            }
            ASSERT_EQ(deck.size(), model.size()) << "round " << round << " step " << step << " op " << op;
            ASSERT_EQ(Str(deck), model) << "round " << round << " step " << step << " op " << op;
            ASSERT_TRUE(deck == DeckOf(model)) << "round " << round << " step " << step;
        }
    }
}

// Iterators store logical positions: an iterator at a nucleotide that survives an operation (no pop removed it) still
// points to that nucleotide afterwards, whatever pushes, pops, wrap-arounds and reallocations happened.
TEST(NuclDeck, IteratorsSurviveChanges) {
    struct Probe {
        NuclDeck::Iterator it;
        size_t idx;
        char c;
    };
    std::mt19937 rng(23);
    for (size_t round = 0; round < 50; round++) {
        NuclDeck deck(1 + rng() % 40);
        std::string model;
        for (size_t step = 0; step < 300; step++) {
            std::vector<Probe> probes;
            for (size_t i = 0; i < 8 && !model.empty(); i++) {
                size_t idx = rng() % model.size();
                probes.push_back({deck.begin() + int(idx), idx, model[idx]});
            }
            size_t op = rng() % 8;
            std::string add = RandomNucls(rng, 1 + rng() % 100);
//            front: nucleotides added (> 0) or removed (< 0) in front; keep: old indices below it survive at the back.
            long front = 0;
            size_t keep = model.size();
            if (op == 0) {
                deck.push_back((unsigned char) dignucl(add[0]));
                model += add[0];
            } else if (op == 1) {
                deck.push_front((unsigned char) dignucl(add[0]));
                model = add[0] + model;
                front = 1;
            } else if (op == 2) {
                deck.push_back(DeckOf(add));
                model += add;
            } else if (op == 3) {
                deck.push_front(DeckOf(add));
                model = add + model;
                front = long(add.size());
            } else if (op == 4) {
                deck.push_back(Sequence(add));
                model += add;
            } else if (op == 5) {
                deck.push_front(Sequence(add));
                model = add + model;
                front = long(add.size());
            } else if (op == 6 && !model.empty()) {
                size_t cnt = rng() % (model.size() + 1);
                deck.pop_back(cnt);
                model.resize(model.size() - cnt);
                keep = model.size();
            } else if (op == 7 && !model.empty()) {
                size_t cnt = rng() % (model.size() + 1);
                deck.pop_front(cnt);
                model = model.substr(cnt);
                front = -long(cnt);
            }
            for (const Probe &p : probes) {
                long new_idx = long(p.idx) + front;
                if (p.idx < keep && new_idx >= 0) {
                    ASSERT_EQ("ACGT"[*p.it], p.c) << "round " << round << " step " << step << " op " << op;
                    ASSERT_EQ(model[new_idx], p.c) << "round " << round << " step " << step << " op " << op;
                    ASSERT_EQ(p.it.getPos(), size_t(new_idx)) << "round " << round << " step " << step << " op " << op;
                }
            }
        }
    }
}

TEST(NuclDeck, EqualityDetectsSingleDifference) {
    std::mt19937 rng(3);
    std::string s = RandomNucls(rng, 200);
    for (size_t pos : {0, 31, 32, 33, 100, 199}) {
        std::string t = s;
        t[pos] = t[pos] == 'A' ? 'C' : 'A';
        EXPECT_FALSE(DeckOf(s) == DeckOf(t)) << "difference at " << pos;
    }
    EXPECT_TRUE(DeckOf(s) == DeckOf(s));
}

// An empty deck holds no buffer until its first push; a moved-from deck is empty and usable again.
TEST(NuclDeck, EmptyAndMovedFromDecks) {
    NuclDeck empty;
    EXPECT_TRUE(empty.empty());
    empty.push_front(DeckOf("ACG"));
    EXPECT_EQ(Str(empty), "ACG");
    NuclDeck a = DeckOf("ACGTACGTTT");
    NuclDeck b(std::move(a));
    EXPECT_EQ(Str(b), "ACGTACGTTT");
    EXPECT_TRUE(a.empty());
    a.push_back(DeckOf("GG"));
    EXPECT_EQ(Str(a), "GG");
    NuclDeck c;
    c = std::move(b);
    EXPECT_EQ(Str(c), "ACGTACGTTT");
    EXPECT_TRUE(b.empty());
    b.push_back((unsigned char) 2);
    EXPECT_EQ(Str(b), "G");
}
