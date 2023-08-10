#pragma once

#include "common/iterator_utils.hpp"
#include <ksw2/ksw_wrapper.hpp>
#include <string>

class AlignmentForm {
private:
    std::vector<CigarPair> cigar;
    size_t qlen;
    size_t tlen;

    void calculateLens();

public:
    struct AlignmentColumn {
        size_t qpos;
        size_t tpos;
        CigarEvent event;
        AlignmentColumn(size_t qpos, size_t tpos, CigarEvent event) : qpos(qpos), tpos(tpos), event(event) {}
    };

    class AlignmentColumnIterator {
        AlignmentForm *alignmentForm;
        size_t cigar_pos;
        size_t block_pos;
        size_t qpos;
        size_t tpos;
    public:
        AlignmentColumnIterator(AlignmentForm &form, size_t cigar_pos, size_t block_pos, size_t qpos, size_t tpos) :
                alignmentForm(&form), cigar_pos(cigar_pos), block_pos(block_pos), qpos(qpos), tpos(tpos) {
        }
        AlignmentColumnIterator(AlignmentForm &form, size_t cigar_pos, size_t block_pos);

        AlignmentColumnIterator &operator++();
        AlignmentColumnIterator operator++(int) const;
        AlignmentColumnIterator &operator--();
        AlignmentColumnIterator operator--(int) const;

        AlignmentColumn operator*() const {return {qpos, tpos, alignmentForm->cigar[cigar_pos].type};}

        bool operator==(const AlignmentColumnIterator &other) const {return cigar_pos == other.cigar_pos && block_pos == other.block_pos;}
        bool operator!=(const AlignmentColumnIterator &other) const {return !(*this == other);}
    };

    AlignmentForm() : qlen(0), tlen(0) {}
    AlignmentForm(std::vector<CigarPair> _cigar) : cigar(std::move(_cigar)), qlen(0), tlen(0) {calculateLens();}
    AlignmentForm(const std::string &s);


    bool empty() const {return cigar.empty();}
    size_t queryLength() const {return qlen;}
    size_t targetLength() const {return tlen;}

    std::vector<CigarPair>::const_iterator begin() const {return cigar.begin();}
    std::vector<CigarPair>::const_iterator end() const {return cigar.end();}

    CigarPair &front() {return cigar.front();}
    const CigarPair &front() const {return cigar.front();}
    CigarPair &back() {return cigar.back();}
    const CigarPair &back() const {return cigar.back();}
    CigarPair &operator[](size_t ind) {return cigar[ind];}
    const CigarPair &operator[](size_t ind) const {return cigar[ind];}

    IterableStorage<AlignmentColumnIterator> columns() {
        return {{*this, 0, 0, 0, 0}, {*this, cigar.size(), 0, qlen, tlen}};
    }

    void operator+=(const AlignmentForm &other);
    void operator+=(CigarEvent e);
    void operator+=(CigarPair p);
    void pop_back() {cigar.pop_back();}
    void pop_front() {cigar = {begin() + 1, end()};}

    AlignmentForm operator+(const AlignmentForm &other) const;
    AlignmentForm RC() const;
    AlignmentForm Reverse() const;
    AlignmentForm Prefix(AlignmentColumnIterator bound);
    AlignmentForm Suffix(AlignmentColumnIterator bound);

    template<class U, class V>
    std::vector<std::string> toString(U &from_seq, V &to_seq,
                                 const std::vector<std::pair<size_t, size_t>> &to_ignore = {}) const {
        size_t from_pos = 0;
        size_t to_pos = 0;
        size_t cur = 0;
        std::vector<char> s1, s2, m;
        for(CigarPair cp: *this) {
            for(size_t t = 0; t < cp.length; t++) {
                if(cur < to_ignore.size() && to_ignore[cur].second == from_pos)
                    cur++;
                if(cur < to_ignore.size() && from_pos < to_ignore[cur].second && from_pos >= to_ignore[cur].first) {
                    m.emplace_back('*');
                } else if(cp.type == 'M' && from_seq[from_pos] == to_seq[to_pos])
                    m.emplace_back('|');
                else
                    m.emplace_back('.');
                if(cp.type == 'D') {
                    s1.emplace_back('-');
                } else {
                    s1.emplace_back(nucl(from_seq[from_pos]));
                    from_pos++;
                }
                if(cp.type == 'I') {
                    s2.emplace_back('-');
                } else {
                    s2.emplace_back(nucl(to_seq[to_pos]));
                    to_pos++;
                }
            }
        }
        std::vector<std::string> res = {std::string(s1.begin(), s1.end()), std::string(m.begin(), m.end()), std::string(s2.begin(), s2.end())};
        return res;
    }
};

class AlignmentHelper {
public:
    template<class U, class V>
    static AlignmentForm::AlignmentColumnIterator
    LastLongMatch(U tseq, V qseq, AlignmentForm &extension, size_t match) {
        size_t cnt = 0;
        auto iter = extension.columns().end();
        while(iter != extension.columns().begin()) {
            --iter;
            auto column = *iter;
            if(column.event == M && tseq[column.tpos] == qseq[column.qpos]) {
                cnt++;
            } else {
                cnt = 0;
            }
            if(cnt == match) {
                for(size_t i = 0; i < match; i++)
                    ++iter;
                break;
            }
        }
        return iter;
    }
};
