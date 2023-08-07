#pragma once

#include <dbg/paths.hpp>
#include <ksw2/ksw_wrapper.hpp>

class AlignmentForm {
private:
    std::vector<CigarPair> cigar;
    size_t qlen;
    size_t tlen;
    struct AlignmentColumn {
        size_t qpos;
        size_t tpos;
        CigarEvent event;
        AlignmentColumn(size_t qpos, size_t tpos, CigarEvent event) : qpos(qpos), tpos(tpos), event(event) {}
    };
    class AlignmentFormIterator {
        AlignmentForm *alignmentForm;
        size_t cigar_pos;
        size_t block_pos;
        size_t qpos;
        size_t tpos;
    public:
        AlignmentFormIterator(AlignmentForm &form, size_t cigar_pos, size_t block_pos, size_t qpos, size_t tpos) :
                alignmentForm(&form), cigar_pos(cigar_pos), block_pos(block_pos), qpos(qpos), tpos(tpos) {
        }
        AlignmentFormIterator(AlignmentForm &form, size_t cigar_pos, size_t block_pos) :
                alignmentForm(&form), cigar_pos(cigar_pos), block_pos(block_pos), qpos(0), tpos(0) {
            for(size_t i = 0; i < cigar_pos; i++) {
                CigarPair p = form.cigar[i];
                if(p.type != CigarEvent::D) {
                    qpos += p.length;
                }
                if(p.type != CigarEvent::I) {
                    tpos += p.length;
                }
            }
            if(block_pos > 0) {
                VERIFY(cigar_pos < form.cigar.size());
                CigarPair p = form.cigar[cigar_pos];
                if(p.type != CigarEvent::D) {
                    qpos += p.length;
                }
                if(p.type != CigarEvent::I) {
                    tpos += p.length;
                }
            }
        }
        AlignmentFormIterator &operator++() {
            block_pos++;
            if(alignmentForm->cigar[cigar_pos].length == block_pos) {
                cigar_pos++;
                block_pos = 0;
            }
            return *this;
        }

        AlignmentFormIterator operator++(int) const {
            AlignmentFormIterator res = *this;
            ++res;
            return res;
        }

        AlignmentColumn operator*() const {
            return {qpos, tpos, alignmentForm->cigar[cigar_pos].type};
        }
    };

    void calculateLens() {
        qlen = 0;
        tlen = 0;
        for(CigarPair p : cigar) {
            if(p.type != CigarEvent::D) {
                qlen += p.length;
            }
            if(p.type != CigarEvent::I) {
                tlen += p.length;
            }
        }
    }
public:
    AlignmentForm(std::vector<CigarPair> _cigar) : cigar(std::move(_cigar)), qlen(0), tlen(0) {
        calculateLens();
    }
    AlignmentForm(const std::string &s) {
        size_t n = 0;
        size_t pos = s.find_last_of(':') + 1;
        for(size_t i = pos; i < s.size(); i++){
            char c = s[i];
            if(c >= '0' && c <= '9') {
                n = n * 10 + c - '0';
            } else {
                if (n == 0)
                    n = 1;
                if (c == '=' || c == 'X') {
                    c = 'M';
                }
                if(!cigar.empty() && c == cigar.back().type) {
                    cigar.back() = {c, cigar.back().length + n};
                } else {
                    cigar.emplace_back(c, n);
                }
                n = 0;
            }
        }
        calculateLens();
    }


    bool empty() const {return cigar.empty();}
    size_t queryLength() const {return qlen;}
    size_t targetLength() const {return tlen;}

    std::vector<CigarPair>::const_iterator begin() const {return cigar.begin();}
    std::vector<CigarPair>::const_iterator end() const {return cigar.end();}

    void operator+=(const AlignmentForm &other) {
        if(empty()) {
            *this = other;
            return;
        }
        for(CigarPair p: other) {
            if(p.type == cigar.back().type) {
                cigar.back().length += p.length;
            } else {
                cigar.emplace_back(p);
            }
        }
        qlen += other.qlen;
        tlen += other.tlen;
    }
    AlignmentForm operator+(const AlignmentForm &other) const {
        AlignmentForm res = *this;
        res += other;
        return std::move(res);
    }

    AlignmentForm RC() const {
        return {std::vector<CigarPair>(cigar.rbegin(), cigar.rend())};
    }
};

template<class Graph>
class GraphAlignment {
private:
    Segment<Contig> segment;
    GraphPath<Graph> path;
    AlignmentForm alignment;
public:
    GraphAlignment(Segment<Contig> segment, GraphPath<Graph> path, AlignmentForm alignment) :
                segment(segment), path(std::move(path)), alignment(std::move(alignment)) {
    }

    void operator+=(const GraphAlignment<Graph> &other) {
        alignment += other.alignment;
        path += other.path;
        segment = segment + other.segment;
    }

    GraphAlignment operator+(const GraphAlignment<Graph> &other) const {
        GraphAlignment<Graph> res = *this;
        res += other;
        return std::move(res);
    }
};