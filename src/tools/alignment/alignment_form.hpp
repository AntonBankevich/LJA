#pragma once

#include "common/iterator_utils.hpp"
#include "sequences/nucl.hpp"
#include <vector>
#include <string>
#include <sstream>

enum CigarEvent : char {
    I = 'I', D = 'D', M = 'M',
};

inline CigarEvent CigarEventFromChar(char c) {
    if(c == 'I')
        return I;
    else if (c == 'D')
        return D;
    else
        return M;
}

struct CigarPair {
    CigarEvent type;
    size_t length;
    CigarPair(CigarEvent type, size_t len) : type(type), length(len) {}
    CigarPair(char type, size_t len) : type(CigarEventFromChar(type)), length(len) {}
    CigarPair Reverse() const {
        switch(type) {
            case M:
                return *this;
            case I:
                return {D, length};
            case D:
                return {I, length};
            default:
                VERIFY(false);
        }
    }
    size_t qlen() const {
        if(type == D)
            return 0;
        return length;
    }
    size_t tlen() const {
        if(type == I)
            return 0;
        return length;
    }
};

inline std::vector<CigarPair> RcCigar(const std::vector<CigarPair> &cigar) {
    return {cigar.rbegin(), cigar.rend()};
}

inline std::string CigarToString(const std::vector<CigarPair> &cigar) {
    std::stringstream ss;
    for(const CigarPair &p : cigar) {
        if(p.length != 1)
            ss << p.length;
        ss << p.type;
    }
    return ss.str();
}

inline std::vector<CigarPair> StringToCigar(const std::string& s) {
    std::vector<CigarPair> cigar;
    size_t n = 0;
    for(char c : s) {
        if(c >= '0' && c <= '9')
            n = n * 10 + c - '0';
        else {
            if(n == 0)
                n = 1;
            if(c == CigarEvent::M || c == CigarEvent::I || c == CigarEvent::D) {
                cigar.emplace_back(c, n);
            }
            n = 0;
        }
    }
    return std::move(cigar);
}

inline std::pair<size_t, size_t> LRskip(const std::string &s) {
    size_t lskip = 0;
    size_t rskip = 0;
    size_t n = 0;
    for(size_t i = 0; i < s.size(); i++) {
        char c = s[i];
        if(c >= '0' && c <= '9')
            n = n * 10 + c - '0';
        else {
            if(n == 0)
                n = 1;
            if(c == 'S' || c == 'H') {
                if(i == s.size() - 1)
                    rskip = n;
                else
                    lskip = n;
            }
            n = 0;
        }
    }
    return {lskip, rskip};
}


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

    class ConstAlignmentColumnIterator {
        const AlignmentForm *alignmentForm;
        size_t cigar_pos;
        size_t block_pos;
        size_t qpos;
        size_t tpos;
    public:
        ConstAlignmentColumnIterator(const AlignmentForm &form, size_t cigar_pos, size_t block_pos, size_t qpos, size_t tpos) :
                alignmentForm(&form), cigar_pos(cigar_pos), block_pos(block_pos), qpos(qpos), tpos(tpos) {
        }
        ConstAlignmentColumnIterator(const AlignmentForm &form, size_t cigar_pos, size_t block_pos);

        ConstAlignmentColumnIterator &operator++();
        ConstAlignmentColumnIterator operator++(int) const;
        ConstAlignmentColumnIterator &operator--();
        ConstAlignmentColumnIterator operator--(int) const;

        AlignmentColumn operator*() const {return {qpos, tpos, alignmentForm->cigar[cigar_pos].type};}

        bool operator==(const ConstAlignmentColumnIterator &other) const {return cigar_pos == other.cigar_pos && block_pos == other.block_pos;}
        bool operator!=(const ConstAlignmentColumnIterator &other) const {return !(*this == other);}
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
    IterableStorage<ConstAlignmentColumnIterator> columns() const {
        return {{*this, 0, 0, 0, 0}, {*this, cigar.size(), 0, qlen, tlen}};
    }
    AlignmentColumnIterator columnByQpos(size_t qpos);
    AlignmentColumnIterator columnByTpos(size_t tpos);


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


    std::string toCigarString() const;

    template<class U, class V>
    std::vector<std::string> toString(const U &from_seq, const V &to_seq,
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
    template<class U, class V>
    std::string toShortRec(const U &from_seq, const V &to_seq) const {
        std::vector<char> res;
        size_t m= 0;
        size_t mm = 0;
        size_t i = 0;
        size_t d = 0;
        size_t cnt = 0;
        for(const AlignmentColumn it : columns()) {
            if(it.event == CigarEvent::M) {
                if(from_seq[it.qpos] == to_seq[it.tpos])
                    m++;
                else
                    mm++;
            } else if(it.event == CigarEvent::I)
                i++;
            else
                d++;
            cnt++;
            if(cnt % 1000 == 0) {
               if(mm + i + d < 10) {
                   res.emplace_back(char('0' + mm + i + d));
               } else if (i >= 30 && d + mm <= 15) {
                   res.emplace_back('I');
               } else if(d >= 30 && i + mm <= 15) {
                   res.emplace_back('D');
               } else
                   res.emplace_back('*');
               m = 0;
               mm = 0;
               i = 0;
               d = 0;
            }
        }
        if(mm + i + d < 10) {
            res.emplace_back(char('0' + mm + i + d));
        } else if (i >= 30 && d + mm <= 15) {
            res.emplace_back('I');
        } else if(d >= 30 && i + mm <= 15) {
            res.emplace_back('D');
        } else
            res.emplace_back('*');
        return {res.begin(), res.end()};
    }
};

class AlignmentHelper {
public:
    template<class U, class V>
    static AlignmentForm::AlignmentColumnIterator
    LastComplexLongMatch(U tseq, V qseq, AlignmentForm &extension, size_t match) {
        size_t cnt = 0;
        size_t nucl_count[4] = {0,0,0,0};
        size_t max2 = match * 3 / 4;
        auto iter = extension.columns().end();
        while(iter != extension.columns().begin()) {
            --iter;
            auto column = *iter;
            if(column.event == M && tseq[column.tpos] == qseq[column.qpos]) {
                nucl_count[dignucl(tseq[column.tpos])]++;
                cnt++;
            } else {
                cnt = 0;
                nucl_count[0] = 0;
                nucl_count[1] = 0;
                nucl_count[2] = 0;
                nucl_count[3] = 0;
            }
            if(cnt >= match) {
                if(nucl_count[0] + nucl_count[1] <= max2 && nucl_count[0] + nucl_count[2] <= max2 && nucl_count[0] + nucl_count[3] <= max2 &&
                        nucl_count[2] + nucl_count[1] <= max2 && nucl_count[3] + nucl_count[1] <= max2 && nucl_count[2] + nucl_count[3] <= max2) {
                    for (size_t i = 0; i < match; i++)
                        ++iter;
                    break;
                } else {
                    VERIFY_MSG(tseq[column.tpos + match - 1] > 0, "Incorrect nucleotide counting");
                    nucl_count[dignucl(tseq[column.tpos + match - 1])]--;
                }
            }
        }
        return iter;
    }
};
