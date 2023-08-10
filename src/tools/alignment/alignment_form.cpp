#include "alignment_form.hpp"

void AlignmentForm::operator+=(const AlignmentForm &other) {
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

void AlignmentForm::operator+=(CigarEvent e) {
    if(e == cigar.back().type) {
        cigar.back().length ++;
    } else {
        cigar.emplace_back(e, 1);
    }
    if(e != CigarEvent::D)
        qlen += 1;
    if(e != CigarEvent::I)
        tlen += 1;
}

void AlignmentForm::operator+=(CigarPair p) {
    if(p.type == cigar.back().type) {
        cigar.back().length += p.length;
    } else {
        cigar.emplace_back(p);
    }
    if(p.type != CigarEvent::D)
        qlen += p.length;
    if(p.type != CigarEvent::I)
        tlen += p.length;
}

AlignmentForm AlignmentForm::operator+(const AlignmentForm &other) const {
    AlignmentForm res = *this;
    res += other;
    return std::move(res);
}

AlignmentForm AlignmentForm::RC() const {
    return {std::vector<CigarPair>(cigar.rbegin(), cigar.rend())};
}

AlignmentForm AlignmentForm::Reverse() const {
    AlignmentForm res;
    for(CigarPair p : cigar) {
        res += p.Reverse();
    }
    return std::move(res);
}

AlignmentForm AlignmentForm::Prefix(AlignmentForm::AlignmentColumnIterator bound) {
    AlignmentForm res;
    for(auto it = columns().begin(); it != bound; ++it)
        res += (*it).event;
    return std::move(res);
}

AlignmentForm AlignmentForm::Suffix(AlignmentForm::AlignmentColumnIterator bound) {
    AlignmentForm res;
    for(auto it = bound; it != columns().end(); ++it)
        res += (*it).event;
    return std::move(res);
}

void AlignmentForm::calculateLens() {
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

AlignmentForm::AlignmentForm(const std::string &s) {
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

std::string AlignmentForm::toCigarString() const {
    std::stringstream ss;
    for(auto p : *this) {
        if(p.length == 1)
            ss << p.type;
        else
            ss << p.length << p.type;
    }
    return ss.str();
}

AlignmentForm::AlignmentColumnIterator::AlignmentColumnIterator(AlignmentForm &form, size_t cigar_pos, size_t block_pos)
        :
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

AlignmentForm::AlignmentColumnIterator &AlignmentForm::AlignmentColumnIterator::operator++() {
    this->qpos += CigarPair(alignmentForm->cigar[cigar_pos].type, 1).qlen();
    this->tpos += CigarPair(alignmentForm->cigar[cigar_pos].type, 1).tlen();
    block_pos++;
    if(alignmentForm->cigar[cigar_pos].length == block_pos) {
        cigar_pos++;
        block_pos = 0;
    }
    return *this;
}

AlignmentForm::AlignmentColumnIterator AlignmentForm::AlignmentColumnIterator::operator++(int) const {
    AlignmentColumnIterator res = *this;
    ++res;
    return res;
}

AlignmentForm::AlignmentColumnIterator &AlignmentForm::AlignmentColumnIterator::operator--() {
    if(block_pos == 0) {
        cigar_pos--;
        block_pos = alignmentForm->cigar[cigar_pos].length;
    }
    block_pos--;
    this->qpos -= CigarPair(alignmentForm->cigar[cigar_pos].type, 1).qlen();
    this->tpos -= CigarPair(alignmentForm->cigar[cigar_pos].type, 1).tlen();
    return *this;
}

AlignmentForm::AlignmentColumnIterator AlignmentForm::AlignmentColumnIterator::operator--(int) const {
    AlignmentColumnIterator res = *this;
    --res;
    return res;
}
