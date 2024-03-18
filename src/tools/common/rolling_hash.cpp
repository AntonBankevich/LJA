#include "rolling_hash.hpp"
#include "iterator_utils.hpp"
namespace hashing {
    IterableStorage <KWHIterator> RollingHash::kmers(const Sequence &seq, size_t from, size_t to) const & {
        return {{*this, seq, from},
                {*this, seq, to == size_t(-1) ? seq.size() - getK() + 1 : to}};
    }

    IterableStorage <KWHIterator> RollingHash::innerKmers(const Sequence &seq) const & {
        return kmers(seq, 1, seq.size() - getK());
    }

    htype RollingHash::hash(const Sequence &seq, size_t pos) const {
        htype hash = 0;
        for (size_t i = pos; i < pos + k; i++) {
            hash = hash * hbase + seq[i];
        }
        return hash;
    }

    htype RollingHash::shiftRight(htype hash, unsigned char first_letter, unsigned char next_letter) const {
        return (hash - kpow * first_letter) * hbase + next_letter;
    }

    htype RollingHash::shiftLeft(htype hash, unsigned char prev_letter, unsigned char last_letter) const {
        return (hash - last_letter) * inv + prev_letter * kpow;
    }

    htype RollingHash::next(const Sequence &seq, size_t pos, htype hash) const {
        return shiftRight(hash, seq[pos], seq[pos + k]);
    }

    htype RollingHash::prev(const Sequence &seq, size_t pos, htype hash) const {
        return shiftLeft(hash, seq[pos - 1], seq[pos + k - 1]);
    }

    std::pair<hashing::htype, size_t> MinimizerCalculator::next() {
        queue.pop((*it).getPos() - w + 1);
        queue.push(*it);
        ++it;
        return queue.get();
    }

    MinimizerCalculator::MinimizerCalculator(Sequence seq, const RollingHash &_hasher, size_t _w) :
            w(_w), it(_hasher, std::move(seq), 0) {
        VERIFY(w >= 1);
        VERIFY(seq.size() >= _hasher.getK() + w - 1)
        for (size_t i = 0; i < w; i++) {
            queue.push(*it);
            ++it;
        }
    }

    std::vector<htype> MinimizerCalculator::minimizerHashs() {
        std::vector<htype> res;
        res.push_back(queue.get().first);
        while (hasNext()) {
            htype val = next().first;
            if (val != res.back()) {
                res.push_back(val);
            }
        }
        return std::move(res);
    }

    bool KWH::operator<(const KWH &other) const {
        if(hash() != other.hash()) {return hash() < other.hash();}
        if(fhash != other.fhash) {return fhash < other.fhash;}
        return rhash < other.hash();
    }

    bool KWH::operator>(const KWH &other) const {
        if(hash() != other.hash()) {return hash() > other.hash();}
        if(fhash != other.fhash) {return fhash > other.fhash;}
        return rhash > other.hash();
    }

    bool KWH::operator<=(const KWH &other) const {
        if(hash() != other.hash()) {return hash() < other.hash();}
        if(fhash != other.fhash) {return fhash < other.fhash;}
        return rhash <= other.hash();
    }

    bool KWH::operator>=(const KWH &other) const {
        if(hash() != other.hash()) {return hash() > other.hash();}
        if(fhash != other.fhash) {return fhash > other.fhash;}
        return rhash >= other.hash();
    }

    KWHIterator &KWHIterator::operator++() {
        VERIFY(!is_end);
        if(kwh.hasNext())
            kwh.moveForward();
        else
            is_end = true;
        return *this;
    }

    KWHIterator &KWHIterator::operator--() {
        VERIFY(is_end || !kwh.isFirst());
        if(is_end)
            is_end = false;
        else
            kwh.moveBackward();
        return *this;
    }

    const MovingKWH &KWHIterator::operator*() const {
        VERIFY(!is_end);
        return kwh;
    }

    MovingKWH &MovingKWH::moveForward() {
        unsigned char first = seq[pos];
        unsigned char next = hasNext() ? seq[pos + hasher->getK()] : 0;
        fhash = hasher->shiftRight(fhash, first, next);
        rhash = hasher->shiftLeft(rhash, complement(next), complement(first));
        pos++;
        return *this;
    }

    MovingKWH &MovingKWH::moveBackward() {
        VERIFY(pos > 0);
        unsigned char prev = seq[pos - 1];
        unsigned char last = seq[pos + hasher->getK() - 1];
        fhash = hasher->shiftLeft(fhash, prev, last);
        rhash = hasher->shiftRight(rhash, complement(last), complement(prev));
        pos--;
        return *this;
    }

    MovingKWH MovingKWH::operator!() const {
        VERIFY(isValid());
        return {*hasher, !seq, seq.size() - getPos() - hasher->getK(), rhash, fhash};
    }

    htype MovingKWH::extendRight(unsigned char c) const {
        return std::min(hasher->extendRight(seq, fhash, c),
                        hasher->extendLeft(!seq, rhash, c ^ 3u));
    }

    htype MovingKWH::extendLeft(unsigned char c) const {
        return std::min(hasher->extendLeft(seq, fhash, c),
                        hasher->extendRight(!seq, rhash, c ^ 3u));
    }

    MovingKWH MovingKWH::next() const {
        VERIFY(hasNext());
        MovingKWH res(*this);
        res.moveForward();
        return std::move(res);
    }

    MovingKWH MovingKWH::prev() const {
        VERIFY(hasPrev());
        MovingKWH res(*this);
        res.moveBackward();
        return std::move(res);
    }

    void MinQueue::push(const KWH &kwh) {
        while (!q.empty() && q.back().first > kwh.hash()) {
            q.pop_back();
        }
        q.emplace_back(kwh.hash(), kwh.getPos());
    }

    void MinQueue::pop(size_t pos) {
        if (!q.empty() && q.front().second < pos) {
            q.pop_front();
        }
    }
}

