#pragma once
//
// Created by anton on 7/20/20.
//

#include "common/hash_utils.hpp"
#include "sequences/sequence.hpp"
#include "iterator_utils.hpp"
#include <deque>
#include <utility>

namespace hashing {
    template<typename T, typename U>
    T pow(T base, U p) {
        if (p == 0)
            return 1;
        T tmp = pow(base, p / 2);
        if (p % 2 == 1)
            return base * tmp * tmp;
        else
            return tmp * tmp;
    }

    class MovingKWH;
    class KWHIterator;
    class RollingHash {
    private:
        size_t k;
        htype hbase;
        htype kpow;
        htype inv;
        static const htype HBASE = 239;
    public:

        explicit RollingHash(size_t _k, htype _hbase = HBASE) : k(_k), hbase(_hbase),
                                               kpow(pow(hbase, k - 1)),
                                               inv(pow(hbase, (htype(1u) << (sizeof(htype) * 8u - 1u)) - 1u)) {
            VERIFY(inv * hbase == htype(1));
        }
        RollingHash extensionHash() const {return RollingHash(k + 1, hbase);}

        size_t getK() const {return k;}
        htype hash(const Sequence &seq, size_t pos) const;
        htype extendRight(const Sequence &seq, htype hash, unsigned char c) const {return hash * hbase + c;}
        htype extendLeft(const Sequence &seq, htype hash, unsigned char c) const {return hash + c * kpow * hbase;}
        htype shiftRight(htype hash, unsigned char first_letter, unsigned char next_letter) const;
        htype shiftLeft(htype hash, unsigned char prev_letter, unsigned char last_letter) const;
        htype next(const Sequence &seq, size_t pos, htype hash) const;
        htype prev(const Sequence &seq, size_t pos, htype hash) const;

        IterableStorage<KWHIterator> kmers(const Sequence &seq, size_t from = 0, size_t to = size_t(-1)) const &;
        IterableStorage<KWHIterator> innerKmers(const Sequence &seq) const &;
    };

//    KeyWithHash (KWH) represents a record of a k-mer with calculated hash and reverse hash.
// Two KWH are equal iff the k-mers they represent are equal. To avoid long comparison calculation we only compare hashes for equality.
    class KWH {
    protected:
        htype fhash;
        htype rhash;
        Sequence seq;
        size_t pos;
    public:
        KWH(Sequence _seq, size_t _pos, htype _fhash, htype _rhash) :
                seq(std::move(_seq)), pos(_pos), fhash(_fhash), rhash(_rhash) {
        }
        size_t getPos() const {return pos;}
        htype hash() const {return std::min(fhash, rhash);}
        htype fHash() const {return fhash;}
        htype rHash() const {return rhash;}
        Sequence getSeq(size_t k) const {return seq.Subseq(getPos(), getPos() + k);}

        bool isFirst() const {return pos == 0;}

        bool operator==(const KWH &other) const {return fhash == other.fhash && rhash == other.rhash;}
        bool operator<(const KWH &other) const;
        bool operator>(const KWH &other) const;
        bool operator<=(const KWH &other) const;
        bool operator>=(const KWH &other) const;

        bool isCanonical() const {return fhash < rhash;}
    };

    class KWHIterator;

    class MovingKWH : public KWH {
        friend class KWHIterator;
    private:
        const RollingHash *hasher;
        MovingKWH(const RollingHash &_hasher, Sequence _seq, size_t _pos, htype _fhash, htype _rhash) :
                hasher(&_hasher), KWH(std::move(_seq), _pos, _fhash, _rhash) {
        }

        //        These two methods avoid copy of Sequence, which contains updates of a synchronized counter
        MovingKWH &moveForward();
        MovingKWH &moveBackward();
    public:
        MovingKWH(const RollingHash &_hasher, const Sequence &_seq, size_t _pos) :
                hasher(&_hasher), KWH(_seq, _pos, _hasher.hash(_seq, _pos), _hasher.hash(!_seq, _seq.size() - _pos - _hasher.getK())) {
        }

        MovingKWH(const MovingKWH &other) = default;
        MovingKWH &operator=(const MovingKWH &other) = default;
        MovingKWH(MovingKWH &&other) = default;
        MovingKWH &operator=(MovingKWH &&other) = default;

        MovingKWH operator!() const;


        Sequence getSeq() const {return KWH::getSeq(hasher->getK());}

        bool isLast() const {return pos + hasher->getK() == seq.size();}

        htype extendRight(unsigned char c) const;
        htype extendLeft(unsigned char c) const;
        bool hasNext() const {return pos + hasher->getK() < seq.size();}
        bool hasPrev() const {return pos > 0;}
        bool isValid() const {return pos + hasher->getK() <= seq.size();}

        MovingKWH next() const;
        MovingKWH prev() const;
    };

    class KWHIterator {
    private:
        bool is_end;
        MovingKWH kwh;
    public:
        KWHIterator(const RollingHash &hasher, Sequence _seq, size_t pos) :
                is_end(pos == _seq.size() - hasher.getK() + 1), kwh(hasher, std::move(_seq), is_end ? pos - 1 : pos) {
        }

        bool isEnd() const {return is_end;}

        KWHIterator &operator++();
        KWHIterator &operator--();

        const MovingKWH &operator*() const;
        const MovingKWH *operator->() const {return &kwh;}

        bool operator==(const KWHIterator &other) {return is_end == other.is_end && kwh.getPos() == other.kwh.getPos();}
        bool operator!=(const KWHIterator &other) {return is_end != other.is_end || kwh.getPos() != other.kwh.getPos();}
    };


    class MinQueue {
        std::deque<std::pair<hashing::htype, size_t>> q;
    public:
        MinQueue() = default;

        bool empty() const {return q.empty();}
        const std::pair<hashing::htype, size_t> &get() const {return q.front();}
        size_t size() const {return q.size();}

        void push(const KWH &kwh);
        void pop(size_t pos);
    };

    class MinimizerCalculator {
    private:
        const size_t w;
        KWHIterator it;
        MinQueue queue;
    public:
        MinimizerCalculator(Sequence seq, const RollingHash &_hasher, size_t _w);
        std::pair<hashing::htype, size_t> next();
        bool hasNext() const {return !it.isEnd();}
        std::vector<htype> minimizerHashs();
    };
}