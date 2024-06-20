#pragma once
#include "common/fire_listeners.hpp"
namespace ag {
    template<class Traits>
    class AlignedReadStorageFire;

    template<class Traits>
    class AlignedReadStorageListener : public AbstractListener {
    public:
        explicit AlignedReadStorageListener(AlignedReadStorageFire<Traits> &fire, const std::string &name) : AbstractListener(fire, name) {}
        AlignedReadStorageListener(AlignedReadStorageListener &&)  noexcept = default;
        AlignedReadStorageListener &operator=(AlignedReadStorageListener &&)  noexcept = default;
        virtual void fireAddRead(const AlignedRead<Traits> &read) {}
        virtual void fireDelayedRerouteRead(AlignedRead<Traits> &read, const std::string &message) {}
        virtual void fireRerouteRead(AlignedRead<Traits> &read) {}
        virtual void fireDelayedInvalidateRead(AlignedRead<Traits> &read, const std::string &message) {}
        virtual void fireInvalidateRead(AlignedRead<Traits> &read) {}
        virtual void fireAppliedCorrections(size_t cnt) {}
        virtual bool fireCheckConsistency() {return true;}
    };

    template<class Traits>
    class AlignedReadStorageFire : public AbstractFire {
    public:
        AlignedReadStorageFire() = default;
        AlignedReadStorageFire(AlignedReadStorageFire &&)  noexcept = default;
        AlignedReadStorageFire &operator=(AlignedReadStorageFire &&)  noexcept = default;

        void fireAddRead(const AlignedRead<Traits> &read) {
            for (AlignedReadStorageListener<Traits> *listener: getListeners<AlignedReadStorageListener<Traits>>()) {
                listener->fireAddRead(read);
            }
        }

        void fireDelayedRerouteRead(AlignedRead<Traits> &read, const std::string &message) {
            for (AlignedReadStorageListener<Traits> *listener: getListeners<AlignedReadStorageListener<Traits>>()) {
                listener->fireDelayedRerouteRead(read, message);
            }
        }

        void fireRerouteRead(AlignedRead<Traits> &read) {
            for (AlignedReadStorageListener<Traits> *listener: getListeners<AlignedReadStorageListener<Traits>>()) {
                listener->fireRerouteRead(read);
            }
        }

        void fireDelayedInvalidateRead(AlignedRead<Traits> &read, const std::string &message) {
            for (AlignedReadStorageListener<Traits> *listener: getListeners<AlignedReadStorageListener<Traits>>()) {
                listener->fireDelayedInvalidateRead(read, message);
            }
        }

        void fireInvalidateRead(AlignedRead<Traits> &read) {
            for (AlignedReadStorageListener<Traits> *listener: getListeners<AlignedReadStorageListener<Traits>>()) {
                listener->fireInvalidateRead(read);
            }
        }

        void fireAppliedCorrections(size_t cnt) {
            for (AlignedReadStorageListener<Traits> *listener: getListeners<AlignedReadStorageListener<Traits>>()) {
                listener->fireAppliedCorrections(cnt);
            }
        }

        bool fireCheckConsistency() {
            for (AlignedReadStorageListener<Traits> *listener: getListeners<AlignedReadStorageListener<Traits>>()) {
                if(!listener->fireCheckConsistency())
                    return false;
            }
            return true;
        }
    };
}