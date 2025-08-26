#pragma once
#include "common/fire_listeners.hpp"
namespace ag {
    class AlignedReadStorageFire;

    class AlignedReadStorageListener : public AbstractListener {
    public:
        explicit AlignedReadStorageListener(AlignedReadStorageFire &fire, const std::string &name);
        AlignedReadStorageListener(AlignedReadStorageListener &&)  noexcept = default;
        AlignedReadStorageListener &operator=(AlignedReadStorageListener &&)  noexcept = default;
        virtual void fireAddRead(const AlignedRead &read) {}
        virtual void fireDelayedRerouteRead(AlignedRead &read, const std::string &message) {}
        virtual void fireRerouteRead(AlignedRead &read) {}
        virtual void fireDelayedInvalidateRead(AlignedRead &read, const std::string &message) {}
        virtual void fireInvalidateRead(AlignedRead &read) {}
        virtual void fireAppliedCorrections(size_t cnt) {}
        virtual bool fireCheckConsistency() {return true;}
    };


    class AlignedReadStorageFire : public AbstractFire {
    public:
        AlignedReadStorageFire() = default;
        AlignedReadStorageFire(AlignedReadStorageFire &&)  noexcept = default;
        AlignedReadStorageFire &operator=(AlignedReadStorageFire &&)  noexcept = default;

        void fireAddRead(const AlignedRead &read);
        void fireDelayedRerouteRead(AlignedRead &read, const std::string &message);
        void fireRerouteRead(AlignedRead &read);
        void fireDelayedInvalidateRead(AlignedRead &read, const std::string &message);
        void fireInvalidateRead(AlignedRead &read);
        void fireAppliedCorrections(size_t cnt);
        bool fireCheckConsistency();
    };

}