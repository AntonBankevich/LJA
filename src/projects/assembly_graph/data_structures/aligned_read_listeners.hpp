#pragma once
#include "common/fire_listeners.hpp"
#include "assembly_graph/aligned_read.hpp"
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
//        Checks the listener's state against the stored reads, reporting every violation to os.
        virtual bool fireCheckConsistency(std::ostream &os) {return true;}
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
//        Runs every listener's check (none is skipped after a failure).
        bool fireCheckConsistency(std::ostream &os);
    };

}