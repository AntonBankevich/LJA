#include <assembly_graph/aligned_read.hpp>
#include "aligned_read_listeners.hpp"
using namespace ag;

AlignedReadStorageListener::AlignedReadStorageListener(AlignedReadStorageFire &fire, const std::string &name) : AbstractListener(fire, name) {}

void AlignedReadStorageFire::fireAddRead(const AlignedRead &read) {
    for (AlignedReadStorageListener *listener: getListeners<AlignedReadStorageListener>()) {
        listener->fireAddRead(read);
    }
}

void AlignedReadStorageFire::fireDelayedRerouteRead(AlignedRead &read, const std::string &message) {
    for (AlignedReadStorageListener *listener: getListeners<AlignedReadStorageListener>()) {
        listener->fireDelayedRerouteRead(read, message);
    }
}

void AlignedReadStorageFire::fireRerouteRead(AlignedRead &read) {
    for (AlignedReadStorageListener *listener: getListeners<AlignedReadStorageListener>()) {
        listener->fireRerouteRead(read);
    }
}

void AlignedReadStorageFire::fireDelayedInvalidateRead(AlignedRead &read, const std::string &message) {
    for (AlignedReadStorageListener *listener: getListeners<AlignedReadStorageListener>()) {
        listener->fireDelayedInvalidateRead(read, message);
    }
}

void AlignedReadStorageFire::fireInvalidateRead(AlignedRead &read) {
    for (AlignedReadStorageListener *listener: getListeners<AlignedReadStorageListener>()) {
        listener->fireInvalidateRead(read);
    }
}

void AlignedReadStorageFire::fireAppliedCorrections(size_t cnt) {
    for (AlignedReadStorageListener *listener: getListeners<AlignedReadStorageListener>()) {
        listener->fireAppliedCorrections(cnt);
    }
}

bool AlignedReadStorageFire::fireCheckConsistency() {
    for (AlignedReadStorageListener *listener: getListeners<AlignedReadStorageListener>()) {
        if(!listener->fireCheckConsistency())
            return false;
    }
    return true;
}
