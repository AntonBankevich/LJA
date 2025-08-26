#include "fire_listeners.hpp"
#include <algorithm>

AbstractFire::AbstractFire(AbstractFire &&other) noexcept {
    *this = std::move(other);
}

AbstractFire &AbstractFire::operator=(AbstractFire &&other) {
    std::swap(listeners, other.listeners);
    for(AbstractListener *listener : listeners) {
        listener->fire = this;
    }
    for(AbstractListener *listener : other.listeners)
        listener->fire = &other;
    return *this;
}

void AbstractFire::replaceListener(AbstractListener &old_listener, AbstractListener &new_listener) { // NOLINT(readability-convert-member-functions-to-static)
    auto it = std::find(listeners.begin(), listeners.end(), &old_listener);
    VERIFY(it != listeners.end());
    *it = &new_listener;
}

void AbstractFire::removeListener(AbstractListener &listener) {
    listeners.erase(std::find(listeners.begin(), listeners.end(), &listener));
}

AbstractListener::AbstractListener(AbstractListener &&other) noexcept: fire(nullptr) {
    *this = std::move(other);
}

AbstractListener::AbstractListener(AbstractFire &fire, const std::string &name) : fire(&fire), name(name) {attach();}

void AbstractListener::attach() {
    VERIFY(fire != nullptr);
    fire->addListener(*this);
}

void AbstractListener::detach() {
    if (fire != nullptr)
        fire->removeListener(*this);
    fire = nullptr;
}

bool AbstractListener::active() const {return fire != nullptr;}

AbstractListener &AbstractListener::operator=(AbstractListener &&other) noexcept {
    if(fire != nullptr)
        fire->replaceListener(*this, other);
    if(other.fire != nullptr)
        other.fire->replaceListener(other, *this);
    std::swap(fire, other.fire);
    std::swap(name, other.name);
    return *this;
}

AbstractListener::~AbstractListener() {
    detach();
}
