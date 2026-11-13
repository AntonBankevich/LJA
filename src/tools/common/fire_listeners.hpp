#pragma once
#include "iterator_utils.hpp"
#include "common/oneline_utils.hpp"
#include "vector"

class AbstractFire;
class AbstractListener {
    friend class AbstractFire;
private:
    AbstractFire *fire = nullptr;
    std::string name;
public:
    explicit AbstractListener(AbstractFire &fire, const std::string &name);
    template<class T>
    T &getFire() {return *((T*)fire);}
    void attach();
    void detach();
    bool active() const;

    const std::string &getName() const {return name;}
    AbstractListener(AbstractListener &&other) noexcept;
    AbstractListener &operator=(AbstractListener &&other)  noexcept;
    AbstractListener(const AbstractListener &other) = delete;
    AbstractListener &operator=(const AbstractListener &other) = delete;
    virtual ~AbstractListener();
};

class AbstractFire {
    friend class AbstractListener;
private:
    std::vector<AbstractListener *> listeners;
protected:
    template<class Listener>
    Listener &getListener() {
        return *(static_cast<Listener *>(listeners.front()));
    }
//    The listeners cast to Listener, iterated in place (no copy): listeners must not attach or detach while an event
//    is being dispatched.
    template<class Listener>
    class ListenerRange {
        typedef std::vector<AbstractListener *>::const_iterator BaseIterator;
        BaseIterator b, e;
    public:
        class Iterator {
            BaseIterator it;
        public:
            explicit Iterator(BaseIterator it) : it(it) {}
            Listener *operator*() const {return static_cast<Listener *>(*it);}
            Iterator &operator++() {++it; return *this;}
            bool operator!=(const Iterator &other) const {return it != other.it;}
        };
        ListenerRange(BaseIterator b, BaseIterator e) : b(b), e(e) {}
        Iterator begin() const {return Iterator(b);}
        Iterator end() const {return Iterator(e);}
    };
    template<class Listener>
    ListenerRange<Listener> getListeners() const {
        return {listeners.begin(), listeners.end()};
    }
public:
    AbstractFire() = default;
    AbstractFire(AbstractFire &&other) noexcept ;
    AbstractFire &operator=(AbstractFire &&other);
    virtual ~AbstractFire() {VERIFY(listeners.empty());}
    void addListener(AbstractListener &listener) {listeners.emplace_back(&listener);}
//        Used when listeners are move-assigned or move-copied
    void replaceListener(AbstractListener &old_listener, AbstractListener &new_listener);
    void removeListener(AbstractListener &listener);
};
