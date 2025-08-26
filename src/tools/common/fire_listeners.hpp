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
    template<class Listener>
    std::vector<Listener *>
//    IterableStorage<TransformingIterator<std::vector<AbstractListener *>::iterator, Listener>>
    getListeners() {
        std::function<Listener *(AbstractListener *&)> transform = [](AbstractListener *&listener) -> Listener *{
            return (static_cast<Listener *>(listener));
        };
        std::vector<Listener*> transformed = oneline::map(listeners.begin(), listeners.end(), transform);
        return std::move(transformed);
//        TransformingIterator<std::vector<AbstractListener *>::iterator, Listener> begin(listeners.begin(), listeners.end(), transform);
//        TransformingIterator<std::vector<AbstractListener *>::iterator, Listener> end(listeners.end(), listeners.end(), transform);
//        return {begin, end};
//        return {{listeners.begin(), listeners.end(), transform}, {listeners.end(), listeners.end(), transform}};
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
