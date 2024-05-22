#include "assembly_graph/assembly_graph.hpp"
#include "supregraph.hpp"
#include "listeners.hpp"
using namespace spg;
ResolutionListener::ResolutionListener(ResolutionFire &fire) : fire(&fire) {
    this->fire->addListener(*this);
}

ResolutionListener::~ResolutionListener() {
    if(fire != nullptr)
        fire->removeListener(*this);
    fire = nullptr;
}

ResolutionListener::ResolutionListener(ResolutionListener &&other) noexcept {
    fire = other.fire;
    other.fire = nullptr;
    fire->replaceListener(other, *this);
}