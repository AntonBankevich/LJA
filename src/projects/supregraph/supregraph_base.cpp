#include "supregraph_base.hpp"

bool spg::SPGVertex::isCore() {
    for(Edge &edge: *this) {
        if(edge.isSuffix())
            return false;
    }
    for(Edge &edge: rc()) {
        if(edge.isSuffix())
            return false;
    }
    return true;
}