#include "graph_listeners.hpp"
#include "assembly_graph.hpp"
using namespace ag;
ResolutionListener::ResolutionListener(ResolutionFire &fire, const std::string &name) : AbstractListener(fire, name) {}

AssemblyGraph &ResolutionListener::getGraph() {
    return static_cast<AssemblyGraph &>(getFire<ResolutionFire>());
}

//    The mirrored call would see the same vertices again if P is its own rc or a removed vertex is.
void ResolutionFire::fireReplacePath(const RAGraphPath &old_path, const RAGraphPath &new_path, size_t dominant_index,
                                     const std::unordered_map<EdgeId, EdgeId> &outer) {
    PathReplacement forward(old_path, new_path, dominant_index, outer);
    RAGraphPath rc_old = old_path.RC();
    RAGraphPath rc_new = new_path.RC();
    bool mirrored = !(rc_old == old_path);
    for (size_t i = 0; i <= old_path.size(); i++)
        mirrored = mirrored && (forward.survives(i) || old_path.getVertex(i) != old_path.getVertex(i).rc());
    std::unique_ptr<PathReplacement> backward;
    if (mirrored)
        backward = std::make_unique<PathReplacement>(rc_old, rc_new, new_path.size() - dominant_index, outer);
    for (ResolutionListener *listener: getListeners<ResolutionListener>()) {
        listener->fireReplacePath(forward);
        if (backward)
            listener->fireReplacePath(*backward);
    }
}
