#include "supregraph.hpp"
#include "multiplexer.hpp"
#include "read_storage.hpp"

using namespace spg;
int main() {
    SupreGraph g;
    PathStorage storage(g);
    AndreyRule rule(storage);
    spg::Multiplexer multiplexer(g, rule);
    multiplexer.fullMultiplex();
}