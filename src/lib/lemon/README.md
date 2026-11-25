Header-only subset of LEMON 1.3.1 (http://lemon.cs.elte.hu, Boost Software License, see LICENSE), unchanged
except for the hand-written lemon/config.h. Only what NetworkSimplex over SmartDigraph needs is included; when
adding another LEMON algorithm, copy its header closure from the same release (g++ -MM lists it).
Include as <lemon/network_simplex.h>; src/lib is on the include path.
