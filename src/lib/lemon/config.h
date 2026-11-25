// Hand-written replacement for LEMON's CMake-generated config.h (LEMON 1.3.1, vendored header-only subset).
// LEMON_ONLY_TEMPLATES makes core.h define INVALID inline, so no LEMON .cc file has to be compiled.
#pragma once
#define LEMON_VERSION "1.3.1"
#define LEMON_HAVE_LONG_LONG 1
#define LEMON_ONLY_TEMPLATES 1
