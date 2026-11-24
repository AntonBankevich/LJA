#include "compression.hpp"
#include "common/string_utils.hpp"
#include <vector>

CompressionParameters CompressionParameters::global;

void CompressionParameters::setDimers(const std::string &s) {
    std::vector<std::string> vals = split(s, ",");
    VERIFY_MSG(vals.size() == 3, "Dimer compression parameters must be <min>,<max>,<step>: " << s);
    min_dimer_to_compress = std::stoull(vals[0]);
    max_dimer_size = std::stoull(vals[1]);
    dimer_step = std::stoull(vals[2]);
}
