#pragma once
#include "error_correction.hpp"
class Precorrector : public AbstractCorrectionAlgorithm {
private:
    double reliable_threshold;
public:
    Precorrector(double reliable_threshold) :
            AbstractCorrectionAlgorithm("Precorrector"), reliable_threshold(reliable_threshold) {}

    std::string correctRead(dbg::GraphAlignment &path) override;
};
