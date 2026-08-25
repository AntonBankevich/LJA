#pragma once
#include "error_correction.hpp"
#include <functional>
class Precorrector : public ag::AbstractCorrectionAlgorithm {
private:
    std::function<bool(const ag::Edge &)> isReliable;
    std::function<bool(const ag::Edge &)> isSuspicious;
public:
    Precorrector(std::function<bool(const ag::Edge &)> isReliable, std::function<bool(const ag::Edge &)> isSuspicious) :
            ag::AbstractCorrectionAlgorithm("Precorrector"), isReliable(std::move(isReliable)), isSuspicious(std::move(isSuspicious)) {}

    std::string correctRead(const std::string &name, ag::GraphPath &path) override;
};
