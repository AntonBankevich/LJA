#pragma once

#include <dbg/paths.hpp>
#include "dbg/sparse_dbg.hpp"

class AbstractReliableFillingAlgorithm {
public:
    virtual size_t Fill(dbg::SparseDBG &dbg) = 0;
    virtual size_t LoggedFill(logging::Logger &logger, dbg::SparseDBG &dbg);
    virtual std::string name() const = 0;
    size_t ReFill(dbg::SparseDBG &dbg);
    size_t LoggedReFill(logging::Logger &logger, dbg::SparseDBG &dbg);
};

class CompositeReliableFiller : public AbstractReliableFillingAlgorithm {
    std::vector<AbstractReliableFillingAlgorithm *> &algorithms;
    std::string _name;
public:
    CompositeReliableFiller(std::vector<AbstractReliableFillingAlgorithm *> &&algorithms);
    std::string name() const override {return _name;}
    size_t Fill(dbg::SparseDBG &dbg) override;
    virtual size_t LoggedFill(logging::Logger &logger, dbg::SparseDBG &dbg) {
        logger.info() << "Running multiple algorithms for reliable edge marking" << std::endl;
        size_t res = 0;
        for(AbstractReliableFillingAlgorithm *alg : algorithms) {
            size_t tmp = alg->Fill(dbg);
            logger.info() << alg->name() << " marked " << tmp << " reliable edges" << std::endl;
            res += tmp;
        }
        logger.info() << "In total " << res << " new reliable edges were marked" << std::endl;
        return res;
    }
};

class CoverageReliableFiller : public AbstractReliableFillingAlgorithm {
private:
    double threshold;
public:
    CoverageReliableFiller(double threshold) :threshold(threshold) {}
    std::string name() const override {return "CoverageReliableFiller";}
    size_t Fill(dbg::SparseDBG &sdbg) override;
};

class LengthReliableFiller : public AbstractReliableFillingAlgorithm {
private:
    size_t min_length;
    double min_rel_cov;
    double max_err_cov;
public:
    LengthReliableFiller(size_t min_length, double min_rel_cov, double max_err_cov) : min_length(min_length), min_rel_cov(min_rel_cov), max_err_cov(max_err_cov) {}
    std::string name() const override {return "LengthReliableFiller";}
    size_t Fill(dbg::SparseDBG &dbg) override;
};

class BridgeReliableFiller : public AbstractReliableFillingAlgorithm {
private:
    size_t max_length;
    std::vector<dbg::Edge *> bridges(dbg::Edge &start);
public:
    BridgeReliableFiller(size_t max_length) : max_length(max_length) {}
    std::string name() const override {return "BridgeReliableFiller";}
    size_t Fill(dbg::SparseDBG &dbg) override;
};

class ConnectionReliableFiller : public AbstractReliableFillingAlgorithm {
private:
    double threshold;
    dbg::Edge *checkBorder(dbg::Vertex &v);
    bool checkInner(dbg::Vertex &v);
public:
    ConnectionReliableFiller(double threshold) : threshold(threshold){}
    std::string name() const override {return "ConnectionReliableFiller";}
    size_t Fill(dbg::SparseDBG &dbg) override;
};


