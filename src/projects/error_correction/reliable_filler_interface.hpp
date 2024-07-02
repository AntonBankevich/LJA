#pragma once

#include <assembly_graph/component.hpp>
#include <dbg/visualization.hpp>
#include "assembly_graph/paths.hpp"
#include "dbg/sparse_dbg.hpp"
#include "correction_utils.hpp"
namespace dbg {
    class AbstractReliableFillingAlgorithm {
    public:
        virtual ~AbstractReliableFillingAlgorithm() = default;

        virtual size_t Fill(dbg::SparseDBG &dbg) = 0;

        virtual size_t LoggedFill(logging::Logger &logger, SparseDBG &dbg);

        virtual std::string name() const = 0;

        size_t ReFill(SparseDBG &dbg);

        size_t LoggedReFill(logging::Logger &logger, SparseDBG &dbg);
    };

    class CompositeReliableFiller : public AbstractReliableFillingAlgorithm {
        std::vector<AbstractReliableFillingAlgorithm *> &algorithms;
        std::string _name;
    public:
        CompositeReliableFiller(std::vector<AbstractReliableFillingAlgorithm *> &&algorithms);

        std::string name() const override { return _name; }

        size_t Fill(SparseDBG &dbg) override;

        virtual size_t LoggedFill(logging::Logger &logger, SparseDBG &dbg) {
            logger.info() << "Running multiple algorithms for reliable edge marking" << std::endl;
            size_t res = 0;
            for (AbstractReliableFillingAlgorithm *alg: algorithms) {
                size_t tmp = alg->Fill(dbg);
                logger.info() << alg->name() << " marked " << tmp << " reliable edges" << std::endl;
                res += tmp;
            }
            logger.info() << "In total " << res << " new reliable edges were marked" << std::endl;
            return res;
        }
    };

}