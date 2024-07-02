#pragma once

#include "bulge_path_marker.hpp"
#include "reliable_filler_interface.hpp"
#include "assembly_graph/paths.hpp"
#include "dbg/sparse_dbg.hpp"
#ifdef USE_LIBTORCH
#include "torch/torch.h"
#endif // USE_LIBTORCH

namespace dbg {
    class CoverageReliableFiller : public AbstractReliableFillingAlgorithm {
    private:
        double threshold;
    public:
        CoverageReliableFiller(double threshold) : threshold(threshold) {}

        std::string name() const override { return "CoverageReliableFiller"; }

        size_t Fill(dbg::SparseDBG &sdbg) override;
    };

    class LengthReliableFiller : public AbstractReliableFillingAlgorithm {
    private:
        size_t min_length;
        double min_rel_cov;
        double max_err_cov;
    public:
        LengthReliableFiller(size_t min_length, double min_rel_cov, double max_err_cov) : min_length(min_length),
                                                                                          min_rel_cov(min_rel_cov),
                                                                                          max_err_cov(max_err_cov) {}

        std::string name() const override { return "LengthReliableFiller"; }

        size_t Fill(dbg::SparseDBG &dbg) override;
    };

    class BridgeReliableFiller : public AbstractReliableFillingAlgorithm {
    private:
        size_t max_length;

        std::vector<dbg::Edge *> bridges(dbg::Edge &start);

    public:
        BridgeReliableFiller(size_t max_length) : max_length(max_length) {}

        std::string name() const override { return "BridgeReliableFiller"; }

        size_t Fill(dbg::SparseDBG &dbg) override;
    };

    class ConnectionReliableFiller : public AbstractReliableFillingAlgorithm {
    private:
        double threshold;

        dbg::Edge *checkBorder(dbg::Vertex &v);

        bool checkInner(dbg::Vertex &v);

    public:
        ConnectionReliableFiller(double threshold) : threshold(threshold) {}

        std::string name() const override { return "ConnectionReliableFiller"; }

        size_t Fill(dbg::SparseDBG &dbg) override;
    };

    inline CompositeReliableFiller
    CreateDefaultReliableFiller(dbg::SparseDBG &dbg, dbg::ReadAlignmentStorage &reads, double reliable_threshold,
                                bool diploid) {
        CoverageReliableFiller cov(reliable_threshold);
        LengthReliableFiller len(20000, 3, 1);
        BridgeReliableFiller bridge(40000);
        ConnectionReliableFiller connect(reliable_threshold);
        dbg::BulgePathMarker bulge(dbg, reads, 60000);
        std::vector<AbstractReliableFillingAlgorithm *> algs = {&len, &cov, &bridge, &connect};
        if (diploid)
            algs.emplace_back(&bulge);
        return {std::move(algs)};
    }

#ifdef USE_LIBTORCH
    std::vector<float> loadInferenceResultMultiplicity(const std::experimental::filesystem::path& container_path) {
    torch::jit::script::Module container = torch::jit::load(container_path.string());
    torch::Tensor multiplicity = container.attr("multiplicity").toTensor();
    std::vector<float> mul_vec(multiplicity.data_ptr<float>(), multiplicity.data_ptr<float>() + multiplicity.numel());
    return mul_vec;
}

std::vector<float> loadInferenceResultProbability(const std::experimental::filesystem::path& container_path) {
    torch::jit::script::Module container = torch::jit::load(container_path.string());
    torch::Tensor probability = container.attr("probability").toTensor();
    std::vector<float> prob_vec(probability.data_ptr<float>(), probability.data_ptr<float>() + probability.numel());
    return prob_vec;
}
#endif // USE_LIBTORCH

    class MLReliableFiller : public AbstractReliableFillingAlgorithm {
    private:
        static size_t cnt;
        std::string mode;
        std::experimental::filesystem::path dir;
        double threshold;
    public:
        MLReliableFiller(std::string mode, std::experimental::filesystem::path dir, double threshold) : mode(std::move(mode)), dir(std::move(dir)), threshold(threshold) {}
        std::string name() const override { return "MLReliableFiller"; }
        size_t Fill(dbg::SparseDBG &dbg) override {
#ifndef USE_LIBTORCH
            VERIFY(false);
#endif // USE_LIBTORCH
#ifdef USE_LIBTORCH
            cnt++;
            std::experimental::filesystem::path cur_path = dir / itos(cnt);
            ensure_dir_existance(cur_path);
            logger.info() << "Exporting torch tensors..." << std::endl;
            printPT(cur_path, Component(dbg));
            std::vector<float> inference_results;
            if(mode == "probability") {
                inference_results = loadInferenceResultProbability(cur_path);
            } else if(mode == "multiplicity") {
                inference_results = loadInferenceResultMultiplicity(cur_path);
            } else
                VERIFY_MSG(false, "Unknown reliability checking mode " << mode);
            int i = 0;
            Component comp(dbg);
            for(dbg::Edge &edge : comp.edges()) {
                if(inference_results[i] >= threshold) {
                    edge.is_reliable = true;
                    edge.rc().is_reliable = true;
                }
                i++;
            }
#endif // USE_LIBTORCH
        }
    };

}