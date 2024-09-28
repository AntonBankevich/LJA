#pragma once

#include "bulge_path_marker.hpp"
#include "reliable_filler_interface.hpp"
#include "assembly_graph/paths.hpp"
#include "dbg/sparse_dbg.hpp"
#ifdef USE_LIBTORCH
#include "torch/torch.h"
#include "torch/script.h"
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
        std::vector<AbstractReliableFillingAlgorithm*> algs;
        algs.push_back(new CoverageReliableFiller(reliable_threshold));
        algs.push_back(new LengthReliableFiller(20000, 3, 1));
        algs.push_back(new BridgeReliableFiller(40000));
        algs.push_back(new ConnectionReliableFiller(reliable_threshold));
        if (diploid)
        {
            algs.push_back(new BulgePathMarker(dbg, reads, 60000));
        }
        return algs;
    }

#ifdef USE_LIBTORCH
inline std::vector<float> loadInferenceResultMultiplicity(const std::experimental::filesystem::path& container_path) {
    torch::jit::script::Module container = torch::jit::load(container_path.string());
    torch::Tensor multiplicity = container.attr("multiplicity").toTensor();
    std::vector<float> mul_vec(multiplicity.data_ptr<float>(), multiplicity.data_ptr<float>() + multiplicity.numel());
    return mul_vec;
}

inline std::vector<float> loadInferenceResultProbability(const std::experimental::filesystem::path& container_path) {
    torch::jit::script::Module container = torch::jit::load(container_path.string());
    torch::Tensor probability = container.attr("probability").toTensor();
    std::vector<float> prob_vec(probability.data_ptr<float>(), probability.data_ptr<float>() + probability.numel());
    return prob_vec;
}

inline void inference2(dbg::SparseDBG &dbg, std::experimental::filesystem::path &cur_path, std::string config_name) {
    /*
    std::string docker_cmd = "docker run -v /tmp/models:/tmp/models -v /data/cconfigs:/configs -v ";
    // and local dir for files
    docker_cmd += cur_path.string();
    docker_cmd += ":";
    docker_cmd += cur_path.string() + " dbg:latest ";
    */
    // run inference using lja_regression_cfg and set data dir to the current dir
    std::string docker_cmd;
    docker_cmd += "PYTHONPATH=/home/marijo/ws/genomic/fictional-pancake/src; ";
    docker_cmd += "PROJECT_ROOT=/home/marijo/ws/genomic/fictional-pancake; ";
    docker_cmd += "mamba run -n fenv python3";
    docker_cmd += " /home/marijo/ws/genomic/fictional-pancake/src/inference.py --config-name=";
    docker_cmd += config_name + " --config-path /data/cconfigs";
    docker_cmd += " paths.data_dir=";
    docker_cmd += cur_path.string();
    //logger.info() << "Calling docker with cmd " << docker_cmd << std::endl;
    auto result = std::system(docker_cmd.c_str());
    //logger.info() << "Docker exited with status " << result << std::endl;

}

// Function to trim leading and trailing whitespaces, backslashes, and quotes
inline std::string trimAndRemoveQuotes(const std::string& str) {
    std::string trimmedStr = str;
    // Remove leading whitespace, backslashes, and quotes
    trimmedStr.erase(trimmedStr.begin(), std::find_if(trimmedStr.begin(), trimmedStr.end(), [](unsigned char ch) {
        return !std::isspace(ch) && ch != '\\' && ch != '"' && ch != '\'';
    }));
    // Remove trailing whitespace, backslashes, and quotes
    trimmedStr.erase(std::find_if(trimmedStr.rbegin(), trimmedStr.rend(), [](unsigned char ch) {
        return !std::isspace(ch) && ch != '\\' && ch != '"' && ch != '\'';
    }).base(), trimmedStr.end());

    return trimmedStr;
}

inline std::string readDockerCmd(std::string config_name) {
    std::string filename;
    if(config_name == "lja_regression_cfg") {
        filename = "/data/docker/regression_cmd.sh";
    } else if(config_name == "lja_classification_cfg") {
        filename = "/data/docker/classification_cmd.sh";
    } else {
        std::cerr << "Failed with wrong config_name: " << config_name << std::endl;
        return "";
    }
    std::ifstream inputFile(filename.c_str()); // Replace with your filename
    if (!inputFile.is_open()) {
        std::cerr << "Failed to open the file: " << filename << std::endl;
        return "";
    }
    std::string line;
    std::stringstream result;
    bool firstLine = true;

    while (std::getline(inputFile, line)) {
        std::string trimmedLine = trimAndRemoveQuotes(line);

        if (!trimmedLine.empty()) {
            if (!firstLine) {
                result << " ";
            }
            result << trimmedLine;
            firstLine = false;
        }
    }

    inputFile.close();

    std::string finalResult = result.str();
    std::cout << "Concatenated string: " << finalResult << std::endl;

    return finalResult;
}

inline void inference(dbg::SparseDBG &dbg, std::experimental::filesystem::path &cur_path, std::string config_name) {
    std::string docker_cmd = readDockerCmd(config_name);
    docker_cmd += " paths.data_dir=";
    docker_cmd += cur_path.string();
    std::cerr << "Calling Docker with cmd " << docker_cmd << std::endl;
    auto result = std::system(docker_cmd.c_str());
    std::cerr << "Docker exit with status " << result << std::endl;
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
            //printPT(cur_path, Component(dbg));
            printPTDBG(cur_path, dbg);

            std::vector<float> inference_results;

            if(mode == "probability") {
                inference(dbg, cur_path, "lja_classification_cfg");
                inference_results = loadInferenceResultProbability(cur_path / "container.pt");
            } else if(mode == "multiplicity") {
                inference(dbg, cur_path, "lja_regression_cfg");
                inference_results = loadInferenceResultMultiplicity(cur_path / "container.pt");
            } else {
                VERIFY_MSG(false, "Unknown reliability checking mode " << mode);
            }
            size_t i = 0;
            size_t i_un = 0;
            std::cerr << "Using mode " << mode << " with ML threshold " << threshold << std::endl;
            size_t index = 0;
            std::ofstream dbg_edges(cur_path / "edge_attrs_post.txt");
            for(dbg::Edge &edge : dbg.edges()) {
                double coverage = edge.getCoverage();
                int size = edge.truncSize();
                dbg_edges << coverage << " " << size << std::endl;
                if(inference_results[index] >= threshold) {
                    edge.is_reliable = true;
                    edge.rc().is_reliable = true;
                    i++;
                } else {
                    edge.is_reliable = false;
                    edge.rc().is_reliable = false;
                    i_un++;
                }
                index++;
            }
            dbg_edges.close();
            std::cerr << "Found total " << i << " reliable edges and " << i_un << " unreliable edges" << std::endl;
            return i;
#endif // USE_LIBTORCH
        }
    };

}
