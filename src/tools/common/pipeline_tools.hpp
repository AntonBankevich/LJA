#pragma once
#include <malloc.h>

#include <utility>
#include <unordered_set>
#include "cl_parser.hpp"
#include "logging.hpp"
#include "verify.hpp"
#include "sequences/seqio.hpp"

class Stage {//This object represents stage objects and its parameter signature
    //All input and output parameters are uniquely named and their values are of type io::Library that is currently just a vector of file paths.
    //Contract: derivative objects should support move constructor.
    //In general objects derived from stage should strive to be unchangable records. When this is not possible there should be a contract that
    //run method should be called only after all changes are finished. After that the object should actually become unchangable.
    //Implementations of innerRun function should generally be separate from actual stage algorithm. Instead it should focus on interfacing:
    //innerRun should translate input from parameterValues into variables of appropriate types, pass them to the actual function that runs the algorithm, then
    //print output to files in the dir directory and return the paths for these files
protected:
    AlgorithmParameters parameters;// This object defines the set of required and default algorithm parameters
    std::unordered_set<std::string> expected_input;// Names of input libraries that will be used as an input
    std::unordered_set<std::string> expected_output;// Names of output files that will be created by this algorithm

    virtual std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                          const std::experimental::filesystem::path &dir, bool debug,
                          const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) = 0;
public:
    Stage(AlgorithmParameters _params, const std::vector<std::string> &expected_input, const std::vector<std::string> &expected_output) :
            parameters(std::move(_params)), expected_input(expected_input.begin(), expected_input.end()), expected_output(expected_output.begin(), expected_output.end()) {
    }

    const std::unordered_set<std::string> &getExpectedInput() const {return expected_input;}
    const std::unordered_set<std::string> &getExpectedOutput() const {return expected_output;}
    const AlgorithmParameters &getParameters() const {return parameters;}
    AlgorithmParameters getStandaloneParameters() const {
        AlgorithmParameters basic({"output-dir=", "threads=16", "help", "debug"},
                                  {expected_input.begin(), expected_input.end()},
                                  "General parameters:\n"
                                  "\t--output-dir <file_name>  Name of output folder. Results will be stored there.\n"
        "\t--threads <int>  Number of threads to be used by parallel processing.\n");
        basic.AddParameters(parameters, "", "");
        return std::move(basic);
    }

    std::unordered_map<std::string, std::experimental::filesystem::path> run(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                                         bool debug, const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input);
    virtual ~Stage() = default;
};

class ComplexStage;

class SubstageRun {
    //This class envelops a Stage object. It redirects input to the output of other objects and stores output inside itself.
private:
    friend class ComplexStage;
    std::unique_ptr<Stage> stage;
    const ComplexStage *super_stage;
    std::string name;
    std::string prefix;
    bool finished = false;
    std::unordered_map<std::string, std::vector<std::pair<std::string, std::string>>> inputBindings = {};
    std::unordered_map<std::string, std::experimental::filesystem::path> output = {};

    template<class StageType>
    SubstageRun(StageType &&stage, const ComplexStage &super_stage, std::string name, std::string prefix) :
            stage(new StageType(std::forward<StageType>(stage))), super_stage(&super_stage), name(std::move(name)), prefix(std::move(prefix)) {}
public:
    SubstageRun() = default;


    const ComplexStage &getSuperStage() const {return *super_stage;}
    AlgorithmParameterValues readParameterValues(const AlgorithmParameterValues &other);
    void bindInput(const std::string &input_name, const std::string &other_stage_name, const std::string &output_name);
    std::experimental::filesystem::path getResult(const std::string &output_name) const;
    bool isFinished() const {return finished;}
    std::string verifyBinding() const;
    std::string verifyOutput() const;
    std::string verifyOutput(const std::unordered_map<std::string, std::experimental::filesystem::path> &loaded_output) const;
    std::unordered_map<std::string, std::experimental::filesystem::path> readOutput(const std::experimental::filesystem::path &report_candidate) const;
    bool loadAttempt(logging::Logger &logger, const std::experimental::filesystem::path &dir);
    void runSubstage(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, bool debug, const AlgorithmParameterValues &parameterValues);
};

class ComplexStage : public Stage {
private:
    std::vector<std::string> stage_order = {};
    std::unordered_map<std::string, SubstageRun> stages = {};
    std::unordered_map<std::string, io::Library> input_values = {};
public:
    explicit ComplexStage(const std::vector<std::string>& expected_input) :
                        Stage(AlgorithmParameters({"continue", "restart-from=none"}, {},
                      "Pipeline parameters:\n  --continue Continue program from the last save point.\n  --restart-from <stage name>  Restart program from stage (see stage names in log).\n"),
                              expected_input, {}) {
    }

    ComplexStage(ComplexStage &&other)  noexcept : Stage(other) {
        stage_order = std::move(other.stage_order);
        stages = std::move(other.stages);
        input_values = std::move(other.input_values);
        for(auto &it : stages) {
            it.second.super_stage = this;
        }
    }

    template<class StageType>
    SubstageRun &addStage(StageType stage, const std::string &name, std::string prefix = "") {
        if(prefix.empty()) {
            prefix = name + ".";
        }
        parameters.AddParameters(stage.getParameters(), name, prefix);
        stage_order.emplace_back(name);
        stages.insert({name, SubstageRun(std::move(stage), *this, name, prefix)});
        expected_output = stage.getExpectedOutput();
        return stages[name];
    }
    const SubstageRun &getStage(const std::string &name) const;
    io::Library getOutput(const std::string &stage_name, const std::string &parameter_name) const;
    std::string verifyBinding() const;
    void verifyReady() const;
    void verifyResult() const;
    std::unordered_map<std::string, std::experimental::filesystem::path> innerRun(logging::Logger &logger, size_t threads,
                      const std::experimental::filesystem::path &dir, bool debug,
                      const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) override;
};

struct OutputItem {
    OutputItem(string name, string stageOutputName, std::string output_file_name) : name(std::move(name)),
                            outout_id(std::move(stageOutputName)), output_file_name(std::move(output_file_name)) {}

    std::string name;
    std::string outout_id;
    string output_file_name;
};

class LoggedProgram {
private:
    std::string name;
    std::unique_ptr<Stage> stage;
    CLParser clParser;
    std::string start;
    std::string finish;
    std::vector<OutputItem> output;
public:
    template<class StageType>
    LoggedProgram(std::string name, StageType stage, CLParser clParser, std::string start, std::string finish, std::vector<OutputItem> output = {}) :
            name(std::move(name)), stage(new StageType(std::move(stage))), clParser(std::move(clParser)),
            start(std::move(start)), finish(std::move(finish)), output(std::move(output)) {
    }

    template<class StageType>
    LoggedProgram(std::string name, StageType _stage, std::vector<OutputItem> output = {}) :
            name(std::move(name)), stage(new StageType(std::forward<StageType>(_stage))), clParser(stage->getParameters(), {}, {}, 0), output(std::move(output)) {
        AlgorithmParameters params = stage->getParameters();
    }
    int run(const std::vector<std::string> &command_line);
};