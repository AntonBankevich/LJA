#include "pipeline_tools.hpp"
#include <sequences/seqio.hpp>
#include "logging.hpp"
#include "cl_parser.hpp"
#include "verify.hpp"
#include <string>
#include <unordered_set>

void
SubstageRun::runSubstage(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir, bool debug, const AlgorithmParameterValues &parameterValues) {
    logger.stage() << "Starting stage " << name << std::endl;
    timespec start{};
    clock_gettime(CLOCK_MONOTONIC, &start);
    std::string message;
    message = verifyBinding();
    VERIFY_MSG(message.empty(), message);
    std::unordered_map<std::string, io::Library> input;
    for(const auto &it : inputBindings) {
        for(const std::pair<std::string, std::string> &link : it.second) {
            io::Library val = super_stage->getOutput(link.first, link.second);
            input[it.first].insert(input[it.first].end(), val.begin(), val.end());
        }
    }
    output = stage->run(logger, threads, dir, debug, parameterValues, input);
    message = verifyOutput();
    VERIFY_MSG(message.empty(), message);
    timespec finish{};
    finished = true;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    auto worktime = size_t(double(finish.tv_sec - start.tv_sec) + double(finish.tv_nsec - start.tv_nsec) / 1000000000.0);
    std::ofstream os;
    os.open(dir / "report.txt");
    for(auto &it : output) {
        os << it.first << " " << it.second.string() <<"\n";
    }
    os << "Finished stage " << name << " in " << (itos(worktime / 60 % 60, 2)) << " minutes\n";
    os.close();
    logger.info() << "Finished stage " << name << std::endl;
}

void SubstageRun::bindInput(const string &input_name, const string &other_stage_name, const string &output_name) {
    VERIFY(!finished);
    inputBindings[input_name].emplace_back(other_stage_name, output_name);
}

AlgorithmParameterValues SubstageRun::readParameterValues(const AlgorithmParameterValues &other) {return stage->getParameters().fillValuesFrom(prefix, other);}

std::experimental::filesystem::path SubstageRun::getResult(const string &output_name) const {
    VERIFY(finished);
    VERIFY_MSG(output.find(output_name) != output.end(), "Request for unknown output named " + output_name);
    return output.find(output_name)->second;
}

std::string SubstageRun::verifyBinding() const {
    std::stringstream result;
    std::unordered_set<std::string> names;
    for(const std::string &input_name : stage->getExpectedInput()) {
        if(inputBindings.find(input_name) == inputBindings.end())
            result <<  "Stage " + name + " error: input parameter " + input_name + " was not properly bound\n";
        names.emplace(input_name);
    }
    for(const auto &p : inputBindings) {
        if(names.find(p.first) == names.end()) {
            result <<  "Stage " + name + " error: extra input parameter " + p.first + "\n";
        }
    }
    return result.str();
}

std::string SubstageRun::verifyOutput() const {
    return verifyOutput(output);
}

std::string SubstageRun::verifyOutput(
        const std::unordered_map<std::string, std::experimental::filesystem::path> &loaded_output) const {
    std::stringstream result;
    for(const std::string &output_name : stage->getExpectedOutput()) {
        if(loaded_output.find(output_name) == output.end()) {
            result << "Stage " << name << " error: output value " << output_name << " was not reported\n";
        } else if(!std::experimental::filesystem::is_regular_file(loaded_output.find(output_name)->second)) {
            result << "Stage " << name << " error: output file " << loaded_output.find(output_name)->second << " containing output named " <<
                            output_name << " does not exist or is corrupted\n";
        }
    }
    for(const auto &output_item : loaded_output) {
        if(stage->getExpectedOutput().find(output_item.first) == stage->getExpectedOutput().end()) {
            result << "Stage " << name << " error: reported unexpected output value " << output_item.first << " that is stored in file " << output_item.second << "\n";
        }
    }
    return result.str();
}

std::unordered_map<std::string, std::experimental::filesystem::path>
SubstageRun::readOutput(const std::experimental::filesystem::path &report_candidate) const {
    std::unordered_map<std::string, std::experimental::filesystem::path> loaded_output;
    std::string tmp1, tmp2;
    std::ifstream is;
    is.open(report_candidate);
    for(size_t i = 0; i < stage->getExpectedOutput().size(); i++) {
        is >> tmp1 >> tmp2;
        loaded_output[tmp1] = tmp2;
    }
    is.close();
    return std::move(loaded_output);
}

bool SubstageRun::loadAttempt(logging::Logger &logger, const std::experimental::filesystem::path &dir) {
    VERIFY_MSG(!finished, "Attampted to load a finished stage");
    logger.info() << "Attempting to load results of stage " << name << std::endl;
    std::experimental::filesystem::path report_candidate = dir / "report.txt";
    if(!std::experimental::filesystem::is_regular_file(report_candidate)) {
        logger.info() << "No stage report found for stage " << name << std::endl;
        return false;
    }
    std::unordered_map<std::string, std::experimental::filesystem::path> loaded_output = readOutput(report_candidate);
    std::string message = verifyOutput(loaded_output);
    if(!message.empty()) {
        logger.info() << "Corrupted output files for stage " << name << std::endl;
        logger << message;
        return false;
    }
    output = loaded_output;
    logger.info() << "Successfully loaded results of stage " << name << std::endl;
    finished = true;
    return true;
}

std::unordered_map<std::string, std::experimental::filesystem::path> Stage::run(logging::Logger &logger, size_t threads, const std::experimental::filesystem::path &dir,
                                                     bool debug, const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) {
    std::string message = parameters.checkMissingValues(parameterValues);
    VERIFY_MSG(message.empty(), message);
    return innerRun(logger, threads, dir, debug, parameterValues, input);
}

const SubstageRun &ComplexStage::getStage(const string &name) const {
    VERIFY(stages.find(name) != stages.end());
    return stages.find(name)->second;
}

io::Library ComplexStage::getOutput(const std::string &stage_name, const std::string &parameter_name) const {
    if(stage_name.empty()) {
        VERIFY_MSG(input_values.find(parameter_name) != input_values.end(), "Unexpected superstage input parameter name: " + parameter_name);
        return input_values.find(parameter_name)->second;
    } else {
        VERIFY_MSG(stages.find(stage_name) != stages.end(), "Unexpected stage name: " + stage_name);
        return {stages.find(stage_name)->second.getResult(parameter_name)};
    }
}

void ComplexStage::verifyReady() const {
    std::string message;
    for(auto &it : stages) {
        message = it.second.verifyBinding();
        VERIFY_MSG(message.empty(), message);
    }
}

void ComplexStage::verifyResult() const {
    std::string message;
    for(auto &it : stages) {
        message = it.second.verifyOutput();
        VERIFY_MSG(message.empty(), message);
    }
}

std::unordered_map<std::string, std::experimental::filesystem::path>
ComplexStage::innerRun(logging::Logger &logger, size_t threads,
                       const std::experimental::filesystem::path &dir, bool debug,
                       const AlgorithmParameterValues &parameterValues, const std::unordered_map<std::string, io::Library> &input) {
    VERIFY_ERROR_MSG(verifyBinding());
    for(const std::string &iname : Stage::expected_input) {
        input_values[iname] = {};
    }
    for(const auto& p : input) {
        VERIFY_MSG(input_values.find(p.first) != input_values.end(), "Unexpected input name: " + p.first);
        input_values[p.first] = p.second;
    }
    std::string restart_from = parameterValues.getValue("restart-from");
    bool continueRun = parameterValues.getCheck("continue") || (restart_from != "none");
    std::string substage = "";
    size_t pos = restart_from.find('.');
    if(pos != size_t(-1)) {
        substage = substage.substr(pos + 1, restart_from.size() - pos - 1);
        restart_from = restart_from.substr(0, pos);
    }
    if(restart_from != "none" && stages.find(restart_from) == stages.end()) {
        std::cerr << restart_from + " is not a name of any of the pipeline stages" << std::endl;
        std::cerr << "Here is the list of all stage names in the pipeline:" << std::endl;
        for(const std::string& sname: stage_order) {
            std::cerr << sname << std::endl;
        }
    }
    VERIFY(restart_from == "none" || stages.find(restart_from) != stages.end());
    verifyReady();
    size_t cnt = -1;
    for(const std::string& sname: stage_order) {
        cnt++;
        std::experimental::filesystem::path sdir = dir / (itos(cnt, 2) + "_" + sname);
        if(continueRun && restart_from != sname) {
            bool success = stages.find(sname)->second.loadAttempt(logger, sdir);
            if (!success) {
                VERIFY_ERROR_MSG("Critical error: could not load stage results for a stage (" + sname + ") before reaching target stage " + restart_from);
            } else
                continue;
        }
        AlgorithmParameterValues sparameters = stages.find(sname)->second.readParameterValues(parameterValues);
        if(continueRun) {
            logger.info() << "Running pipeline starting from stage " << sname << std::endl;
            continueRun = false;
            if(!substage.empty()) {
                sparameters.addValue("restart-from", substage);
            }
        }
        stages.find(sname)->second.runSubstage(logger, threads, sdir, debug, sparameters);
    }
    verifyResult();
    malloc_trim(0);
    if(!stages.empty()) {
        return stages.find(stage_order.back())->second.output;
    }
    return {};
}

std::string ComplexStage::verifyBinding() const {
    std::stringstream ss;
    for(const auto &p : stages) {
        ss << p.second.verifyBinding();
    }
    return ss.str();
}

int LoggedProgram::run(const std::vector<std::string> &command_line) {
    AlgorithmParameterValues params = clParser.parseCL(command_line);
    if (params.hasCheck("help") && params.getCheck("help")) {
        std::cout << params.helpMessage() << std::endl;
        return 0;
    }
    std::string message = clParser.getParameters().checkMissingValues(params);
    VERIFY_MSG(message.empty(), message);
    VERIFY(params.hasValue("output-dir"));
    bool debug = false;
    if (params.hasCheck("debug") && params.getCheck("debug"))
        debug = true;
    size_t threads = 1;
    if (params.hasValue("threads"))
        threads = std::stoull(params.getValue("threads"));
    const std::experimental::filesystem::path dir(params.getValue("output-dir"));
    ensure_dir_existance(dir);
    logging::LoggerStorage ls(dir, name);
    logging::Logger logger;
    logger.addLogFile(ls.newLoggerFile(), debug ? logging::debug : logging::trace);
    logger.info() << start << std::endl;
    logger.trace() << "Command line:" << std::endl;
    for (const std::string &s : command_line) {
        logger << s << " ";
    }
    logger << std::endl;
    logging::logGit(logger, dir / "version.txt");
    std::unordered_map<std::string, io::Library> input;
    for(const std::string &iname : stage->getExpectedInput()) {
        input[iname] = oneline::initialize<std::experimental::filesystem::path>(params.getListValue(iname));
    }
    std::unordered_map<std::string, std::experimental::filesystem::path> stage_output = stage->run(logger, threads, dir, debug, params, input);
    VERIFY(stage_output.size() == stage->getExpectedOutput().size());
    if(!output.empty()) {
        logger.info() << "Final results can be found in the following file(s): " << std::endl;
        for(auto &rec : output) {
            std::experimental::filesystem::path path = stage_output.at(rec.outout_id);
            std::experimental::filesystem::remove(dir / rec.output_file_name);
            std::experimental::filesystem::copy(path, dir / rec.output_file_name);
        }
        for(auto &rec : output) {
            logger.info() << rec.name << ": " << dir / rec.output_file_name << std::endl;
        }
    } else if(!stage->getExpectedOutput().empty()) {
        logger.info() << "Final results can be found in the following file(s): " << std::endl;
        for(auto &it: stage_output) {
            logger.info() << it.first << ": " << it.second << std::endl;
        }
    }
    logger.info() << finish << std::endl;
    return 0;
}


