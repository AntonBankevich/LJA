//
// Created by anton on 08.07.2020.
//

#pragma once
#include "oneline_utils.hpp"
#include "string_utils.hpp"
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <sstream>
class AlgorithmParameters;
class AlgorithmParameterValues;
class AlgorithmParameters {
protected:
    std::string help_message;
    std::map<std::string, std::string> values = {};
    std::map<std::string, bool> checks = {};
    std::string delim = "^";
public:
    AlgorithmParameters() = default;
    AlgorithmParameters(const std::vector<std::string>& one_value_parameters, const std::vector<std::string>& lib_params, std::string help_message);
    AlgorithmParameters(const std::vector<std::pair<std::string, AlgorithmParameters>> &to_combine, const std::string &help_message);

    void AddParameters(const AlgorithmParameters &other, const std::string &name, const std::string &prefix);

    const std::string &helpMessage() const {return help_message;}
    virtual ~AlgorithmParameters() = default;
    bool hasCheck(const std::string &s) const {return checks.find(s) != checks.end();}
    bool hasValue(const std::string &s) const {return values.find(s) != values.end();}
    bool isList(const std::string &s) const {return values.find(s) != values.end() && !values.find(s)->second.empty() && values.find(s)->second[0] == delim[0];}
    AlgorithmParameterValues fillValuesFrom(const std::string &prefix, const AlgorithmParameterValues & parameterValues) const;
    std::string checkMissingValues(const AlgorithmParameterValues &other) const;
    std::string str() const {
        std::stringstream ss;
        ss << "Values:";
        for(auto &it : values) {
            ss << " " << it.first;
        }
        ss << " Checks:";
        for(auto &it : checks) {
            ss << " " << it.first;
        }
        return ss.str();
    }
};

class AlgorithmParameterValues : public AlgorithmParameters {
public:
    explicit AlgorithmParameterValues(const AlgorithmParameters &default_values) : AlgorithmParameters(default_values) {
    }

    void addValue(const std::string &name, const std::string &val, bool strict = true);
    void addCheck(const std::string &name, bool strict = true);
    const std::string &getValue(const std::string &s) const;
    bool getCheck(const std::string &s) const;
    std::vector<std::string> getListValue(const std::string &s) const {return split(getValue(s), delim);}
    std::string checkMissingValues() const;
    std::string str() const {
        std::stringstream ss;
        ss << "Values:\n";
        for(auto &it : values) {
            ss << it.first << " " << it.second << "\n";
        }
        ss << "Checks:\n";
        for(auto &it : checks) {
            ss << it.first << " " << it.second << "\n";
        }
        return ss.str();
    }
};

class CLParser {
private:
    AlgorithmParameters parameters;
    size_t max_start_size = 0;
    std::map<char, std::string> short_param_map;
    std::unordered_map<std::string, std::vector<std::string>> long_param_map;
    std::vector<std::string> start;
    std::string command_line;
public:
    CLParser(const AlgorithmParameters& parameters, const std::vector<std::string>& short_params, const std::vector<std::string>& long_params = {},
             size_t max_start_size = 1);
    CLParser(CLParser &&other) = default;
    CLParser &operator=(CLParser &&other) = default;

//    TODO: check what happens with quotes
    AlgorithmParameterValues parseCL(const std::vector<std::string>& args, bool strict = true);
    AlgorithmParameterValues parseCL(int argc, char **argv);
    const AlgorithmParameters &getParameters() const {return parameters;}
    const std::vector<std::string> &getStart() const {return start;}
    const std::string &getCL() const {return command_line;}
};