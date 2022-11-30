//
// Created by anton on 08.07.2020.
//
#include "cl_parser.hpp"
#include <iostream>
#include <common/verify.hpp>
#include <utility>
#include <unordered_set>

AlgorithmParameterValues CLParser::parseCL(const std::vector <std::string> &args, bool strict) {
    AlgorithmParameterValues result(parameters);
    std::stringstream ss;
    for(const std::string & arg : args) {
        ss << arg;
        ss << " ";
    }
    command_line = command_line + ss.str();
    bool isstart = true;
    std::string name;
    for(const std::string &s : args) {
        if (!name.empty()) {
            if(long_param_map.find(name) != long_param_map.end()) {
                for(const std::string &actual_name : long_param_map[name]) {
                    result.addValue(actual_name, s, strict);
                }
            } else {
                result.addValue(name, s, strict);
            }
            name = "";
        } else if(s[0] == '-') {
            isstart = false;
            if (s[1] == '-') {
                name = s.substr(2, s.size() - 2);
            } else {
                VERIFY_MSG(s.size() == 2, "Long command line parameters should start with --");
                VERIFY_MSG(!strict || short_param_map.find(s[1]) != short_param_map.end(), "Unknown command line parameter: " + s);
                if(short_param_map.find(s[1]) == short_param_map.end()) {
                    name = "";
                    continue;
                }
                name = short_param_map[s[1]];
            }
            VERIFY_MSG(!strict || result.hasCheck(name) || result.hasValue(name), "Unknown command line parameter: " + s);
            if (result.hasCheck(name)) {
                if(long_param_map.find(name) != long_param_map.end()) {
                    for(const std::string &actual_name : long_param_map[name]) {
                        result.addCheck(actual_name, strict);
                    }
                } else {
                    result.addCheck(name, strict);
                }
                name = "";
            }
        } else {
            VERIFY_MSG(!strict || isstart, "Incorrect command line: " << name << " is neither a value nor a parameter. Did you forget \"--\" in front of an option name?");
            if(isstart)
                start.push_back(s);
        }
    }
    VERIFY_MSG(!strict || start.size() <= max_start_size, "Incorrect command line: " << start[max_start_size] << " is neither a value nor a parameter. Did you forget \"--\" in fromt of an option name?");
    return std::move(result);
}

AlgorithmParameterValues CLParser::parseCL(int argc, char **argv) {
    return parseCL(oneline::initialize<std::string, char*>(argv, argv + argc));
}

CLParser::CLParser(const AlgorithmParameters& parameters, const std::vector<std::string>& short_params,
                   const std::vector<std::string>& long_params, size_t max_start_size) :
        parameters(parameters), max_start_size(max_start_size) {
    for(const std::string& s : short_params) {
        short_param_map[s[0]] = s.substr(2, s.size() - 2);
    }
    for(const std::string& s : long_params) {
        size_t pos = s.find('=');
        long_param_map[s.substr(0, pos)].emplace_back(s.substr(pos + 1, s.size() - pos - 1));
    }
}

AlgorithmParameters::AlgorithmParameters(const std::vector<std::string>& one_value_parameters, const std::vector<std::string>& lib_params,
                                         std::string help_message) : help_message(std::move(help_message)) {
    for(const std::string& s : one_value_parameters) {
        size_t pos = s.find('=');
        if (pos != -1) {
            values[s.substr(0, pos)] = s.substr(pos + 1, s.size() - pos - 1);
        } else {
            checks[s] = false;
        }
    }
    for(const std::string& s : lib_params) {
        size_t pos = s.find('=');
        if (pos != size_t(-1)) {
            values[s.substr(0, pos)] = s.substr(pos + 1, s.size() - pos - 1);
        } else {
            values[s] = delim;
        }
    }
}

AlgorithmParameters::AlgorithmParameters(const std::vector<std::pair<std::string, AlgorithmParameters>> &to_combine,
                                         const std::string &help_message) {
    std::unordered_set<std::string> all_params;
    std::unordered_set<std::string> all_libs;
    std::unordered_set<std::string> all_checks;
    std::unordered_set<std::string> prefixes;
    std::stringstream message;
    message << help_message;
    for(auto &p : to_combine) {
        VERIFY_MSG(prefixes.find(p.first) == prefixes.end(), "Duplicate stage prefix: " + p.first);
        prefixes.insert(p.first);
        VERIFY_MSG(!p.first.empty(), "Empty prefix for subalgorithm " + p.first);
        message << "\n" << "Parameters of " << p.first << "\n" << p.second.help_message;
        for(auto &pvalue : p.second.values) {
            values[p.first + pvalue.first] = pvalue.second;
            if(pvalue.second.empty() || pvalue.second[0] != delim[0])
                all_params.insert(pvalue.first);
            else
                all_checks.insert(pvalue.first);
        }
        for(auto &pvalue : p.second.checks) {
            checks[p.first + pvalue.first] = pvalue.second;
            all_checks.insert(pvalue.first);
        }
    }
    for(const std::string &pname : all_params) {
        values[pname] = "";
    }
    for(const std::string &pname : all_checks) {
        if(all_params.find(pname) == all_params.end())
            checks[pname] = false;
    }
    for(const std::string &pname : all_libs) {
        if(all_params.find(pname) == all_params.end() && all_checks.find(pname) == all_checks.end())
            values[pname] = delim;
    }
    this->help_message = message.str();
}

AlgorithmParameterValues AlgorithmParameters::fillValuesFrom(const std::string &prefix, const AlgorithmParameterValues &parameterValues) const {
    AlgorithmParameterValues result(*this);
    for(const auto & key : values) {
        VERIFY_MSG(parameterValues.hasValue(prefix + key.first), prefix + " " + key.first);
//        VERIFY_MSG(parameterValues.hasValue(key.first), key.first);
        if (!parameterValues.getValue(prefix + key.first).empty()) {
            result.addValue(key.first, parameterValues.getValue(prefix + key.first));
        } else if(parameterValues.hasValue(key.first)) {
            result.addValue(key.first, parameterValues.getValue(key.first));
        }
    }
    for(const auto & key : checks) {
        VERIFY(parameterValues.hasCheck(prefix + key.first));
//        VERIFY(parameterValues.hasCheck(key.first));
        if(parameterValues.getCheck(prefix + key.first) || (parameterValues.hasCheck(key.first) && parameterValues.getCheck(key.first)))
            result.addCheck(key.first);
    }
    return std::move(result);
}

std::string AlgorithmParameters::checkMissingValues(const AlgorithmParameterValues &other) const {
    std::stringstream result;
    for(const auto & key : values) {
        if(!other.hasValue(key.first)) {
            result <<  "Missing parameter " << key.first << "\n";
        } else if(other.getValue(key.first).empty()) {
            result <<  "Missing parameter value " << key.first << "\n";
        }
    }
    for(const auto & key : checks) {
        if(!other.hasCheck(key.first)) {
            result <<  "Missing check " << key.first << "\n";
        }
    }
    return result.str();
}

AlgorithmParameters AlgorithmParameters::AddParameters(const AlgorithmParameters &other, const std::string &name,
                                                       const std::string &prefix) const {
    AlgorithmParameters res = *this;
    res.help_message = res.help_message + "\nParameters of " + name + "\n" + other.help_message;

    std::unordered_set<std::string> all_params;
    std::unordered_set<std::string> all_libs;
    std::unordered_set<std::string> all_checks;
    std::unordered_set<std::string> prefixes;
    std::stringstream message;
    message << help_message;
    for(auto &pvalue : other.values) {
        std::string pname = prefix + pvalue.first;
        VERIFY_MSG(res.values.find(pname) == res.values.end(), "Duplicate parameter " + pname);
        VERIFY_MSG(res.checks.find(pname) == res.checks.end(), "Duplicate parameter " + pname);
        res.values[pname] = pvalue.second;
    }
    for(auto &pvalue : other.checks) {
        std::string pname = prefix + pvalue.first;
        VERIFY_MSG(res.values.find(pname) == res.values.end(), "Duplicate parameter " + pname);
        VERIFY_MSG(res.checks.find(pname) == res.checks.end(), "Duplicate parameter " + pname);
        res.checks[pname] = pvalue.second;
    }
    return std::move(res);
}

void AlgorithmParameterValues::addValue(const std::string &name, const std::string &val, bool strict) {
    VERIFY(!strict || hasValue(name));
    if(!hasValue(name))
        return;
    if(isList(name)) {
        size_t start = 0;
        size_t end = val.size();
        while(start < end && val[start] == delim[0])
            start++;
        while(end > start && val[end - 1] == delim[0])
            end--;
        if(end > start)
            values[name] += val.substr(start, end - start) + delim;
    } else {
        values[name] = val;
    }
}

void AlgorithmParameterValues::addCheck(const std::string &name, bool strict) {
    VERIFY(!strict || hasCheck(name));
    if(hasCheck(name))
        checks[name] = true;
}

bool AlgorithmParameterValues::getCheck(const std::string &s) const {
    VERIFY_MSG(checks.find(s) != checks.end(), s);
    return checks.find(s)->second;
}

const std::string &AlgorithmParameterValues::getValue(const std::string &s) const {
    auto it = values.find(s);
    if (it == values.end()) {
        std::cerr << "Missing parameter " << s << std::endl;
        exit(1);
    }
    return it->second;
}

std::string AlgorithmParameterValues::checkMissingValues() const {
    return AlgorithmParameters::checkMissingValues(*this);
}