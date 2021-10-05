//
// Created by anton on 08.07.2020.
//

#pragma once
#include "oneline_utils.hpp"
#include "string_utils.hpp"
#include <vector>
#include <string>
#include <map>
#include <sstream>

class CLParser {
private:
    const std::vector<std::string> long_params;
    const std::vector<std::string> list_params;
    const std::vector<std::string> short_params;
    std::map<std::string, std::string> values;
    std::map<std::string, bool> checks;
    std::map<char, std::string> short_to_long;
    std::vector<std::string> start;
    std::vector<std::string> extra;
    std::vector<std::string> errors;
    const static std::string emptystring;
    std::string command_line;
    std::string _message;
    std::string delim = "^";
public:
    CLParser(std::vector<std::string> _long_params, std::vector<std::string> _list_params,
             std::vector<std::string> _short_params, std::string message_ = "");

//    TODO: check what happens with quotes
//    TODO: make failsafe
    void parseCL(const std::vector<std::string>& args);

    void parseCL(int argc, char **argv);

    std::string check() {
        for(const auto & key : values) {
            if (key.second.empty()) {
                return key.first + " missing";
            }
        }
        if(!errors.empty()) {
            return errors[0];
        }
        if(!extra.empty()) {
            return "Unknown parameter(s): " + join(" ", extra);
        }
        return "";
    }

    const std::string &message() const {
        return _message;
    }

    const std::string & getValue(const std::string &s) const;

    std::vector<std::string> getListValue(const std::string &s) const {
        return split(getValue(s), delim);
    }

    bool getCheck(const std::string &s) const;

    const std::vector<std::string> &getStart() const {
        return start;
    }

    const std::vector<std::string> &getExtra() const {
        return extra;
    }

    const std::vector<std::string> &getErrors() const {
        return errors;
    }

    const std::string &getCL() const {
        return command_line;
    }
};