//
// Created by anton on 08.07.2020.
//
#include <iostream>
#include <common/verify.hpp>
#include <utility>
#include "cl_parser.hpp"

const std::string CLParser::emptystring;

void CLParser::parseCL(const std::vector <std::string> &args) {
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
            if(!values[name].empty() && values[name][0] == delim[0]) {
                values[name] += delim + s;
            } else {
                values[name] = s;
            }
            name = "";
        } else if(s[0] == '-') {
            isstart = false;
            if (s[1] == '-') {
                name = s.substr(2, s.size() - 2);
            } else if (s.size() == 2 && short_to_long.count(s[1]) == 1) {
                name = short_to_long[s[1]];
            } else {
                name = "";
            }
            if (checks.count(name)) {
                checks[name] = true;
                name = "";
            } else if (values.count(name) == 0) {
                errors.push_back("Unknown option " + s);
                name = "";
            }
        } else if (isstart) {
            start.push_back(s);
        } else {
            extra.push_back(s);
        }
    }
}

void CLParser::parseCL(int argc, char **argv) {
    parseCL(oneline::initialize<std::string, char*>(argv, argv + argc));
}

const std::string &CLParser::getValue(const std::string &s) const {
    auto it = values.find(s);
    if (it == values.end()) {
        std::cerr << "Missing parameter " << s << std::endl;
        exit(1);
    }
    VERIFY(it != values.end());
    if(it == values.end()) {
        return emptystring;
    } else {
        return it->second;
    }
}

bool CLParser::getCheck(const std::string &s) const {
    VERIFY_MSG(checks.find(s) != checks.end(), s);
    return checks.find(s)->second;
}

CLParser::CLParser(std::vector<std::string> _long_params, std::vector<std::string> _list_params,
                   std::vector<std::string> _short_params, std::string message_) :
        long_params(std::move(_long_params)), list_params(std::move(_list_params)),
        short_params(std::move(_short_params)), _message(std::move(message_)) {
    for(const std::string& s : long_params) {
        size_t pos = s.find('=');
        if (pos != -1) {
            values[s.substr(0, pos)] = s.substr(pos + 1, s.size() - pos - 1);
        } else {
            checks[s] = false;
        }
    }
    for(const std::string& s : list_params) {
        size_t pos = s.find('=');
        if (pos != size_t(-1)) {
            values[s.substr(0, pos)] = s.substr(pos + 1, s.size() - pos - 1);
        } else {
            values[s] = delim;
        }
    }
    for(const std::string& s : short_params) {
        short_to_long[s[0]] = s.substr(2, s.size() - 2);
    }
}
