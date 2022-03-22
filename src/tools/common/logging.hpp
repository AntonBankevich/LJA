//
// Created by anton on 7/27/20.
//
#pragma once
#include "output_utils.hpp"
#include "dir_utils.hpp"
#include "string_utils.hpp"
#include "malloc.h"
#include <omp.h>
#include "sys/sysinfo.h"
#include <sys/resource.h>
#include <experimental/filesystem>
#include <string>
#include <sstream>
#include <utility>
#include <vector>
#include <iostream>
#include <ostream>
#include <fstream>
#include <chrono>

namespace logging {

    const std::string endl = "\n";

    class TimeSpace {
    private:
        std::chrono::time_point<std::chrono::system_clock> start;
    public:
        TimeSpace() : start(std::chrono::system_clock::now()) {
        }

        std::string get() const {
            std::chrono::time_point<std::chrono::system_clock> finish = std::chrono::system_clock::now();
            std::chrono::duration<double> seconds = finish - start;
            auto full_seconds = size_t(seconds.count());
            struct sysinfo memInfo;
            sysinfo (&memInfo);
            struct rusage usage;
            getrusage(RUSAGE_SELF, &usage);
            double all = double(usage.ru_maxrss) / 1024;
            std::string t_used = "Mb";
            std::string t_all = "Mb";
            if (all > 700) {
                all = all / 1024;
                t_all = "Gb";
            }
            all = size_t(all * 100) * 0.01;
            std::stringstream ss;
            ss << itos(full_seconds / 60 / 60, 2) << ":" << itos(full_seconds / 60 % 60, 2) << ":"
                    << itos(full_seconds % 60, 2) << " "  << all << t_all << " ";
            return ss.str();
        }
    };

    class LoggerStorage {
    private:
        const std::experimental::filesystem::path dir;
        const std::experimental::filesystem::path logFile;
        const std::experimental::filesystem::path backupDir;
    public:
        explicit LoggerStorage(std::experimental::filesystem::path _dir, const std::string &_programName):
                dir(std::move(_dir)), logFile(dir / (_programName + ".log")), backupDir(dir / "old_logs") {
        }

        std::experimental::filesystem::path backup() const {
            if(!std::experimental::filesystem::is_regular_file(logFile)) {
                return {};
            }
            ensure_dir_existance(backupDir);
            size_t max = 0;
            for (const std::experimental::filesystem::path & file : std::experimental::filesystem::directory_iterator(backupDir)) {
                std::string fname = file.filename().string();
                if (fname.size() < 5 || fname.substr(fname.size() - 4) != ".log")
                    continue;
                try {
                    max = std::max<size_t>(max, std::stoi(fname.substr(0, fname.size() - 4)));
                } catch (const std::invalid_argument& ia) {
                }
            }
            std::experimental::filesystem::path backup = backupDir / (itos(max + 1, 2) + ".log");
            std::experimental::filesystem::copy_file(logFile, backup);
            std::experimental::filesystem::remove(logFile);
            return std::move(backup);
        }

        std::experimental::filesystem::path newLoggerFile() const {
            backup();
            return logFile;
        }
    };

    //info goes to console, trace goes to log file. Debug goes to log file if debug is enabled
    //stage like info but for large stage declaration
    enum LogLevel {stage, info, trace, debug};

    class Logger : public std::streambuf , public std::ostream {
    private:
        struct LogStream {
            std::ofstream * os;
            LogLevel level;
            LogStream(const std::experimental::filesystem::path &fn, LogLevel level) :
                    os(new std::ofstream()), level(level) {
                os->open(fn);
            }

            ~LogStream() {
                os->close();
                delete os;
                os = nullptr;
            }
        };

        std::vector<LogStream> oss;
        TimeSpace time;
        Logger *empty_logger = nullptr;
        LogLevel curlevel;
        bool add_cout;
    public:
        explicit Logger(bool _add_cout = true) :
                    std::ostream(this), curlevel(LogLevel::trace), add_cout(_add_cout) {
        }

        Logger(const Logger &) = delete;

        void addLogFile(const std::experimental::filesystem::path &fn, LogLevel level = LogLevel::trace) {
            oss.emplace_back(fn, level);
        }

    //    template<class T>
    //    DummyLogger &operator<<(const T &val) {
    //        std::cout << time.get() << val;
    //        for(std::ofstream *os : oss) {
    //            *os << time.get() << val;
    //        }
    //        return dummyLogger;
    //    }

        int overflow(int c) override {
            if(curlevel <= LogLevel::info)
                std::cout << char(c);
            for(LogStream &os : oss) {
                if(curlevel <= os.level)
                    *os.os << char(c);
            }
            if(c == '\n') {
                forceFlush();
            }
            return 0;
        }

        void forceFlush() {
            if(curlevel <= LogLevel::info)
                std::cout.flush();
            for(LogStream &os : oss) {
                if(curlevel <= os.level)
                    os.os->flush();
            }
        }

        Logger & stage() {
            curlevel = LogLevel::stage;
            *this << "==============================================================\n";
            *this << time.get() << " NEW STAGE: ";
            return *this;
        }

        Logger & info() {
            curlevel = LogLevel::info;
            *this << time.get() << " INFO: ";
            return *this;
        }

        Logger & trace() {
            curlevel = LogLevel::trace;
            *this << time.get() << " TRACE: ";
            return *this;
        }

        Logger & debug() {
            curlevel = LogLevel::debug;
            *this << time.get() << " DEBUG: ";
            return *this;
        }

        ~Logger() override {
            delete empty_logger;
        }
    };

    class ProgressBar {
    private:
        Logger &logger;
        size_t cur;
        size_t size;
        size_t bar_num;
        size_t mod;
        std::vector<size_t> progress;
        size_t cnt;
    public:
        ProgressBar(Logger &logger, size_t size, size_t threads, size_t bar_num = 20) : logger(logger), cur(1),
                            size(size), bar_num(std::min(bar_num, size)),
                            mod(std::max<size_t>(1, size / bar_num / threads / 5)), progress(threads), cnt(0) {

        }

        void tick() {
            size_t thread = omp_get_thread_num();
            progress[thread]++;
#pragma omp atomic
            cnt++;
            if(progress[thread] % mod == 0) {
#pragma omp critical
                {
                    size_t tmp = cnt;
                    if(cnt * bar_num >= size * cur) {
                        logger << "|";
                    }
                }
            }
        }

        void finish() {
            logger << std::endl;
        }
    };

    inline void logGit(Logger &logger, const std::experimental::filesystem::path &out) {
        int code = system(("bash -c \"git rev-parse HEAD;git diff\" 2> /dev/null > " + out.string()).c_str());
        std::ifstream is;
        is.open(out);
        std::string res;
        is >> res;
        if(!res.empty())
            logger.info() << res << std::endl;
        res = "";
        logger.trace() << "Git diff:" << std::endl;
        while (getline(is,res)) {
            logger << res << '\n';
        }
        is.close();
    }
}
