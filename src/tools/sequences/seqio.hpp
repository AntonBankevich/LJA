#pragma once

#include "common/string_utils.hpp"
#include "stream.hpp"
#include "contigs.hpp"
#include <experimental/filesystem>
#include <iterator>
#include <string>
#include <utility>
#include <vector>
#include <functional>

namespace io{
    typedef std::vector<std::experimental::filesystem::path> Library;

    class IContigReader;


    class ContigIterator{
    private:
        IContigReader &reader;
        bool isend;
    public:
        typedef StringContig value_type;
        typedef const StringContig& reference;
        ContigIterator(IContigReader &_reader, bool _isend);
        void operator++();
        StringContig operator*();
        bool operator==(const ContigIterator &other) const;
        bool operator!=(const ContigIterator &other) const;
    };


    class IContigReader {
    protected:
        StringContig next{};
    public:
        virtual ~IContigReader();
        const StringContig& get();
        bool eof();
        virtual void inner_read() = 0;
        ContigIterator begin();
        ContigIterator end();
        virtual void reset() = 0;
        StringContig read();
    };

    class IContigFromFileReader: public IContigReader {
    protected:
        std::istream* stream = nullptr;
    public:
        explicit IContigFromFileReader(const std::experimental::filesystem::path& _file_name);
        ~IContigFromFileReader() override;
        void reset() override;
    };

    class SeqReader final: public IContigReader{
    private:
        const Library lib;
        Library::const_iterator file_it;
        size_t max_subread_size;
        size_t min_read_size;
        size_t overlap;
        size_t cur_start = 0;
        size_t cur_end = 0;
        IContigReader* subreader = nullptr;

        void nextFile();

    public:
        explicit SeqReader(Library _lib, size_t _min_read_size = size_t(-1) / 2, size_t _overlap = size_t(-1) / 8);
        explicit SeqReader(const std::experimental::filesystem::path & file_name,
                           size_t _min_read_size = size_t(-1) / 2, size_t _overlap = size_t(-1) / 8);
        void inner_read() override;
        void reset() override;
    };


    class FASTQReader final: public IContigFromFileReader{
    public:
        explicit FASTQReader(const std::experimental::filesystem::path& _file_name);
        void inner_read() override;
    };


    class FASTAReader final : public IContigFromFileReader{
    public:
        explicit FASTAReader(const std::experimental::filesystem::path &_file_name);
        void inner_read() override;
    };


    class GFAReader final: public IContigFromFileReader{
    public:
        explicit GFAReader(const std::experimental::filesystem::path& _file_name);
        void inner_read() override;
    };
}

inline io::Library operator+(const io::Library &lib1, const io::Library &lib2) {
    io::Library res = lib1;
    res.insert(res.end(), lib2.begin(), lib2.end());
    return std::move(res);
}
