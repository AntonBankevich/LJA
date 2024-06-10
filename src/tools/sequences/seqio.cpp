#include "seqio.hpp"

void io::ContigIterator::operator++() {
    reader.inner_read();
    if (reader.eof()) {
        isend = true;
    }
}

io::ContigIterator::ContigIterator(io::IContigReader &_reader, bool _isend) :
    reader(_reader), isend(_isend){
    if (reader.eof()) {
        isend = true;
    }
}

StringContig io::ContigIterator::operator*() {
    return reader.get();
}

bool io::ContigIterator::operator==(const ContigIterator &other) const {
    return isend == other.isend;
}

bool io::ContigIterator::operator!=(const ContigIterator &other) const {
    return isend != other.isend;
}

io::IContigReader::~IContigReader() = default;

const StringContig & io::IContigReader::get() {
    return next;
}

std::vector<StringContig> io::IContigReader::readAll() {
    std::vector<StringContig> res;
    for (StringContig ctg : *this) {
        res.push_back(std::move(ctg));
    }
    return std::move(res);
}

bool io::IContigReader::eof() {
    return next.isNull();
}

io::ContigIterator io::IContigReader::begin() {
    return {*this, false};
}

io::ContigIterator io::IContigReader::end() {
    return {*this, true};
}

StringContig io::IContigReader::read() {
    StringContig string_contig = next;
    inner_read();
    return std::move(string_contig);
}

io::IContigFromFileReader::IContigFromFileReader(const std::experimental::filesystem::path &_file_name) {
    if (endsWith(_file_name.string(), ".gz")) {
        stream = new gzstream::igzstream(_file_name.c_str());
    } else {
        stream = new std::ifstream(_file_name);
    }
}

io::IContigFromFileReader::~IContigFromFileReader() {
    if (stream != nullptr) {
        delete stream;
        stream = nullptr;
    }
}

void io::ISeqReader::nextFile() {
    cur_start = 0;
    cur_end = 0;
    next = {};
    if (subreader != nullptr){
        delete subreader;
        subreader = nullptr;
    }
    if (file_it == lib.end()) {
        return;
    }
    const std::experimental::filesystem::path file_name = *file_it;
    VERIFY_MSG(std::experimental::filesystem::is_regular_file(file_name), "not a regular file: " << file_name);
    initReader(file_name);
    ++file_it;
}


io::ISeqReader::ISeqReader(Library _lib, size_t _min_read_size, size_t _overlap):
    lib(std::move(_lib)), file_it(lib.begin()),
    max_subread_size(2 * _min_read_size - _overlap),
    min_read_size(_min_read_size), overlap(_overlap) {
    VERIFY(min_read_size >= overlap * 2);
}

void io::ISeqReader::inner_read() {
    VERIFY(subreader!=nullptr)
    if (! subreader->eof() && subreader->get().size() == cur_end) {
        subreader->inner_read();
        cur_start = 0;
        cur_end = 0;
    }
    while (subreader->eof()) {
        nextFile();
        if (subreader == nullptr) return;
    }
    size_t cur_pos = (cur_end != 0) ? cur_end - overlap : 0;
    if (subreader->get().size() > max_subread_size + cur_pos) {
        cur_start = cur_pos;
        cur_end = cur_start + min_read_size;
    } else {
        cur_start = cur_pos;
        cur_end = subreader->get().size();
    }
    string id_ext = (cur_start == 0 && cur_end == subreader->get().size()) ? "" : "_" + std::to_string(cur_start);
    next = StringContig(subreader->get().seq.substr(cur_start, cur_end - cur_start),
                        subreader->get().id + id_ext);
}


void io::SeqReader::initReader(const std::experimental::filesystem::path &file_name) {
    if (endsWith(file_name, std::vector<std::string>{".gfa", ".gfa.gz"})) {
        subreader = new GFAReader(file_name);
    } else if (endsWith(file_name, std::vector<std::string>{".fastq", ".fastq.gz", ".fq", "fq.gz"})){
        subreader = new FASTQReader(file_name);
    } else if (endsWith(file_name, std::vector<std::string>{".fasta", ".fasta.gz", ".fa", ".fa.gz"})) {
        subreader = new FASTAReader(file_name);
    } else {
        VERIFY_MSG(false, "unkwn file ext: " << file_name);
    }
}

io::SeqReader::SeqReader(Library _lib, size_t _min_read_size, size_t _overlap):
io::ISeqReader::ISeqReader(_lib, _min_read_size, _overlap) {
    nextFile();
    SeqReader::inner_read();
}

io::SeqReader::SeqReader(const std::experimental::filesystem::path &file_name, size_t _min_read_size, size_t _overlap):
    io::SeqReader::SeqReader(Library({file_name}), _min_read_size, _overlap){}


io::FASTQReader::FASTQReader(const std::experimental::filesystem::path &_file_name):
    IContigFromFileReader(_file_name){
    FASTQReader::inner_read();
}

void io::FASTQReader::inner_read() {
    std::string id, seq;
    std::getline(*stream, id);
    std::getline(*stream, seq);
    if (!id.empty() and !seq.empty()) {
        // verify here
        std::stringstream ss;
        ss << seq;
        while(stream->peek() != EOF && stream->peek() != '+') {
            std::getline(*stream, seq);
            if (seq.empty()) {
                next = {};
                break;
            }
            ss << seq;
        }
        next = {ss.str(), std::move(trim(id.substr(1, id.size() - 1)))};
        std::getline(*stream, id);
        size_t qlen= 0;
        while(!stream->eof() && qlen < next.size()) {
            std::getline(*stream, seq);
            trim(seq);
            qlen += seq.size();
            if (seq.empty())
                break;
        }
        return;
    }
    next = {};
}


io::FASTAReader::FASTAReader(const std::experimental::filesystem::path &_file_name):
    IContigFromFileReader(_file_name) {
    FASTAReader::inner_read();
}

void io::FASTAReader::inner_read() {
    std::string id, seq;
    std::getline(*stream, id);
    std::getline(*stream, seq);
    //trim(seq);
    if (!id.empty() and !seq.empty()) {
        std::stringstream ss;
        ss << seq;
        while(stream->peek() != EOF && stream->peek() != '>') {
            std::getline(*stream, seq);
            VERIFY(seq.empty() || seq[seq.size()-1] != '\n');
            //trim(seq);
            if (seq.empty()) {
                next = {};
                break;
            }
            ss << seq;
        }
        next = {ss.str(), std::move(trim(id.substr(1, id.size() - 1)))};
        return;
    }
    next = {};
}

io::GFAReader::GFAReader(const std::experimental::filesystem::path &_file_name):
    IContigFromFileReader(_file_name){
    GFAReader::inner_read();
}

void io::GFAReader::inner_read() {
    while(stream->peek() != EOF) {
        std::string line;
        std::getline(*stream, line);
        if(line.empty()) {
            break;
        }
        trim(line);
        if (!line.empty() && line[0] != 'L') {
            if(line[0] == 'S') {
                size_t pos = line.find('\t', 2);
                size_t pos2 = line.find('\t', pos + 1);
                if(pos2 == size_t(-1))
                    pos2 = line.size();
                next = {line.substr(pos + 1, pos2 - pos - 1), line.substr(2, pos - 2)};
                return;
            }
        }
    }
    next = {};
}

