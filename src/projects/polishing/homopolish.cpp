#include "homopolish.hpp"
#include <sequences/contigs.hpp>
#include <common/cl_parser.hpp>
#include <common/logging.hpp>
#include <sequences/seqio.hpp>
#include <common/omp_utils.hpp>
#include <common/zip_utils.hpp>
#include <unordered_set>
#include <utility>
#include <vector>
#include <iostream>
#include <array>
#include <spoa/spoa.hpp>
#include <ksw2/ksw_wrapper.hpp>

using std::vector;
using std::pair;
using std::array;
using std::sort;
using logging::Logger;


struct AlignmentInfo {
    string read_id;
    string contig_id;
    size_t read_start = 0;
    size_t read_end = 0;
    size_t alignment_start = 0;
    size_t alignment_end = 0;
    bool rc = false;

    size_t length() const {
        return alignment_end - alignment_start;
    }

    string str() const {
        std::stringstream ss;
        ss << read_id << " " << contig_id << " " <<alignment_start << " " << alignment_end << " "<< rc;
        return ss.str();
    }
};


struct dinucleotide {
    size_t start, multiplicity;
    //do we need sequence?
    string seq;
    static const size_t MIN_DINUCLEOTIDE_REPEAT = 10;

    dinucleotide(size_t start, size_t multiplicity, string seq): start(start), multiplicity(multiplicity), seq(std::move(seq)) {
    }

    string str() const {
        std::stringstream ss;
        ss << start << " " << multiplicity << " " <<seq;
        return ss.str();
    }
};

//complex region positions sorted
template <class Contig>
vector<dinucleotide> GetDinucleotideRepeats(const Contig& compressed_contig) {
    size_t len = compressed_contig.size();
    size_t start = 0;
    vector<dinucleotide> res;
    while (start < len) {
        size_t multiplicity = 1;
        while ((start + 2 * multiplicity + 1 < len) &&
               compressed_contig[start] == compressed_contig[start + 2 * multiplicity ] &&
               compressed_contig[start + 1] == compressed_contig[start + 2 * multiplicity + 1]) {
            multiplicity++;
        }

        if (multiplicity >= dinucleotide::MIN_DINUCLEOTIDE_REPEAT)
        {
            std::stringstream ss;
            ss << compressed_contig[start] << compressed_contig[start+1];
            res.emplace_back(start, multiplicity, ss.str());
//Dirty hack to avoid reporting overlapping dinucleotide repeats
            start ++;
        }
        start = start + 2 * multiplicity - 1;
    }
    return res;
}

struct ContigInfo {
    string sequence;
    string name;
    size_t len = 0;
    vector<uint8_t > quantity;
//array size. possibly switch to []
    static const size_t VOTES_STORED = 21;
    //To avoid multiple resizes real length stored in first element of each vector
    vector<array<uint8_t, VOTES_STORED + 1>> amounts;
    vector<uint16_t> sum;
    size_t zero_covered = 0;
//neighbourhoud for complex regions;
    static const size_t COMPLEX_EPS = 5;
    static const size_t MAX_CONSENSUS_COVERAGE = 25;
    static constexpr double MAX_ALLOWED_MSA_LENGTH_VARIATION = 1.1;
    //Complex_regions: map from start to length;

    vector<pair<size_t, size_t>> complex_regions;
    std::map<size_t, vector<string>> complex_strings;
//TODO more efficient data structure

/*    vector <dinucleotide> dinucleotide_coords;
    vector <int> num_dinucleotide;
    vector<vector<size_t>> dinucleotide_read_voting;
*/
    void FillComplex() {
        size_t current_id = 0;
        auto dinucleotide_coords = GetDinucleotideRepeats(sequence);
        size_t total_d = dinucleotide_coords.size();
        while (current_id < total_d) {
            size_t cur_start = get_start(dinucleotide_coords[current_id]);
            size_t cur_finish = get_finish(dinucleotide_coords[current_id]);
            size_t next_id = current_id + 1;
            while ((next_id) < total_d && get_start(dinucleotide_coords[next_id]) <= cur_finish) {
                cur_finish = std::max(cur_finish, get_finish(dinucleotide_coords[next_id]));
                next_id ++;

            }
            if (cur_finish > len)
                cur_finish = len;
            complex_regions.emplace_back(std::make_pair(cur_start, cur_finish - cur_start));
            complex_strings[cur_start] = vector<string>();
            current_id = next_id;
        }

    }

    ContigInfo(const string& sequence, string name) :
                                sequence(sequence), name(std::move(name)) {
        len = sequence.length();
        quantity.resize(len);
        sum.resize(len);
        amounts.resize(len);
        for (size_t i = 0; i < len; i ++){
            sum[i] = 0;
            quantity[i] = 0;
//            amounts[i].resize(VOTES_STORED + 1);
            amounts[i][0] = 0;
        }
        FillComplex();
    }
    ContigInfo() = default;

    size_t get_finish(const dinucleotide& d) {
        return d.start + COMPLEX_EPS + d.multiplicity * 2;
    }

    size_t get_start(const dinucleotide& d) {
        if (d.start < COMPLEX_EPS) return 0;
        else return d.start - COMPLEX_EPS;
    }


//    void AddRead() {
//
//    }

    bool validLength(size_t med_len, size_t some) {
        return (med_len < some *MAX_ALLOWED_MSA_LENGTH_VARIATION && med_len * MAX_ALLOWED_MSA_LENGTH_VARIATION > some);
    }

    string MSAConsensus(vector<string> &s, Logger & logger) {
//Magic consts from spoa default settings
        auto alignment_engine = spoa::AlignmentEngine::Create(
// -8 in default for third parameter(gap) opening, -6 for forth(gap extension)
                spoa::AlignmentType::kNW, 10, -8, -8, -1);  // linear gaps
        spoa::Graph graph{};
        if (s.size() == 0) {
#pragma OMP critical
            logger.trace() << "WARNING: zero strings were provided for consensus counting" << endl;
            return "";
        }
        size_t cov = 0;
        vector<size_t> all_len;
        for (const auto& it: s) {
            all_len.push_back(it.length());
        }
        sort(all_len.begin(), all_len.end());
        size_t med_len = all_len[all_len.size()/2];
        size_t good_cons = 0;
        for (size_t i = 0; i < all_len.size(); i++) {
            if (med_len < all_len[i] *1.1 && med_len * 1.1> all_len[i]){
                good_cons ++;
            }
        }
        bool only_good = true;
        if (good_cons < 5) {
            only_good = false;
            logger.debug() << "Few good consensus seqs: " << good_cons << " of " << all_len.size() << endl;
        }

        for (const auto& it : s) {
            if (only_good && !validLength(med_len, it.length())) {
                continue;
            }
            std::int32_t score = 0;
            auto alignment = alignment_engine->Align(it, graph, &score);
            graph.AddAlignment(alignment, it);
            cov ++;
            if (cov == MAX_CONSENSUS_COVERAGE) {
                break;
            }
        }

        vector<uint32_t > coverages;
        string consensus = graph.GenerateConsensus(&coverages);
        size_t pref_remove = 0;
        int suf_remove = int(coverages.size()) - 1;
        while (pref_remove < coverages.size() && coverages[pref_remove] < cov / 2 )
            pref_remove ++;
        while (suf_remove >= 0 && coverages[suf_remove] < cov / 2 )
            suf_remove --;
        if (pref_remove > suf_remove) {
            return "";
        }
        return consensus.substr(pref_remove, suf_remove - pref_remove + 1);

    }

    string checkMSAConsensus(string s, vector<string> &all) {
        size_t cons_len = s.length();
//consensus in last
        vector<size_t> all_len(all.size() - 1);
        for (size_t i = 0; i < all.size() - 1; i ++) {
            all_len[i] = all[i].length();
        }
        sort(all_len.begin(), all_len.end());
        s.erase(std::unique(s.begin(), s.end()), s.end());
        if (s.length() <= 2 * COMPLEX_EPS + 2) return "TOO SHORT";
        if (cons_len != all_len[all_len.size()/2]) {
            return "CONSENSUS LENGTH DIFFER FROM MEDIAN";
        }
        return "";
//What are other suspicious cases? Since we can glue two dimeric regions, commented case is  actually OK
/*        size_t len = s.length();
        for (size_t i = COMPLEX_EPS + 2; i < len - COMPLEX_EPS; i++) {
            if (s[i] != s[COMPLEX_EPS] && s[i] != s[COMPLEX_EPS + 1])
                return "On position " + to_string(i) + " of " + to_string(len) + " not dimeric nucleo";
        } */
    }

    string GenerateConsensus(Logger & logger){
        std::stringstream ss;
        vector<int> quantities(256);
        std::ofstream debug;
        size_t total_count = 0 ;
        size_t cur_complex_ind = 0;
#pragma omp parallel for default(none) shared(logger)
        for (size_t i = 0; i < complex_regions.size(); i++) {
            size_t start_pos = complex_regions[i].first;
            auto consensus = MSAConsensus(complex_strings[start_pos], logger);
            complex_strings[start_pos].push_back(consensus);
        }
        logger.debug() << " Consenus for contig " << name << " calculated "<< endl;
        string consensus;


        for (size_t i = 0; i < len; ) {
            if (!complex_regions.empty() && complex_regions[cur_complex_ind].first == i ) {
                consensus = complex_strings[i][complex_strings[i].size() - 1];
                ss << consensus;
                auto check = checkMSAConsensus(consensus, complex_strings[i]);
 //            logger.info() << "consensus of " << complex_strings[start_pos].size() << ": " << consensus.length() << endl << "At position " <<start_pos << endl;
                if (!check.empty()){
                    logger.debug() << "Problematic consensus starting on decompressed position " << total_count <<" " << check <<" of " <<complex_strings[i].size() - 1 << " sequences "<< endl;
                    logger.debug() << "Position " << complex_regions[cur_complex_ind].first << " len " << complex_regions[cur_complex_ind].second << endl;
                    std::stringstream debug_l;
                    debug_l << "lengths: ";
                    for (size_t j = 0; j < complex_strings[i].size() - 1; j++) {
                        debug_l << complex_strings[i][j].length() << " ";
                    }
                    debug_l <<" : " << consensus.length() << endl;
                    logger.debug() << debug_l.str();
                    for (size_t j = 0; j < complex_strings[i].size() - 1; j++) {
                        logger.debug() << complex_strings[i][j] << endl;
                    }
                    logger.debug() << endl;
                    logger.debug() << consensus << endl;
                }
                i += complex_regions[cur_complex_ind].second;
                if (cur_complex_ind + 1< complex_regions.size())
                    cur_complex_ind ++;
            } else {
                int real_cov = 1;
                int cov = 1;
                if (quantity[i] == 0) {
                    zero_covered++;
                } else {
                    size_t real_len =  amounts[i][0];
                    sort(amounts[i].begin() + 1, amounts[i].begin() + 1 + real_len);
                    real_cov = amounts[i][(1 + real_len)/2];

                    cov = int(round(sum[i] * 1.0 / quantity[i]));
                    if (real_cov != cov) {
                        cov = real_cov;
                    }
                }

                quantities[quantity[i]]++;
                for (size_t j = 0; j < cov; j++) {
                    ss << sequence[i];
                    total_count++;
                }
                i++;
            }
        }
        logger.trace() <<"Contig " << name << " uncompressed length: " << sequence.length() << " processed." << endl;
        logger.trace() << "Zero covered (after filtering) " << zero_covered << endl;
//Constants;
        for (size_t i = 0; i < 20; i++) {
            logger.debug() << i << " " << quantities[i] << endl;
        }
        return ss.str();
    }
};

struct AssemblyInfo {
    std::map<string, ContigInfo> contigs;
//dinucleotide repeats of larger length will be compressed
    size_t compression_length;

    static const size_t SW_BANDWIDTH = 10;

//We do not believe matches on the ends of match region, DO WE?
    static const size_t MATCH_EPS = 0;

    static const size_t BATCH_SIZE = 100000;

    explicit AssemblyInfo (logging::Logger &logger,
                           const std::experimental::filesystem::path &contigs_file,
                           size_t dicompress) {
//TODO paths;
        std::vector<Contig> assembly = io::SeqReader(contigs_file).readAllContigs();

        for (const auto& contig: assembly) {
//TODO switch to Contig()
            contigs.emplace(contig.id, ContigInfo(contig.seq.str(), contig.id));
            logger.debug() << contig.id << endl;
        }
        compression_length= dicompress;
    }

//    void AddRead(const string& contig_name){
//        contigs[contig_name].AddRead();
//    }


    AlignmentInfo readAlignment(std::ifstream &ss){
        AlignmentInfo res;
        if (ss.eof())
            return res;
        ss >> res.read_id >> res.read_start >> res.read_end >>res.contig_id >> res.alignment_start>>res.alignment_end;

        if (res.contig_id[0] == '-') {
            res.contig_id = res.contig_id.substr(1);
            res.rc = true;
        } else {
            res.rc = false;
        }

        return res;
    }

    void RC(AlignmentInfo & aln, size_t read_len){
        size_t contig_len = contigs[aln.contig_id].len;
        size_t tmp = contig_len - aln.alignment_end;
        aln.alignment_end = contig_len - aln.alignment_start;
        aln.alignment_start = tmp;
        aln.read_end = std::min(aln.read_end, read_len);
        tmp = read_len - aln.read_end;
        aln.read_end = read_len - aln.read_start;
        aln.read_start = tmp;

    }

    size_t getUncompressedPos(const string &s, size_t pos){
        size_t real_pos = 0;
        while (pos > 0) {
            while (s[real_pos] == s[real_pos + 1]) {
                real_pos ++;
            }
            pos --;
            real_pos ++;
        }
        return real_pos;
    }

    string uncompress (size_t from, size_t to, const string &s) {
//        logger.info() << from << " " << to << " " << s.length() << endl;
        size_t real_from = getUncompressedPos(s, from);
//can be optimised
        size_t real_to = getUncompressedPos(s, to);
        VERIFY(real_to >= real_from);
        return s.substr(real_from,  real_to -real_from);
    }

    string uncompressCoords (size_t from, size_t to, const string &s, const vector<size_t> & coords) {
//        logger.info() << from << " " << to << " " << s.length() << endl;
        size_t real_from = coords[from];
//can be optimised
        size_t real_to = coords[to];
        VERIFY_MSG(real_to >= real_from, "Something went WRONG in uncompression");
        return s.substr(real_from,  real_to -real_from);
    }


    //Some strange cigars are output
    bool verifyCigar(vector<cigar_pair> &cigars, int bandwidth ) {
//TODO constant??
        if (cigars.size() > 200) {
            return false;
        }
        int shift = 0;
        for(size_t i = 0; i + 1 < cigars.size(); i++) {
            if (cigars[i].type != 'M' && cigars[i+1].type != 'M') {
                return false;
            }
            if (cigars[i].type == 'D') {
                shift -= cigars[i].length;
            } else if (cigars[i].type == 'I') {
                shift += cigars[i].length;
            }
            if (shift > bandwidth || shift < -bandwidth)
                return false;
        }
        return true;
    }

    string str(vector<cigar_pair> &cigars) {
        std::stringstream ss;
        for(size_t i = 0; i < cigars.size(); i++) {
            ss << cigars[i].length << cigars[i].type;
        }
        return ss.str();
    }

    size_t matchedLength(vector<cigar_pair> &cigars) {
        size_t res = 0;
        for(size_t i = 0; i < cigars.size(); i++) {
            if (cigars[i].type == 'M')
                res += cigars[i].length;
        }
        return res;
    }

    std::vector<cigar_pair> getFastAln(logging::Logger &logger, AlignmentInfo& aln, const char * contig, const char *read) {

        size_t cur_bandwidth = SW_BANDWIDTH;
//strings, match, mismatch, gap_open, gap_extend, width
        auto cigars = align_ksw(contig, read, 1, -5, 5, 2, cur_bandwidth);
        auto str_cigars = str(cigars);
        size_t matched_l = matchedLength(cigars);
        bool valid_cigar = true;
//TODO: consts
        while ((matched_l < strlen(read) * 0.9 || !(valid_cigar = verifyCigar(cigars, cur_bandwidth)))) {
//Do we really need this?
            if (matched_l < 50) {
                logger.debug() << aln.read_id << " ultrashort alignmnent, doing nothing" << endl;
                logger.debug() <<str_cigars<< endl;
//                logger.trace() << string(contig) << endl;
//                logger.trace() << string (read) << endl;
                break;
            } else {
//We allow large indels now, so we'll have to increase bandwidth significantly to overcome them, otherwise iterative process will be stopped
                if (cur_bandwidth == 10)
                    cur_bandwidth *= 5;
                else
                    cur_bandwidth *= 2;
                if (cur_bandwidth > 260) {
                    break;
                }
                logger.debug() << aln.read_id << endl << str(cigars) << endl << "aln length " << aln.length()
                               << " read length " << strlen(read)
                               << " matched length " << matched_l << endl;
                cigars = align_ksw(contig, read, 1, -5, 5, 2, cur_bandwidth);
                size_t new_matched_len = matchedLength(cigars);
                logger.debug() << aln.read_id << " alignment replaced using bandwindth " << cur_bandwidth << endl
                               << str(cigars) << endl;
                if (new_matched_len <= matched_l && valid_cigar) {
                    logger.debug() << aln.read_id << " alignmnent length did not improve after moving to bandwidth " << cur_bandwidth << endl;
                    matched_l = new_matched_len;
                    break;
                }
                matched_l = new_matched_len;
            }
        }
        return cigars;
    }

    string compressRead(const string& read, vector<size_t>& uncompressed_positions) {
//TODO:: is it fast?
//        logger.trace() << read << endl;
        string res;
        size_t current_coord = 0;
        uncompressed_positions.resize(0);
        for (size_t i = 0; i < read.length(); i++) {
            bool compressed = false;
            if (i > 0 && read[i] == read[i - 1]) {
                compressed = true;
            } else if (current_coord > 2*compression_length) {
                bool in_repeat = (read[i] ==res[current_coord - 2] || read[i] == res[current_coord - 1]);
                if (in_repeat) {
                    for (size_t j = 1; j < compression_length; j++) {
                        if (res[current_coord  - 2 * j] != res[current_coord -2]) {
                            in_repeat = false;
                            break;
                        }
                    }
                }
                if (in_repeat) {
                    for (size_t j = 1; j < compression_length; j++) {
                        if (res[current_coord -1 - 2 * j] != res[current_coord -1]) {
                            in_repeat = false;
                            break;
                        }
                    }
                }
                compressed |= in_repeat;
            }
            if (!compressed) {
                res = std::move(res) + read[i];
                uncompressed_positions.push_back(i);
                current_coord ++;
            }
        }
//        logger.trace() << res.str() << endl;
        return res;
    }

    void processReadPair (logging::Logger &logger, string& read, AlignmentInfo& aln) {
//        logger.info() << read.id << endl;
        if (contigs.find(aln.contig_id) == contigs.end())
            return;
        Sequence uncompressed_read_seq (read);
        size_t rlen = read.length();
        vector<size_t> compressed_read_coords;
        string compressed_read;
        if (aln.rc) {
            uncompressed_read_seq = !uncompressed_read_seq;
            compressed_read = compressRead(uncompressed_read_seq.str(), compressed_read_coords);
            RC(aln, compressed_read.length());
//            read = read.RC();
        } else {
            compressed_read = compressRead(uncompressed_read_seq.str(), compressed_read_coords);
        }
        logger.trace() << aln.read_id << " "<<  aln.alignment_start << " " << aln.alignment_end << endl;
        ContigInfo& current_contig = contigs[aln.contig_id];
//        compressed_read.erase(std::unique(compressed_read.begin(), compressed_read.end()), compressed_read.end());
        string contig_seq = current_contig.sequence.substr(aln.alignment_start , aln.alignment_end - aln.alignment_start);
        string read_seq = compressed_read.substr(aln.read_start, aln.read_end - aln.read_start);
        auto cigars = getFastAln(logger, aln, contig_seq.c_str(), read_seq.c_str());
        if (matchedLength(cigars) < 50) {
            logger.trace()<< "Read " << aln.read_id << " not aligned " << endl;
            return;
        }
        int cur_ind = 0;
        std:vector<size_t> quantities;
        quantities.resize(compressed_read.length());

        for (size_t i = 0; i < compressed_read.size(); i++){
            size_t count = 0;
            cur_ind = compressed_read_coords[i];
            while (cur_ind < uncompressed_read_seq.size() && nucl(uncompressed_read_seq[cur_ind]) == compressed_read[i]) {
                count ++;
                cur_ind ++;
            }
            quantities[i] = count;
        }
        size_t cont_coords = 0; //minimapaln[0].seg_to.left;
        size_t read_coords = 0; //minimapaln[0].seg_from.left;
        read_coords = aln.read_start;
        size_t matches = 0;
        size_t mismatches = 0;
        size_t indels = 0;
        size_t complex_start = -1;
        size_t complex_fragment_finish = -1;
        size_t complex_id = -1;
        for (auto it = cigars.begin(); it != cigars.end(); ++it) {
            if ((*it).type == 'I') {
                read_coords += (*it).length;
                indels += (*it).length;
            } else if ((*it).type == 'D') {
                cont_coords += (*it).length;
//TODO do not forget add complex regions here
                indels += (*it).length;
            } else {
                size_t cur_complex_coord = -1;
                auto complex_regions_iter = lower_bound(current_contig.complex_regions.begin(),
                                                        current_contig.complex_regions.end(),
                                                        std::make_pair(cont_coords + aln.alignment_start, size_t(0)));
                if (complex_regions_iter !=  current_contig.complex_regions.end()) {
                    cur_complex_coord = complex_regions_iter->first;
                }
//TODO possibly move critical here;
#pragma omp critical
                for (size_t i = MATCH_EPS; i + MATCH_EPS< (*it).length; i++) {
                    size_t coord = cont_coords + aln.alignment_start + i;


//                    logger.info() << current_contig.sequence[coord]<< compressed_read[read_coords + i] << endl;
                    if (current_contig.sequence[coord] == nucl(compressed_read[read_coords + i]) && current_contig.quantity[coord] != 255 && current_contig.sum[coord] < 60000) {
                        {
                            current_contig.quantity[coord]++;
                            current_contig.sum[coord] += quantities[read_coords + i];
                            if (current_contig.amounts[coord][0] < ContigInfo::VOTES_STORED)
//#pragma omp critical
                                current_contig.amounts[coord][++current_contig.amounts[coord][0]] = quantities[
                                        read_coords + i];
                        }
                        matches ++;
                    } else {
                        mismatches ++;
                    }
//Only complete traversion of complex regions taken in account;
//                    if (current_contig.complex_regions.find(coord) != current_contig.complex_regions.end()) {
                    if (coord == cur_complex_coord) {
                        size_t complex_len = complex_regions_iter->second;
                        if (read_coords + complex_len < matchedLength(cigars)) {
                            complex_start = read_coords + i;
                            complex_fragment_finish = coord + complex_len;
                            complex_id = coord;
                        }
                        complex_regions_iter ++;
                        if (complex_regions_iter != current_contig.complex_regions.end())
                            cur_complex_coord = complex_regions_iter->first;
                    }
                    if (coord == complex_fragment_finish) {
//#pragma omp critical
                        current_contig.complex_strings[complex_id].push_back(uncompressCoords(complex_start, read_coords + i, uncompressed_read_seq.str(), compressed_read_coords));
                        complex_fragment_finish = -1;
                    } else if (coord > complex_fragment_finish) {
                        logger.debug() << "Read " << aln.read_id << " missed fragment finish " << complex_fragment_finish << endl;
                        complex_fragment_finish = -1;
                    }
                }
                read_coords += (*it).length;
                cont_coords += (*it).length;
            }
        }
        if (matches < mismatches * 3)
            logger.debug()<<"For read too much mismatches " << aln.read_id << " matches/MM: " << matches << "/" << mismatches << endl;
    }

    void processBatch(logging::Logger &logger, vector<string>& batch, vector<AlignmentInfo>& alignments){
        size_t len = batch.size();
#pragma omp parallel for default(none) shared(logger, len, batch, alignments)
        for (size_t i = 0; i < len; i++) {
            processReadPair(logger, batch[i], alignments[i]);
        }
    }

    vector<Contig> process(logging::Logger &logger, const io::Library &lib,
                           const std::experimental::filesystem::path &alignmens_file) {
        std::ifstream compressed_reads;
        std::ofstream corrected_contigs;
        compressed_reads.open(alignmens_file);
        io::SeqReader reader(lib);
        logger.trace() << "Initialized\n";
        if (compressed_reads.eof()) {
            logger.info() << "NO ALIGNMENTS AVAILABLE!";
            exit(1);
        }
        AlignmentInfo cur_align = readAlignment(compressed_reads);
        string cur_compressed = cur_align.read_id;
        string cur_read;
        string cur_seq;
        StringContig cur;
        logger.info() << "Reading and processing initial reads from " << lib << "\n";
        size_t reads_count = 0;
        size_t aln_count = 1;
        vector<AlignmentInfo> align_batch;
        vector<string> contig_batch;
        while (!compressed_reads.eof()) {
            bool reads_over = false;
            while (cur.id != cur_compressed) {
                if (reader.eof()) {
                    logger.info() << "Reads are over\n";
                    reads_over = true;
                    break;
                }
                cur = reader.read();
//Some tools clip read names after whitespaces
                cur.id = split(cur.id)[0];
                reads_count ++;
                if (reads_count % 1000 == 0) {
                    logger.trace() << "Processed " << reads_count << " original reads " << endl;
                }
            }
            if (reads_over) {
                break;
            }
            align_batch.push_back(cur_align);
            contig_batch.push_back(cur.seq);
//TODO:: appropriate logic for multiple alignment
            do {
                cur_align = readAlignment(compressed_reads);
                aln_count ++;
                if (aln_count % BATCH_SIZE == 0) {
                    logger.trace() << "Batch of size " <<BATCH_SIZE <<" created, processing" << endl;
                    processBatch(logger, contig_batch, align_batch);

                    logger.trace() << "Processed " << aln_count << " compressed mappings " << endl;
                    contig_batch.resize(0);
                    align_batch.resize(0);
                    //exit(0);
                }

                if (cur_compressed == cur_align.read_id) {
                    align_batch.push_back(cur_align);
                    contig_batch.push_back(cur.seq);
                }
            } while (cur_compressed == cur_align.read_id);
            cur_compressed = cur_align.read_id;
        }
        processBatch(logger, contig_batch, align_batch);
        logger.trace() << "Processed final batch of " << align_batch.size() << " compressed reads " << endl;
        vector<Contig> res;
        logger.info() << "Uncompressing homopolymers in contigs" << endl;
        for (auto& contig: contigs){
            logger.trace() << "Generating consensus for contig " << contig.first << endl;
            res.emplace_back(Sequence(contig.second.GenerateConsensus(logger)), contig.first);
        }
        size_t total_zero_covered = 0;
        for (auto & contig: contigs) {
            total_zero_covered += contig.second.zero_covered;
        }
        logger.info() << "Total zero covered nucleotides "  << total_zero_covered << endl;
        return std::move(res);
    }
};


std::experimental::filesystem::path Polish(logging::Logger &logger, size_t threads,
                                           const std::experimental::filesystem::path &output_file,
                                           const std::experimental::filesystem::path &contigs_file,
                                           const std::experimental::filesystem::path &alignments,
                                           const io::Library &reads, size_t dicompress) {
    omp_set_num_threads(threads);
    AssemblyInfo assemblyInfo(logger, contigs_file, dicompress);
    std::vector<Contig> res = assemblyInfo.process(logger, reads, alignments);
    std::ofstream os;
    os.open(output_file);
    for(Contig &contig : res) {
        os << ">" << contig.getId() << "\n" << contig.seq << "\n";
    }
    os.close();
    return output_file;
}

