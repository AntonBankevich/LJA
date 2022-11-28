#pragma once

#include "multi_graph.hpp"
#include "gap_closing.hpp"
#include "polishing/perfect_alignment.hpp"
#include "error_correction/mult_correction.hpp"
#include "error_correction/mitochondria_rescue.hpp"
#include "error_correction/manyk_correction.hpp"
#include "error_correction/precorrection.hpp"
#include "sequences/seqio.hpp"
#include "dbg/dbg_construction.hpp"
#include "common/rolling_hash.hpp"
#include "common/dir_utils.hpp"
#include "common/cl_parser.hpp"
#include "common/logging.hpp"

#include <iostream>
#include <queue>
#include <wait.h>
#include <error_correction/dimer_correction.hpp>
#include "common/rolling_hash.hpp"
#include <utility>
#include <ksw2/ksw_wrapper.hpp>

namespace pipeline {

    struct LJAPipeline {

        std::vector<Contig> ref;
        size_t stage_num;

        LJAPipeline(io::Library ref_lib) {
            ref = io::SeqReader(ref_lib).readAllContigs();
            stage_num = 0;
        }

    };
}