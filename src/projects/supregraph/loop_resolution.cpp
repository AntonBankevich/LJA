#include "loop_resolution.hpp"
#include <cmath>

using namespace ag;
using spg::Supregraph;
using spg::MergeUnbranching;

namespace {
    bool Distinct(const std::vector<Vertex *> &vertices) {
        for (size_t i = 0; i < vertices.size(); i++)
            for (size_t j = i + 1; j < vertices.size(); j++)
                if (*vertices[i] == *vertices[j] || *vertices[i] == vertices[j]->rc())
                    return false;
        return true;
    }

//    Read support of the pair in -> out at in's finish, if the suffix tracker keeps a record for in.
    std::optional<size_t> PairSupport(SuffixTracker &suffixes, Edge &in, Edge &out) {
        if (!suffixes.hasSuffixRecord(in))
            return {};
        return suffixes.getSuffixRecord(in).countStartsWith(GraphPath(out));
    }

    double Coverage(const Vertex &v) {
        return v.coverage_info.weight > 0 ? v.getSPGCoverage() : -1;
    }

    VertexId FindVertex(Supregraph &spg, Vertex::id_type id) {
        for (Vertex &v : spg.vertices())
            if (v.getId().innerId() == id)
                return v.getId();
        return {};
    }

//    Cleans the graph around a resolution step the way the multiplexer does: the unbranching paths through the created
//    vertices and through their non-junction neighbours are merged (MergeUnbranching: merge, then flips). Returns the
//    vertex the first created vertex was merged into. Ids, not handles: merges remove vertices.
    VertexId MergeAround(logging::Logger &logger, Supregraph &spg, const std::vector<VertexId> &created) {
        std::vector<Vertex::id_type> seeds;
        for (VertexId v : created)
            seeds.push_back(v.innerId());
        for (VertexId v : created) {
            for (Edge &e : *v)
                seeds.push_back(e.getFinish().getId().innerId());
            for (Edge &e : v->incoming())
                seeds.push_back(e.getStart().getId().innerId());
        }
        Vertex::id_type first = seeds.empty() ? 0 : seeds.front();
        for (size_t i = 0; i < seeds.size(); i++) {
            VertexId vid = FindVertex(spg, seeds[i]);
            if (vid.valid() && !vid->isJunction()) {
                Vertex::id_type res = MergeUnbranching(logger, spg, *vid).getId().innerId();
                spg.removeMarked();
                if (i == 0)
                    first = res;
            }
        }
        return seeds.empty() ? VertexId() : FindVertex(spg, first);
    }

//    Resolves the core of loop: the entrance continues into the body and the body into the exit; with keep_loop the
//    body also continues into itself, a single new vertex. After merging around the step, that vertex's merge is the
//    new body, and its neighbour the new core (returned if the loop is simple).
    std::optional<spg::SimpleLoop> UnrollStep(logging::Logger &logger, Supregraph &spg, const spg::SimpleLoop &loop,
                                              bool keep_loop) {
        Vertex &core = *loop.core;
        VertexResolutionPlan plan(core);
        plan.add(*loop.e_in, *loop.l_out);
        if (keep_loop)
            plan.add(*loop.l_in, *loop.l_out);
        plan.add(*loop.l_in, *loop.e_out);
        VertexResolutionResult result = spg.resolveVertex(core, plan);
        std::vector<VertexId> created;
        for (Vertex &v : result.newVertices())
            created.push_back(v.getId());
        spg.removeMarked();
//        With keep_loop the entrance and exit pairs share an edge with the loop pair, so resolveVertex connects them by
//        direct edges; the loop pair is the only new vertex.
        VERIFY(!keep_loop || created.size() == 1);
        VertexId body = MergeAround(logger, spg, created);
        if (!keep_loop || !body.valid())
            return {};
        std::optional<spg::SimpleLoop> next = spg::DetectSimpleLoop(body->frontVertex());
        if (next && next->body != body)
            return {};
        return next;
    }
}

std::optional<spg::SimpleLoop> spg::DetectSimpleLoop(Vertex &core) {
    if (core.marked() || core.inDeg() != 2 || core.outDeg() != 2)
        return {};
    for (Edge &l_out : core) {
        Vertex &body = l_out.getFinish();
        if (body.inDeg() == 1 && body.outDeg() == 1 && body.frontVertex() == core) {
            Edge &l_in = body.front();
            Edge *e_in = nullptr;
            Edge *e_out = nullptr;
            for (Edge &e : core.incoming())
                if (e != l_in)
                    e_in = &e;
            for (Edge &e : core)
                if (e != l_out)
                    e_out = &e;
            if (e_in == nullptr || e_out == nullptr)
                return {};
            Vertex &u = e_in->getStart();
            Vertex &w = e_out->getFinish();
            if (core == core.rc() || body == body.rc() || !Distinct({&core, &body, &u, &w}))
                return {};
            return SimpleLoop{core.getId(), body.getId(), e_in->getId(), l_out.getId(), l_in.getId(), e_out->getId()};
        }
    }
    return {};
}

double spg::EstimateLoopCopies(const SimpleLoop &loop) {
    double single = (Coverage(loop.e_in->getStart()) + Coverage(loop.e_out->getFinish())) / 2;
    if (Coverage(loop.e_in->getStart()) <= 0 || Coverage(loop.e_out->getFinish()) <= 0)
        return -1;
    double body_estimate = Coverage(*loop.body) / single;
    double core_estimate = Coverage(*loop.core) / single - 1;
    double body_len = loop.body->size();
    double core_len = loop.core->size();
    return (body_estimate * body_len + core_estimate * core_len) / (body_len + core_len);
}

size_t spg::ResolveLoops(logging::Logger &logger, Supregraph &spg, SuffixTracker &suffixes, size_t max_copies) {
    logger.info() << "Resolving simple loops" << std::endl;
//    Ids, not handles: unrolling one loop removes vertices, possibly a later candidate.
    std::vector<Vertex::id_type> candidates;
    for (Vertex &v : spg.verticesUnique())
        if (DetectSimpleLoop(v))
            candidates.push_back(v.getId().innerId());
    size_t resolved = 0;
    for (Vertex::id_type id : candidates) {
        VertexId vid = FindVertex(spg, id);
        std::optional<SimpleLoop> loop;
        if (vid.valid())
            loop = DetectSimpleLoop(*vid);
        if (!loop) {
            logger.trace() << "Loop at " << id << " is gone" << std::endl;
        } else {
            double estimate = EstimateLoopCopies(*loop);
            size_t copies = estimate < 0.5 ? 0 : size_t(std::lround(estimate));
            logger.trace() << "Loop at " << id << ": core " << loop->core->size() << "bp cov " << Coverage(*loop->core)
                           << ", body " << loop->body->getId() << " " << loop->body->size() << "bp cov " << Coverage(*loop->body)
                           << ", neighbours cov " << Coverage(loop->e_in->getStart()) << " " << Coverage(loop->e_out->getFinish())
                           << ", estimated copies " << estimate << std::endl;
            if (copies == 0 || copies > max_copies) {
                logger.trace() << "Loop at " << id << " left: estimated copies out of range" << std::endl;
                loop.reset();
            }
            size_t made = 0;
            while (loop) {
                std::optional<size_t> bypass = PairSupport(suffixes, *loop->e_in, *loop->e_out);
                std::optional<size_t> around = PairSupport(suffixes, *loop->l_in, *loop->l_out);
                if (!bypass || !around) {
                    logger.trace() << "Loop at " << id << " left after " << made << " copies: no read record" << std::endl;
                    loop.reset();
                } else if (*bypass > 0) {
                    logger.trace() << "Loop at " << id << " left after " << made << " copies: reads bypass it" << std::endl;
                    loop.reset();
                } else if (made + 1 > max_copies) {
                    logger.trace() << "Loop at " << id << " left after " << made << " copies: too many copies" << std::endl;
                    loop.reset();
                } else {
                    bool keep_loop = made + 1 < copies || *around > 0;
                    if (keep_loop && made + 1 >= copies)
                        logger.trace() << "Loop at " << id << ": reads go around it more than " << copies << " times" << std::endl;
                    made++;
                    loop = UnrollStep(logger, spg, *loop, keep_loop);
                    if (!keep_loop) {
                        logger.trace() << "Loop at " << id << " unrolled into " << made << " copies" << std::endl;
                        resolved++;
                    } else if (!loop) {
                        logger.trace() << "Loop at " << id << " left after " << made << " copies: lost its shape" << std::endl;
                    }
                }
            }
        }
    }
    logger.info() << "Unrolled " << resolved << " of " << candidates.size() << " simple loops" << std::endl;
    return resolved;
}
