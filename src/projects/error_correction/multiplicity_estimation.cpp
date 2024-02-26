#include <common/simple_computation.hpp>
#include "diploidy_analysis.hpp"
#include "multiplicity_estimation.hpp"
#include "correction_utils.hpp"

using namespace dbg;
size_t BoundRecord::inf = 1000000000000ul;

MappedNetwork::MappedNetwork(const Component &component, const std::function<bool(const dbg::Edge &)> &unique,
                             double rel_coverage, double unique_coverage, double double_coverage) {
    dbg::SparseDBG &graphr = component.graph();
    for(dbg::Vertex &v : component.vertices()) {
        vertex_mapping[&v] = addVertex();
    }
    std::vector<std::pair<int, dbg::Edge *>> edgeids;
    for(dbg::Vertex &v : component.vertices()) {
        for(dbg::Edge &edge : v) {
            if (!unique(edge)) {
                size_t min_flow = (!edge.is_reliable || (edge.getCoverage() < rel_coverage && edge.truncSize() < 20000)) ? 0 : 1;
                size_t max_flow = 1000000000;
                if(edge.truncSize() > 1000 && edge.getCoverage() < double_coverage)
                    max_flow = 2;
                if(edge.truncSize() > 1000 && edge.getCoverage() < unique_coverage)
                    max_flow = 1;
                int eid = addEdge(vertex_mapping[&v], vertex_mapping[&edge.getFinish()], min_flow, max_flow);
                VERIFY(eid > 0);
                edgeids.emplace_back(eid, &edge);
            } else {
                addSink(vertex_mapping[&v], 1);
                addSource(vertex_mapping[&v.rc()], 1);
            }
        }
    }
    if(edgeids.empty())
        return;
    edge_mapping.resize(std::max_element(edgeids.begin(), edgeids.end())->first + 1, nullptr);
    for(std::pair<int, dbg::Edge *> &p : edgeids) {
        edge_mapping[p.first] = p.second;
    }
}

std::vector<dbg::Edge *> MappedNetwork::getUnique(logging::Logger &logger) {
    std::vector<dbg::Edge*> res;
    std::unordered_map<int, size_t> multiplicities = findFixedMultiplicities();
    for (auto &rec : multiplicities) {
//        logger << "Edge " << edge_mapping[rec.first]->start()->hash() << edge_mapping[rec.first]->start()->isCanonical()
//               << "ACGT"[edge_mapping[rec.first]->seq[0]]
//               << " has fixed multiplicity " << rec.second << std::endl;
        if(rec.second == 1 && rec.first < edge_mapping.size() && edge_mapping[rec.first] != nullptr)
            res.emplace_back(edge_mapping[rec.first]);
    }
    return std::move(res);
}

std::unordered_map<dbg::Edge *, std::pair<size_t, size_t>> MappedNetwork::findBounds() {
    std::unordered_map<dbg::Edge *, std::pair<size_t, size_t>> res;
    for(const auto &rec : Network::findBounds()) {
        if(rec.first < edge_mapping.size() && edge_mapping[rec.first] != nullptr) {
            res.emplace(edge_mapping[rec.first], rec.second);
        }
    }
    return std::move(res);
}

void UniqueClassificator::markPseudoHets() const {
    for(Edge &edge : dbg.edges()) {
        Vertex &start = edge.getFinish();
        if(!isUnique(edge) || start.outDeg() != 2 || start.inDeg() != 1)
            continue;
        Edge &correct = start.front().getCoverage() > start.back().getCoverage() ? start.front() : start.back();
        Edge &incorrect = start.front().getCoverage() <= start.back().getCoverage() ? start.front() : start.back();
        incorrect.is_reliable = false;
        incorrect.rc().is_reliable = false;
        dbg::GraphPath cor_ext = reads_storage.getRecord(start).
                getFullUniqueExtension(correct.truncSeq().Subseq(0, 1), 1, 0).unpack();
        dbg::GraphPath incor_ext = reads_storage.getRecord(start).
                getFullUniqueExtension(incorrect.truncSeq().Subseq(0, 1), 1, 0).unpack();
        for(Segment<Edge> seg : incor_ext) {
            bool found= false;
            for(Segment<Edge> seg1 : cor_ext) {
                if(seg.contig() == seg1.contig()) {
                    found = true;
                    break;
                }
            }
            if(found)
                break;
            seg.contig().is_reliable = false;
            seg.contig().rc().is_reliable = false;
        }
    }
}

void UniqueClassificator::classify(logging::Logger &logger, size_t unique_len,
                                   const std::experimental::filesystem::path &dir) {
    logger.info() << "Looking for unique edges" << std::endl;
    if(debug)
        recreate_dir(dir);
    size_t cnt = 0;
    for(Edge &edge : dbg.edges()) {
        edge.is_reliable = true;
    }
    if(diploid) {
        SetUniquenessStorage duninque = BulgePathFinder(dbg).uniqueEdges(unique_len);
        for (Edge &edge : dbg.edges()) {
            if (duninque.isUnique(edge)) {
                updateBounds(edge, 1, 1);
                cnt++;
            }
        }
    } else {
        for (Edge &edge : dbg.edges()) {
            if(this->isUnique(edge))
                continue;
            dbg::GraphPath al = FindLongestCoveredExtension(edge, 3, 1);
            if(al.truncLen() > unique_len) {
                for(Segment<Edge> seg : al) {
                    updateBounds(seg.contig(), 1, 1);
                }
            } else if(edge.getStart().inDeg() == 0 && edge.truncSize() > unique_len / 3) {
                updateBounds(edge, 1, 1);
                cnt++;
            }
        }
    }
    logger.info() << "Marked " << cnt << " long edges as unique" << std::endl;
    logger.info() << "Marking extra edges as unique based on read paths" << std::endl;
    std::vector<Edge *> extra_unique;
    for(Edge &edge : dbg.edges()) {
        if(isUnique(edge) || edge.getFinish().outDeg() > 1) {
            continue;
        }
        const VertexRecord &rec = reads_storage.getRecord(edge.getStart());
        CompactPath unique_extension = rec.getFullUniqueExtension(edge.truncSeq().Subseq(0, 1), 1, 0);
        dbg::GraphPath path = unique_extension.unpack();
        size_t len = 0;
        for(size_t i = 1; i < path.size(); i++) {
            if(isUnique(path[i].contig())) {
                path = path.subPath(0, i + 1);
                break;
            }
            len += path[i].size();
        }
        if(!isUnique(path.backEdge()) ||len > 3000 || rec.countStartsWith(CompactPath(path).cpath()) < 4)
            continue;
        CompactPath back_unique = reads_storage.getRecord(path.finish().rc()).getFullUniqueExtension(
                path.backEdge().rc().truncSeq().Subseq(0, 1), 1, 0);
        if(back_unique.size() >= path.size()) {
            extra_unique.emplace_back(&edge);
            cnt++;
            logger.trace() << "Found extra unique edge " << edge.getInnerId() << " " << edge.truncSize() << " " << edge.getCoverage() << std::endl;
        }
    }
    for(Edge *eit : extra_unique) {
        updateBounds(*eit, 1, 1);
    }
    logger.info() << "Marked " << cnt << " edges as unique" << std::endl;
    logger.trace() << "Marking bulges to collapse" << std::endl;
    markPseudoHets();
    logger.info() << "Splitting graph with unique edges" << std::endl;
    std::vector<Component> split = UniqueSplitter(*this).split(Component(dbg));
    logger.info() << "Processing " << split.size() << " components" << std::endl;
    size_t component_cnt = 0;
    for(Component &component : split) {
        component_cnt += 1;
        if(debug)
            printDot(dir / (std::to_string(component_cnt) + ".dot"), component, reads_storage.labeler());
        //TODO make parallel trace
        logger.trace() << "Component parameters: size=" << component.uniqueSize() << " border=" << component.countBorderEdges() <<
                      " tips=" << component.countTips() <<
                      " subcomponents=" << component.realCC() << " acyclic=" << component.isAcyclic() <<std::endl;
        if(component.uniqueSize() > 2 && component.countBorderEdges() == 2 &&component.countTips() == 0 &&
           component.realCC() == 2 && component.isAcyclic()) {
            processSimpleComponent(logger, component);
        }
        cnt += processComponent(logger, component);
        if(debug) {
            logger.trace() << "Printing component to " << (dir / (std::to_string(component_cnt) + ".dot")) << std::endl;
            printDot(dir / (std::to_string(component_cnt) + ".dot"), component,
                     this->labeler() + reads_storage.labeler(), this->colorer());
        }
    }
    logger.info() << "Finished unique edges search. Found " << cnt << " unique edges" << std::endl;
    logger.info() << "Analysing repeats of multiplicity 2 and looking for additional unique edges" << std::endl;
    std::function<bool(const dbg::Edge &)> mult2 = [this](const dbg::Edge &edge) {
        return MultiplicityBounds::lowerBound(edge) == 2 && MultiplicityBounds::lowerBound(edge) == 2;
    };
    split = ConditionSplitter(mult2).splitGraph(dbg);
    cnt = 0;
    for(Component &component : split) {
        if(component.uniqueSize() != 4 || !component.isAcyclic() || component.realCC() != 2 || component.borderVertices().size() != 2) {
            continue;
        }
        logger.trace() << "Found suspicious repeat of multiplicity 2: ";
        for(Vertex &vertex : component.verticesUnique()) {
            logger << " " << vertex.getShortId();
        }
        logger << std::endl;
        bool success = processSimpleRepeat(component);
        if(success) {
            cnt++;
            logger.trace() << "Found erroneous edge" << std::endl;
        }
    }
    logger.info() << "Finished processing of repeats of multiplicity 2. Found " << cnt << " erroneous edges." << std::endl;
}

bool UniqueClassificator::processSimpleRepeat(const Component &component) {
    std::vector<Vertex *> border = component.borderVertices();
    VERIFY(border.size() == 2);
    Vertex &start = border[0]->rc();
    Vertex &end = *border[1];
    if(start.outDeg() !=2 || end.inDeg() != 2)
        return false;
    std::vector<Edge *> bad_candidates;
    if(start.front().getStart() == end || start.back().getStart() == end) {
        size_t ind = start.front().getStart() == end ? 0 : 1;
        Edge &e1 = start.front().getStart() == end ? start.front() : start.back();
        Edge &e21 = start.front().getStart() == end ? start.back() : start.front();
        if(e21.getFinish().outDeg() != 2)
            return false;
        Edge &e221 = e21.getFinish().front();
        Edge &e222 = e21.getFinish().back();
        if(e221.getFinish() != e222.getFinish())
            return false;
        if(e221.getFinish().outDeg() != 1)
            return false;
        Edge &e23 = e221.getFinish().front();
        if(e23.getFinish() != end)
            return false;
        if(!component.contains(e21.getFinish()) || !component.contains(e23.getStart()))
            return false;
        bad_candidates = {&e221, &e222};
    } else {
        size_t ind1 = start.front().truncSize() > start.back().truncSize() ? 0 : 1;
        Edge &e11 = start.front().truncSize() > start.back().truncSize() ? start.front() : start.back();
        Edge &e12 = start.front().truncSize() > start.back().truncSize() ? start.back() : start.front();
        VERIFY(e11.truncSize() >= e12.truncSize());
        size_t ind2 = end.rc().front().truncSize() > end.rc().back().truncSize() ? 1 : 0;
        Edge &e31 = end.rc().front().truncSize() > end.rc().back().truncSize() ? end.rc().back().rc() : end.rc().front().rc();
        Edge &e32 = end.rc().front().truncSize() > end.rc().back().truncSize() ? end.rc().front().rc() : end.rc().back().rc();
        VERIFY(e31.truncSize() <= e32.truncSize());
        if(e12.getFinish().outDeg() != 2 || e31.getStart().inDeg() != 2 || e11.getFinish() != e31.getStart() ||
                e12.getFinish() != e32.getStart())
            return false;
        Edge &e2 = e12.getFinish().front().getFinish() == e31.getStart() ? e12.getFinish().front() : e12.getFinish().back();
        if(!component.contains(e2.getStart()) || ! component.contains(e2.getFinish()))
            return false;
        bad_candidates = {&e11, &e2, &e32};
    }
    std::function<double(Edge* const &)> f = [](Edge * const &edge){return edge->getCoverage();};
    Edge &bad = *bad_candidates[MinIndex<Edge *, double>(bad_candidates, f)];
    std::unordered_set<Edge *> nonUnique = {&bad, &bad.rc()};
    if(bad.getStart() == start) {
        nonUnique.emplace(&start.front());
        nonUnique.emplace(&start.back());
        nonUnique.emplace(&start.front().rc());
        nonUnique.emplace(&start.back().rc());
    }
    if(bad.getFinish() == end) {
        nonUnique.emplace(&end.rc().front());
        nonUnique.emplace(&end.rc().back());
        nonUnique.emplace(&end.rc().back().rc());
        nonUnique.emplace(&end.rc().back().rc());
    }
    for(Edge &edge : component.edgesInner()) {
        if(nonUnique.find(&edge) == nonUnique.end())
            updateBounds(edge, 1, 1);
    }
    return true;
}

size_t UniqueClassificator::ProcessUsingCoverage(logging::Logger &logger,
                                          const Component &subcomponent,
                                          const std::function<bool(const dbg::Edge &)> &is_unique,
                                          double rel_coverage) {
    size_t ucnt = 0;
    std::pair<double, double> tmp = minmaxCov(subcomponent, reads_storage, is_unique);
    double min_cov = tmp.first;
    double max_cov = tmp.second;
    double threshold = std::max(min_cov * 1.4, max_cov * 1.2);
    double double_threshold = 2.4 * min_cov;
    double adjusted_rel_coverage = std::min(min_cov * 0.9, max_cov * 0.7);
    logger.trace() << "Processing component: ";
    for(Vertex &vertex : subcomponent.verticesUnique()) {
        logger << " " << vertex.getShortId();
    }
    logger << std::endl;
    logger.trace() << "Attempting to use coverage for multiplicity estimation with coverage threshold " << threshold << std::endl;
//    rel_coverage = 0;
    MappedNetwork net2(subcomponent, is_unique, rel_coverage, threshold);
    bool res = net2.fillNetwork();
    if (res) {
        logger.trace() << "Succeeded to use coverage for multiplicity estimation" << std::endl;
        for(auto rec : net2.findBounds()) {
            if(!MultiplicityBounds::isUnique(*rec.first) && rec.second.first == 1 && rec.second.second == 1) {
                ucnt++;
            }
            updateBounds(*rec.first, rec.second.first, rec.second.second);
        }
    } else {
        logger.trace() << "Failed to use coverage for multiplicity estimation" << std::endl;
    }
    if(!res && rel_coverage == 0) {
        logger.trace() << "Attempting to adjust reliable edge threshold from " << rel_coverage << " to " <<
                           adjusted_rel_coverage << std::endl;
        rel_coverage = adjusted_rel_coverage;
        MappedNetwork net3(subcomponent, is_unique, rel_coverage, threshold);
        res = net3.fillNetwork();
        if (res) {
            logger.trace() << "Succeeded to use coverage for multiplicity estimation" << std::endl;
            for(auto rec : net3.findBounds()) {
                if(!MultiplicityBounds::isUnique(*rec.first) && rec.second.first == 1 && rec.second.second == 1) {
                    ucnt++;
                }
                updateBounds(*rec.first, rec.second.first, rec.second.second);
            }
        } else {
            logger.trace() << "Failed to use adjusted reliable coverage for multiplicity estimation" << std::endl;
        }
    }
    if(!res) {
        logger.trace() << "Setting unique threshold to 0" << std::endl;
        threshold = 0;
        MappedNetwork net4(subcomponent, is_unique, rel_coverage, threshold);
        res = net4.fillNetwork();
        if (res) {
            logger.trace() << "Succeeded to use coverage for multiplicity estimation" << std::endl;
            for(auto rec : net4.findBounds()) {
                if(!MultiplicityBounds::isUnique(*rec.first) && rec.second.first == 1 && rec.second.second == 1) {
                    ucnt++;
                }
                updateBounds(*rec.first, rec.second.first, rec.second.second);
            }
        } else {
            logger.trace() << "Failed to use coverage for multiplicity estimation" << std::endl;
        }
    }
    if(!res)
        return ucnt;
    logger.trace() << "Attempting to set additional condition for edges of multiplicity 2" << std::endl;
    MappedNetwork net5(subcomponent, is_unique, rel_coverage, threshold, double_threshold);
    res = net5.fillNetwork();
    if(res) {
        for(auto rec : net5.findBounds()) {
            if(!MultiplicityBounds::isUnique(*rec.first) && rec.second.first == 1 && rec.second.second == 1) {
                ucnt++;
            }
            updateBounds(*rec.first, rec.second.first, rec.second.second);
        }
    } else {
        double_threshold = 0;
    }
    logger.trace() << "Final multiplicity estimation with reliable coverage threshold = " << rel_coverage
                   << " unique coverage threshold = " << threshold << " double coverage threshold = " << double_threshold << std::endl;
    return ucnt;
}

std::pair<double, double> minmaxCov(const Component &subcomponent, const RecordStorage &reads_storage,
                               const std::function<bool(const dbg::Edge &)> &is_unique) {
    double max_cov= 0;
    double min_cov= 100000;
    for(Edge &edge : subcomponent.edges()) {
        if(is_unique(edge) && subcomponent.contains(edge.getFinish())) {
            const VertexRecord & record = reads_storage.getRecord(edge.getStart());
            std::string s = edge.truncSeq().Subseq(0, 1).str();
            size_t cnt = record.countStartsWith(Sequence(s + "A")) +
                         record.countStartsWith(Sequence(s + "C")) +
                         record.countStartsWith(Sequence(s + "G")) +
                         record.countStartsWith(Sequence(s + "T"));
            if(edge.truncSize() < 20000) {
                min_cov = std::min<double>(min_cov, std::min<double>(cnt, edge.getCoverage()));
                max_cov = std::max<double>(max_cov, std::min<double>(cnt, edge.getCoverage()));
            } else {
                min_cov = std::min<double>(min_cov, cnt);
                max_cov = std::max<double>(max_cov, cnt);
            }
        }
    }
    return {min_cov, max_cov};
}

Edge &getStart(const Component &component) {
    for(Edge &edge : component.edges()) {
        if(!component.contains(edge.getFinish()))
            return edge.rc();
    }
    VERIFY(false);
}

std::vector<Vertex *> topSort(const Component &component) {
    std::vector<Vertex *> stack;
    std::vector<Vertex *> res;
    std::unordered_set<Vertex *> visited;
    stack.emplace_back(&getStart(component).getFinish());
    while(!stack.empty()) {
        Vertex *cur = stack.back();
        if(visited.find(cur) != visited.end()) {
            stack.pop_back();
            continue;
        }
        bool ok = true;
        for(Edge &e : *cur) {
            if(component.contains(e.getFinish()) && visited.find(&e.getFinish()) == visited.end())
                ok = false;
        }
        if(ok) {
            visited.emplace(cur);
            stack.pop_back();
            res.emplace_back(&cur->rc());
        } else {
            for(Edge &e : *cur) {
                VERIFY(component.contains(e.getFinish()));
                if(visited.find(&e.getFinish()) == visited.end())
                    stack.emplace_back(&e.getFinish());
            }
        }
    }
    return std::move(res);
}
void UniqueClassificator::processSimpleComponent(logging::Logger &logger, const Component &component) const {
    logger.trace() << "Collapsing acyclic component" << std::endl;
    std::vector<Vertex *> order = topSort(component);
    if(order.size() != component.uniqueSize()) {
        logger.trace() << "Failed to collapse acyclic component" << std::endl;
        return;
    }
    VERIFY(order.front()->inDeg() == 1);
    VERIFY(!component.contains(order.front()->rc().front().getFinish()));
    VERIFY(order.back()->outDeg() == 1);
    VERIFY(!component.contains(order.back()->front().getFinish()));
    std::unordered_map<Vertex *, std::pair<Edge *, size_t>> prev;
    for(Vertex *cur : order) {
        Edge *pedge = nullptr;
        size_t val = 0;
        for(Edge &edge : cur->rc()) {
            if(!component.contains(edge.getFinish())) {
                break;
            }
            VERIFY(prev.find(&edge.getFinish().rc()) != prev.end());
            size_t score = edge.intCov() + prev[&edge.getFinish().rc()].second;
            if(score > val) {
                pedge = &edge.rc();
                val = score;
            }
        }
        VERIFY(val > 0 || pedge == nullptr);
        prev[cur] = {pedge, val};
    }
    for(Edge &edge : component.edgesInner()) {
        edge.is_reliable = false;
    }
    Vertex *end = order.back();
    size_t cnt = 0;
    while(prev[end].first != nullptr) {
        prev[end].first->is_reliable = true;
        prev[end].first->rc().is_reliable = true;
        end = &prev[end].first->getStart();
        cnt++;
        VERIFY(cnt < 50);
    }
}

size_t UniqueClassificator::processComponent(logging::Logger &logger, const Component &component) {
    size_t ucnt = 0;
    logger.trace() << "Component: ";
    for(Vertex &vertex : component.verticesUnique()) {
        logger << " " << vertex.getShortId();
    }
    logger << std::endl;
    double rel_coverage = initial_rel_coverage;
    std::function<bool(const dbg::Edge &)> is_unique = [this](const dbg::Edge &edge) {
        return isUnique(edge);
    };
    MappedNetwork net(component, is_unique, rel_coverage);
    bool res = net.fillNetwork();
    if(res) {
        logger.trace() << "Found unique edges in component" << std::endl;
        for(Edge * edge : net.getUnique(logger)) {
            updateBounds(*edge, 1, 1);
            ucnt++;
        }
    } else {
        logger.trace() << "Could not find unique edges in component" << std::endl;
        logger.trace() << "Relaxing flow conditions" << std::endl;
        rel_coverage = std::max<size_t>(15, rel_coverage);
        MappedNetwork net1(component, is_unique, rel_coverage);
        res = net1.fillNetwork();
        if(res) {
            logger.trace() << "Found unique edges in component" << std::endl;
            for(Edge * edge : net1.getUnique(logger)) {
                updateBounds(*edge, 1, 1);
                ucnt++;
            }
        } else {
            logger.trace() << "Could not find unique edges with relaxed conditions in component" << std::endl;
        }
    }
    if(res) {
        std::vector<Component> subsplit = ConditionSplitter(this->asFunction()).split(component);
        logger.trace() << "Component was split into " << subsplit.size() << " subcompenents" << std::endl;
        for(Component &subcomponent : subsplit) {
            ucnt += ProcessUsingCoverage(logger, subcomponent, this->asFunction(), rel_coverage);
        }
    }
    return ucnt;
}

std::pair<Edge *, Edge *> CheckLoopComponent(const Component &component) {
    if(component.uniqueSize() != 2)
        return {nullptr, nullptr};
    Vertex *vit = &*component.vertices().begin();
    if(vit->inDeg() != 2)
        vit = &vit->rc();
    Vertex &start = *vit;
    if(start.inDeg() != 2 || start.outDeg() != 1)
        return {nullptr, nullptr};
    Edge &forward_edge = start.front();
    Vertex &end = forward_edge.getFinish();
    if(start == end || start == end.rc())
        return {nullptr, nullptr};
    if(end.inDeg() != 1 || end.outDeg() != 2)
        return {nullptr, nullptr};
    Edge &back_edge = (end.front().getFinish() == start) ? end.front() : end.back();
    if(forward_edge.truncSize() > 30000 || back_edge.truncSize() > 50000)
        return {nullptr, nullptr};
    if(back_edge.getStart() != end || back_edge.getFinish() != start)
        return {nullptr, nullptr};
    Edge &out = end.front() == back_edge ? end.back() : end.front();
    Edge &in = start.rc().front().rc() == back_edge ? start.rc().back().rc() : start.rc().front().rc();
    if(component.contains(out.getFinish()) || component.contains(in.getStart()))
        return {nullptr, nullptr};
    return {&forward_edge, &back_edge};
}

RecordStorage ResolveLoops(logging::Logger &logger, size_t threads, SparseDBG &dbg, RecordStorage &reads_storage,
                           const AbstractUniquenessStorage &more_unique) {
    RecordStorage res(dbg, 0, 10000000000ull, threads, reads_storage.getLogger(), false, reads_storage.log_changes);
    for(const Component &comp : UniqueSplitter(more_unique).splitGraph(dbg)) {
        std::pair<Edge *, Edge *> check = CheckLoopComponent(comp);
        if(check.first == nullptr)
            continue;
        Edge &forward_edge = *check.first;
        Edge &back_edge = *check.second;
        Vertex &start = forward_edge.getStart();
        Vertex &end = forward_edge.getFinish();
        Edge &out = end.front() == back_edge ? end.back() : end.front();
        Edge &in = start.rc().front().rc() == back_edge ? start.rc().back().rc() : start.rc().front().rc();
        std::pair<double, double> tmp = minmaxCov(comp, reads_storage, more_unique.asFunction());
        double min_cov = tmp.first;
        double max_cov = tmp.second;
        double med_cov = (max_cov + min_cov) / 2;
        double dev = std::max<double>(4, max_cov - min_cov) / med_cov / 2;
        size_t vote1 = floor(forward_edge.getCoverage() / med_cov + 0.5);
        size_t vote2 = floor(back_edge.getCoverage() / med_cov + 0.5);
        if(vote1 * dev * 2 > med_cov || vote1 != vote2 + 1)
            continue;
        dbg::GraphPath longest = reads_storage.getRecord(in.getStart()).
                getFullUniqueExtension(in.truncSeq().Subseq(0, 1), 1, 0).unpack();
        size_t pos = longest.find(out);
        if(pos != size_t(-1)) {
            VERIFY(pos % 2 == 0 && pos >= 2);
            if(pos / 2 != vote1) {
                logger.trace() << "Coverage contradicts bridging read. Skipping loop " << forward_edge.getInnerId() << " "
                               << back_edge.getInnerId() << " with size " << forward_edge.truncSize() + back_edge.truncSize()
                        << " and multiplicity " << vote2 << " vs "  << pos / 2 - 1 << std::endl;
            }
            continue;
        }
        Sequence bad = (back_edge.truncSeq().Subseq(0, 1) + forward_edge.truncSeq().Subseq(0, 1)) * (vote2 + 1);
        if(reads_storage.getRecord(end).countStartsWith(bad) > 0) {
            logger.trace() << "Coverage contradicts circling read. Skipping loop " << forward_edge.getInnerId() << " "
                           << back_edge.getInnerId() << " with size " << forward_edge.truncSize() + back_edge.truncSize()
                        << " and multiplicity " << vote2 << std::endl;
        }
        dbg::GraphPath alignment;
        alignment += Segment<Edge>(in, in.truncSize() - std::min<size_t>(in.truncSize(), 1000), in.truncSize());
        alignment += forward_edge;
        for(size_t i = 0; i < vote2; i++) {
            alignment += back_edge;
            alignment += forward_edge;
        }
        alignment += Segment<Edge>(out, 0, std::min<size_t>(out.truncSize(), 1000));
        res.addRead(AlignedRead(back_edge.getInnerId().str() + "_" + itos(vote2), alignment));
        logger.trace() << "Resolved loop " << forward_edge.getInnerId() << " " << back_edge.getInnerId() <<
                       " with size " << forward_edge.truncSize() + back_edge.truncSize() << " and multiplicity " << vote2 << std::endl;
    }
    return std::move(res);
}
