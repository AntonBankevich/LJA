#include <common/simple_computation.hpp>
#include "diploidy_analysis.hpp"
#include "multiplicity_estimation.hpp"
using namespace dbg;
size_t BoundRecord::inf = 1000000000000ul;

MappedNetwork::MappedNetwork(const Component &component, const std::function<bool(const dbg::Edge &)> &unique,
                             double rel_coverage, double unique_coverage, double double_coverage) {
    dbg::SparseDBG &graphr = component.graph();
    for(dbg::Vertex &v : component.vertices()) {
        vertex_mapping[&v] = addVertex();
    }
    for(dbg::Vertex &v : component.vertices()) {
        for(dbg::Edge &edge : v) {
            if (!unique(edge)) {
                size_t min_flow = (!edge.is_reliable || edge.getCoverage() < rel_coverage) ? 0 : 1;
                size_t max_flow = 1000000000;
                if(edge.getCoverage() < double_coverage)
                    max_flow = 2;
                if(edge.getCoverage() < unique_coverage)
                    max_flow = 1;
                int eid = addEdge(vertex_mapping[&v], vertex_mapping[edge.end()], min_flow, max_flow);
                edge_mapping[eid] = &edge;
            } else {
                addSink(vertex_mapping[&v], 1);
                addSource(vertex_mapping[&v.rc()], 1);
            }
        }
    }
}

std::vector<dbg::Edge *> MappedNetwork::getUnique(logging::Logger &logger) {
    std::vector<dbg::Edge*> res;
    std::unordered_map<int, size_t> multiplicities = findFixedMultiplicities();
    for (auto &rec : multiplicities) {
//        logger << "Edge " << edge_mapping[rec.first]->start()->hash() << edge_mapping[rec.first]->start()->isCanonical()
//               << "ACGT"[edge_mapping[rec.first]->seq[0]]
//               << " has fixed multiplicity " << rec.second << std::endl;
        if(rec.second == 1)
            res.emplace_back(edge_mapping[rec.first]);
    }
    return std::move(res);
}

size_t MappedNetwork::addTipSinks() {
    size_t res = 0;
    for(Edge &edge : edges) {
        if(edge.min_flow == 0)
            continue;
        Vertex end = vertices[edge.end];
        if(end.out.empty()) {
            addSink(end.id, 1);
            res += 1;
        }
        Vertex start = vertices[edge.start];
        if(start.out.empty()) {
            addSource(start.id, 1);
            res += 1;
        }
    }
    return res;
}

MultiplicityBoundsEstimator::MultiplicityBoundsEstimator(SparseDBG &dbg,
                                                         const AbstractUniquenessStorage &uniquenessStorage) : dbg(dbg){
    for(auto & it : dbg) {
        for(auto v_it : {&it.second, &it.second.rc()}) {
            for(Edge &edge : *v_it) {
                if (uniquenessStorage.isUnique(edge)) {
                    bounds.updateBounds(edge, 1, 1);
                }
            }
        }
    }
}

bool MultiplicityBoundsEstimator::updateComponent(logging::Logger &logger, const Component &component,
                                                  const AbstractUniquenessStorage &uniquenessStorage,
                                                  double rel_coverage, double unique_coverage) {
    std::unordered_set<const dbg::Edge *> unique_in_component;
    std::function<bool(const dbg::Edge &)> is_unique =
            [&uniquenessStorage, unique_coverage, rel_coverage](const dbg::Edge &edge) {
                return uniquenessStorage.isUnique(edge) || (edge.getCoverage() >= rel_coverage && edge.getCoverage() < unique_coverage);
            };
    MappedNetwork net(component, is_unique, rel_coverage);
    bool res = net.fillNetwork();
    if(!res) {
        logger << "Initial flow search failed. Adding tip sinks." << std::endl;
        size_t tips = net.addTipSinks();
        if(tips > 0)
            res = net.fillNetwork();
    }
    if(res) {
        logger << "Found multiplicity bounds in component" << std::endl;
        for(auto rec : net.findBounds()) {
            this->bounds.updateBounds(*net.edge_mapping[rec.first], rec.second.first, rec.second.second);
        }
        return true;
    }
    logger << "Flow search failed. Multiplicity bounds were not updated." << std::endl;
    return false;
}

void MultiplicityBoundsEstimator::update(logging::Logger &logger, double rel_coverage,
                                         const std::experimental::filesystem::path &dir) {
    ensure_dir_existance(dir);
    size_t cnt = 0;
    for(const Component &component : UniqueSplitter(bounds).splitGraph(dbg)) {
        if(component.size() <= 2)
            continue;
        cnt += 1;
        std::ofstream os1;
        os1.open(dir / (std::to_string(cnt) + "_before.dot"));
        printDot(os1, component, bounds.labeler(), bounds.colorer());
        os1.close();
        updateComponent(logger, component, bounds, rel_coverage);
        std::experimental::filesystem::path out_file = dir / (std::to_string(cnt) + ".dot");
        logger << "Printing component to " << out_file << std::endl;
        std::ofstream os;
        os.open(out_file);
        printDot(os, component, bounds.labeler(), bounds.colorer());
        os.close();
    }
}

void UniqueClassificator::markPseudoHets() const {
    for(Edge &edge : dbg.edges()) {
        if(!isUnique(edge) || edge.end()->outDeg() != 2 || edge.end()->inDeg() != 1)
            continue;
        Vertex &start = *edge.end();
        Edge &correct = start[0].getCoverage() > start[1].getCoverage() ? start[0] : start[1];
        Edge &incorrect = start[0].getCoverage() <= start[1].getCoverage() ? start[0] : start[1];
        incorrect.is_reliable = false;
        incorrect.rc().is_reliable = false;
        GraphAlignment cor_ext = reads_storage.getRecord(start).
                getFullUniqueExtension(correct.seq.Subseq(0, 1), 1, 0).getAlignment();
        GraphAlignment incor_ext = reads_storage.getRecord(start).
                getFullUniqueExtension(incorrect.seq.Subseq(0, 1), 1, 0).getAlignment();
        for(Segment<Edge> &seg : incor_ext) {
            bool found= false;
            for(Segment<Edge> &seg1 : cor_ext) {
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
        SetUniquenessStorage duninque = BulgePathAnalyser(dbg, unique_len).uniqueEdges();
        for (Edge &edge : dbg.edges()) {
            if (duninque.isUnique(edge)) {
                updateBounds(edge, 1, 1);
                cnt++;
            }
        }
    } else {
        for (Edge &edge : dbg.edges()) {
            if(edge.size() > unique_len || (edge.start()->inDeg() == 0 && edge.size() > unique_len / 3)) {
                updateBounds(edge, 1, 1);
                cnt++;
            }
        }
    }
    logger.info() << "Marked " << cnt << " long edges as unique" << std::endl;
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
        logger.trace() << "Component parameters: size=" << component.size() << " border=" << component.countBorderEdges() <<
                      " tips=" << component.countTips() <<
                      " subcomponents=" << component.realCC() << " acyclic=" << component.isAcyclic() <<std::endl;
        if(component.size() > 2 && component.countBorderEdges() == 2 &&component.countTips() == 0 &&
           component.realCC() == 2 && component.isAcyclic()) {
            processSimpleComponent(logger, component);
        }
        cnt += processComponent(logger, component);
        if(debug)
            logger.trace() << "Printing component to " << (dir / (std::to_string(component_cnt) + ".dot")) << std::endl;
        if(debug)
            printDot(dir / (std::to_string(component_cnt) + ".dot"), component,
                     this->labeler() + reads_storage.labeler(), this->colorer());
    }
    logger.info() << "Finished unique edges search. Found " << cnt << " unique edges" << std::endl;
    logger.info() << "Analysing repeats of multiplicity 2 and looking for additional unique edges" << std::endl;
    std::function<bool(const dbg::Edge &)> mult2 = [this](const dbg::Edge &edge) {
        return MultiplicityBounds::lowerBound(edge) == 2 && MultiplicityBounds::lowerBound(edge) == 2;
    };
    split = ConditionSplitter(mult2).splitGraph(dbg);
    cnt = 0;
    for(Component &component : split) {
        if(component.size() != 8 || !component.isAcyclic() || component.realCC() != 2 || component.borderVertices().size() != 2) {
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
    Vertex &start = border[0]->rc();
    Vertex &end = *border[1];
    if(start.outDeg() !=2 || end.inDeg() != 2)
        return false;
    std::vector<Edge *> bad_candidates;
    if(start[0].end() == &end || start[1].end() == &end) {
        size_t ind = start[0].end() == &end ? 0 : 1;
        Edge &e1 = start[ind];
        Edge &e21 = start[1 - ind];
        if(e21.end()->outDeg() != 2)
            return false;
        Edge &e221 = e21.end()->operator[](0);
        Edge &e222 = e21.end()->operator[](1);
        if(e221.end() != e222.end())
            return false;
        if(e221.end()->outDeg() != 1)
            return false;
        Edge &e23 = e221.end()->operator[](0);
        if(e23.end() != &end)
            return false;
        if(!component.contains(*e21.end()) || !component.contains(*e23.start()))
            return false;
        bad_candidates = {&e221, &e222};
    } else {
        size_t ind1 = start[0].size() > start[1].size() ? 0 : 1;
        Edge &e11 = start[ind1];
        Edge &e12 = start[1 - ind1];
        if(end.rc().inDeg() != 2)
            return false;
        size_t ind2 = end.rc()[0].size() > end.rc()[1].size() ? 1 : 0;
        Edge &e31 = start[ind1];
        Edge &e32 = start[1 - ind1];
        if(e12.end()->outDeg() != 2 || e31.start()->inDeg() != 2 || e11.end() != e31.start() || e12.end() != e32.start())
            return false;
        Edge &e2 = e12.end()->operator[](0).end() == e31.start() ? e12.end()->operator[](0) : e12.end()->operator[](1);
        if(!component.contains(*e2.start()) || ! component.contains(*e2.end()))
            return false;
        bad_candidates = {&e11, &e2, &e32};
    }
    std::function<double(Edge* const &)> f = [](Edge * const &edge){return edge->getCoverage();};
    Edge &bad = *bad_candidates[MinIndex<Edge *, double>(bad_candidates, f)];
    for(Edge &edge : component.edgesInner()) {
        if(edge != bad && edge != bad.rc())
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
    logger.trace() << "Attempting to use coverage for multiplicity estimation with coverage threshold " << threshold << std::endl;
    logger.trace() << "Component: ";
    for(Vertex &vertex : subcomponent.verticesUnique()) {
        logger << " " << vertex.getShortId();
    }
    logger << std::endl;
//    rel_coverage = 0;
    MappedNetwork net2(subcomponent, is_unique, rel_coverage, threshold);
    bool res = net2.fillNetwork();
    if (res) {
        logger.trace() << "Succeeded to use coverage for multiplicity estimation" << std::endl;
        for(auto rec : net2.findBounds()) {
            if(!MultiplicityBounds::isUnique(*net2.edge_mapping[rec.first]) && rec.second.first == 1 && rec.second.second == 1) {
                ucnt++;
            }
            updateBounds(*net2.edge_mapping[rec.first], rec.second.first, rec.second.second);
        }
    } else {
        logger.trace() << "Failed to use coverage for multiplicity estimation" << std::endl;
        if(rel_coverage == 0) {
            logger.trace() << "Adjusted reliable edge threshold from " << rel_coverage << " to " <<
                           adjusted_rel_coverage << std::endl;
            rel_coverage = adjusted_rel_coverage;
            MappedNetwork net3(subcomponent, is_unique, rel_coverage, threshold);
            res = net3.fillNetwork();
            if (res) {
                logger.trace() << "Succeeded to use coverage for multiplicity estimation" << std::endl;
                for(auto rec : net3.findBounds()) {
                    if(!MultiplicityBounds::isUnique(*net3.edge_mapping[rec.first]) && rec.second.first == 1 && rec.second.second == 1) {
                        ucnt++;
                    }
                    updateBounds(*net3.edge_mapping[rec.first], rec.second.first, rec.second.second);
                }
            } else {
                logger.trace() << "Failed to use coverage for multiplicity estimation" << std::endl;
            }
        }
    }
    MappedNetwork net4(subcomponent, is_unique, rel_coverage, threshold, double_threshold);
    res = net4.fillNetwork();
    if(res) {
        for(auto rec : net4.findBounds()) {
            if(!MultiplicityBounds::isUnique(*net4.edge_mapping[rec.first]) && rec.second.first == 1 && rec.second.second == 1) {
                ucnt++;
            }
            updateBounds(*net4.edge_mapping[rec.first], rec.second.first, rec.second.second);
        }
    } else {
        double_threshold = 0;
    }
    logger.trace() << "Succeeded to use coverage for multiplicity estimation with rel_coverage = " << rel_coverage
                   << " threshold = " << threshold << " double_threshold = " << double_threshold << std::endl;
    return ucnt;
}

std::pair<double, double> minmaxCov(const Component &subcomponent, const RecordStorage &reads_storage,
                               const std::function<bool(const dbg::Edge &)> &is_unique) {
    double max_cov= 0;
    double min_cov= 100000;
    for(Edge &edge : subcomponent.edges()) {
        if(is_unique(edge) && subcomponent.contains(*edge.end())) {
            const VertexRecord & record = reads_storage.getRecord(*edge.start());
            std::string s = edge.seq.Subseq(0, 1).str();
            size_t cnt = record.countStartsWith(Sequence(s + "A")) +
                         record.countStartsWith(Sequence(s + "C")) +
                         record.countStartsWith(Sequence(s + "G")) +
                         record.countStartsWith(Sequence(s + "T"));
            if(edge.size() < 20000) {
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
        if(!component.contains(*edge.end()))
            return edge.rc();
    }
    VERIFY(false);
}

std::vector<Vertex *> topSort(const Component &component) {
    std::vector<Vertex *> stack;
    std::vector<Vertex *> res;
    std::unordered_set<Vertex *> visited;
    stack.emplace_back(getStart(component).end());
    while(!stack.empty()) {
        Vertex *cur = stack.back();
        if(visited.find(cur) != visited.end()) {
            stack.pop_back();
            continue;
        }
        bool ok = true;
        for(Edge &e : *cur) {
            if(component.contains(*e.end()) && visited.find(e.end()) == visited.end())
                ok = false;
        }
        if(ok) {
            visited.emplace(cur);
            stack.pop_back();
            res.emplace_back(&cur->rc());
        } else {
            for(Edge &e : *cur) {
                VERIFY(component.contains(*e.end()));
                if(visited.find(e.end()) == visited.end())
                    stack.emplace_back(e.end());
            }
        }
    }
    return std::move(res);
}
void UniqueClassificator::processSimpleComponent(logging::Logger &logger, const Component &component) const {
    logger.trace() << "Collapsing acyclic component" << std::endl;
    std::vector<Vertex *> order = topSort(component);
    if(order.size() != component.size()) {
        logger.trace() << "Failed to collapse acyclic component" << std::endl;
        return;
    }
    VERIFY(order.front()->inDeg() == 1);
    VERIFY(!component.contains(*order.front()->rc()[0].end()));
    VERIFY(order.back()->outDeg() == 1);
    VERIFY(!component.contains(*order.back()->operator[](0).end()));
    std::unordered_map<Vertex *, std::pair<Edge *, size_t>> prev;
    for(Vertex *cur : order) {
        Edge *pedge = nullptr;
        size_t val = 0;
        for(Edge &edge : cur->rc()) {
            if(!component.contains(*edge.end())) {
                break;
            }
            VERIFY(prev.find(&edge.end()->rc()) != prev.end());
            size_t score = edge.intCov() + prev[&edge.end()->rc()].second;
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
    while(prev[end].first != nullptr) {
        prev[end].first->is_reliable = true;
        prev[end].first->rc().is_reliable = true;
        end = prev[end].first->start();
    }
}

size_t UniqueClassificator::processComponent(logging::Logger &logger, const Component &component) {
    size_t ucnt = 0;
    logger.trace() << "Component: ";
    for(Vertex &vertex : component.verticesUnique()) {
        logger << " " << vertex.getShortId();
    }
    logger << std::endl;
    double rel_coverage = 0;
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
        rel_coverage = 15;
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
    if(component.size() != 2)
        return {nullptr, nullptr};
    Vertex *vit = &*component.vertices().begin();
    if(vit->inDeg() != 2)
        vit = &vit->rc();
    if(vit->inDeg() != 2 || vit->outDeg() != 1)
        return {nullptr, nullptr};
    Vertex &start = *vit;
    Edge &forward_edge = start[0];
    Vertex &end = *forward_edge.end();
    if(start == end || start == end.rc())
        return {nullptr, nullptr};
    if(end.inDeg() != 1 || end.outDeg() != 2)
        return {nullptr, nullptr};
    Edge &back_edge = (end[0].end() == &start) ? end[0] : end[1];
    if(forward_edge.size() > 30000 || back_edge.size() > 50000)
        return {nullptr, nullptr};
    if(back_edge.start() != &end || back_edge.end() != &start)
        return {nullptr, nullptr};
    Edge &out = end[0] == back_edge ? end[1] : end[0];
    Edge &in = start.rc()[0].rc() == back_edge ? start.rc()[1].rc() : start.rc()[0].rc();
    if(component.contains(*out.end()) || component.contains(*in.start()))
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
        Vertex &start = *forward_edge.start();
        Vertex &end = *forward_edge.end();
        Edge &out = end[0] == back_edge ? end[1] : end[0];
        Edge &in = start.rc()[0].rc() == back_edge ? start.rc()[1].rc() : start.rc()[0].rc();
        std::pair<double, double> tmp = minmaxCov(comp, reads_storage, more_unique.asFunction());
        double min_cov = tmp.first;
        double max_cov = tmp.second;
        double med_cov = (max_cov + min_cov) / 2;
        double dev = std::max<double>(4, max_cov - min_cov) / med_cov / 2;
        size_t vote1 = floor(forward_edge.getCoverage() / med_cov + 0.5);
        size_t vote2 = floor(back_edge.getCoverage() / med_cov + 0.5);
        if(vote1 * dev * 2 > med_cov || vote1 != vote2 + 1)
            continue;
        GraphAlignment longest = reads_storage.getRecord(*in.start()).
                getFullUniqueExtension(in.seq.Subseq(0, 1), 1, 0).getAlignment();
        size_t pos = longest.find(out);
        if(pos != size_t(-1)) {
            VERIFY(pos % 2 == 0 && pos >= 2);
            if(pos / 2 != vote1) {
                logger.trace() << "Coverage contradicts bridging read. Skipping loop " << forward_edge.getId() << " "
                        << back_edge.getId() << " with size " << forward_edge.size() + back_edge.size()
                        << " and multiplicity " << vote2 << " vs "  << pos / 2 - 1 << std::endl;
            }
            continue;
        }
        Sequence bad = (back_edge.seq.Subseq(0, 1) + forward_edge.seq.Subseq(0, 1)) * (vote2 + 1);
        if(reads_storage.getRecord(end).countStartsWith(bad) > 0) {
            logger.trace() << "Coverage contradicts circling read. Skipping loop " << forward_edge.getId() << " "
                        << back_edge.getId() << " with size " << forward_edge.size() + back_edge.size()
                        << " and multiplicity " << vote2 << std::endl;
        }
        GraphAlignment alignment;
        alignment += Segment<Edge>(in, in.size() - std::min<size_t>(in.size(), 1000), in.size());
        alignment += forward_edge;
        for(size_t i = 0; i < vote2; i++) {
            alignment += back_edge;
            alignment += forward_edge;
        }
        alignment += Segment<Edge>(out, 0, std::min<size_t>(out.size(), 1000));
        res.addRead(AlignedRead(back_edge.getId() + "_" + itos(vote2), alignment));
        logger.trace() << "Resolved loop " << forward_edge.getId() << " " << back_edge.getId() <<
                " with size " << forward_edge.size() + back_edge.size() << " and multiplicity " << vote2 << std::endl;
    }
    return std::move(res);
}
