#include "graph.hh"
#include "solve_mcs.hh"

#include <algorithm>
#include <iostream>
#include <set>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

std::atomic<bool> abort_due_to_timeout;

namespace {
    struct Bidomain
    {
        // start indices of left and right sets
        int l;
        int r;

        int left_len;
        int right_len;
    };

    struct DomainStore
    {
        vector<int> left;
        vector<int> right;
        vector<Bidomain> domains;
    };

    class MCS
    {
        const Graph & linegraph0;
        const Graph & linegraph1;
        const std::vector<Edge> edges0;
        const std::vector<Edge> edges1;
        const Graph & originalgraph0;
        const Graph & originalgraph1;
        const Params params;
        McsStats stats;
        vector<Assignment> incumbent;

        auto print_state(const vector<Assignment>& current, const DomainStore &domain_store) -> void
        {
            cout << "Nodes: " << stats.nodes << endl;
            cout << "Length of current assignment: " << current.size() << endl;
            cout << "Current assignment:";
            for (unsigned int i=0; i<current.size(); i++) {
                cout << "  (" << current[i].v << " -> " << current[i].w << ")";
            }
            cout << endl;
            for (unsigned int i=0; i<domain_store.domains.size(); i++) {
                struct Bidomain bd = domain_store.domains[i];
                cout << "Left  ";
                for (int j=0; j<bd.left_len; j++)
                    cout << domain_store.left[bd.l + j] << " ";
                cout << endl;
                cout << "Right  ";
                for (int j=0; j<bd.right_len; j++)
                    cout << domain_store.right[bd.r + j] << " ";
                cout << endl;
            }
            cout << "\n" << endl;
        }

        auto calc_bound(const vector<Bidomain>& domains) -> int
        {
            int bound = 0;
            for (const Bidomain &bd : domains)
                bound += std::min(bd.left_len, bd.right_len);
            return bound;
        }

        // precondition: len > 0
        auto find_min_value(const vector<int>& arr, int start_idx, int len) -> int
        {
            return *std::min_element(arr.begin() + start_idx, arr.begin() + start_idx + len);
        }

        auto bd_is_adjacent(const Bidomain & bd, const std::vector<int> & vtx_counts0,
                const DomainStore & domain_store) -> bool
        {
            int v = domain_store.left[bd.l];
            Edge e = edges0[v];
            return vtx_counts0[e.v] || vtx_counts0[e.w];
        }

        auto select_bidomain(const DomainStore & domain_store, int current_matching_size,
                const std::vector<int> & vtx_counts0) -> int
        {
            // Select the bidomain with the smallest max(leftsize, rightsize), breaking
            // ties on the smallest vertex index in the left set
            std::pair<int, int> best_score { std::numeric_limits<int>::max(), std::numeric_limits<int>::max() };
            int best = -1;
            for (unsigned int i=0; i<domain_store.domains.size(); i++) {
                const Bidomain &bd = domain_store.domains[i];
                if (params.connected && current_matching_size > 0 &&
                        !bd_is_adjacent(bd, vtx_counts0, domain_store)) {
                    continue;
                }
                int len = params.heuristic == Heuristic::min_max ?
                        std::max(bd.left_len, bd.right_len) :
                        bd.left_len * bd.right_len;
                int tie_breaker = find_min_value(domain_store.left, bd.l, bd.left_len);
                auto score = std::make_pair( len, tie_breaker );
                if (score < best_score) {
                    best_score = score;
                    best = i;
                }
            }
            return best;
        }

        // Returns length of left half of array
        auto partition(vector<int>& all_vv, int start, int len, const vector<unsigned int> & adjrow) -> int
        {
            auto it = std::partition(
                    all_vv.begin() + start,
                    all_vv.begin() + start + len,
                    [&](const int elem){ return 0 != adjrow[elem]; });
            return it - (all_vv.begin() + start);
        }

        auto refined_domains(const DomainStore & ds, int v, int w) -> DomainStore
        {
            DomainStore new_ds { ds.left, ds.right, {} };
            new_ds.domains.reserve(ds.domains.size());
            for (const Bidomain &old_bd : ds.domains) {
                int l = old_bd.l;
                int r = old_bd.r;
                // After these two partitions, left_len and right_len are the lengths of the
                // arrays of vertices with edges from v or w
                int left_len = partition(new_ds.left, l, old_bd.left_len, linegraph0.adjmat[v]);
                int right_len = partition(new_ds.right, r, old_bd.right_len, linegraph1.adjmat[w]);
                int left_len_noedge = old_bd.left_len - left_len;
                int right_len_noedge = old_bd.right_len - right_len;
                if (left_len_noedge && right_len_noedge)
                    new_ds.domains.push_back({l+left_len, r+right_len, left_len_noedge, right_len_noedge});
                if (left_len && right_len)
                    new_ds.domains.push_back({l, r, left_len, right_len});
            }
            return new_ds;
        }

        auto remove_vtx_from_left_domain(vector<int> & left, Bidomain & bd, int v) -> void
        {
            int i = 0;
            while(left[bd.l + i] != v) i++;
            std::swap(left[bd.l+i], left[bd.l+bd.left_len-1]);
            bd.left_len--;
        }

        auto remove_bidomain(vector<Bidomain>& domains, int idx) -> void
        {
            domains[idx] = domains[domains.size()-1];
            domains.pop_back();
        }

        auto update_incumbent(vector<Assignment> & current) -> void
        {
            incumbent = current;

            stats.time_of_last_incumbent_update = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::steady_clock::now() - params.start_time).count();

            if (!params.quiet) cout << "New incumbent " << incumbent.size() << " at time " <<
                    stats.time_of_last_incumbent_update << " ms" << endl;
        }

        enum class Search
        {
            Aborted,
            Done,
            ReachedGoal
        };

        auto update_vtx_counts(
                std::vector<int> & vtx_counts0,
                std::vector<int> & vtx_counts1,
                int v,
                int w) -> bool
        {
            Edge e0 = edges0[v];
            // new_count0 is the number of new non-zeroes in vtx_counts0
            int new_count0 = ((0==vtx_counts0[e0.v]) + (0==vtx_counts0[e0.w]));
            ++vtx_counts0[e0.v];
            ++vtx_counts0[e0.w];
            Edge e1 = edges1[w];
            // new_count1 is the number of new non-zeroes in vtx_counts1
            int new_count1 = ((0==vtx_counts1[e1.v]) + (0==vtx_counts1[e1.w]));
            ++vtx_counts1[e1.v];
            ++vtx_counts1[e1.w];
            return new_count0 == new_count1;
        }

        auto downdate_vtx_counts(
                std::vector<int> & vtx_counts0,
                std::vector<int> & vtx_counts1,
                int v,
                int w) -> void
        {
            Edge e0 = edges0[v];
            --vtx_counts0[e0.v];
            --vtx_counts0[e0.w];
            Edge e1 = edges1[w];
            --vtx_counts1[e1.v];
            --vtx_counts1[e1.w];
        }

        auto search(
                vector<Assignment> & current,
                DomainStore & domain_store,
                unsigned int matching_size_goal,
                std::vector<int> & vtx_counts0,
                std::vector<int> & vtx_counts1) -> Search
        {
            if (abort_due_to_timeout)
                return Search::Aborted;

            ++stats.nodes;

            if (params.verbose)
                print_state(current, domain_store);

            if (current.size() > incumbent.size()) {
                update_incumbent(current);
                if (incumbent.size() == matching_size_goal)
                    return Search::ReachedGoal;
            }

            unsigned int bound = current.size() + calc_bound(domain_store.domains);
            if (bound <= incumbent.size() || (params.mcsplit_down && bound < matching_size_goal))
                return Search::Done;

            int bd_idx = select_bidomain(domain_store, current.size(), vtx_counts0);
            if (bd_idx == -1)  // this might occur if params.connected is true
                return Search::Done;

            Bidomain &bd = domain_store.domains[bd_idx];

            int v = find_min_value(domain_store.left, bd.l, bd.left_len);

            // Try assigning v to each vertex w in the colour class beginning at bd.r, in turn
            auto right_label_class_begin = domain_store.right.begin() + bd.r;
            auto right_label_class_end = right_label_class_begin + bd.right_len;
            int & right_label_class_last_elem = *std::prev(right_label_class_end);
            std::sort(right_label_class_begin, right_label_class_end);

            remove_vtx_from_left_domain(domain_store.left, bd, v);
            bd.right_len--;

            for (auto it=right_label_class_begin; it!=right_label_class_end; ++it) {
                int w = *it;

                Search search_result = Search::Done;
                if (update_vtx_counts(vtx_counts0, vtx_counts1, v, w)) {
                    std::swap(*it, right_label_class_last_elem); // swap w to the end of its label class
                    current.push_back({v, w});
                    auto new_domain_store = refined_domains(domain_store, v, w);
                    search_result = search(current, new_domain_store, matching_size_goal,
                            vtx_counts0, vtx_counts1);
                    current.pop_back();
                    std::swap(*it, right_label_class_last_elem); // swap w back to its correct place in the sorted label class
                }
                downdate_vtx_counts(vtx_counts0, vtx_counts1, v, w);

                switch (search_result)
                {
                case Search::Aborted:     return Search::Aborted;
                case Search::ReachedGoal: return Search::ReachedGoal;
                case Search::Done:        break;
                }
            }

            // Try assigning v to \bot (i.e. not using it in our mapping)
            bd.right_len++;
            if (bd.left_len == 0)
                remove_bidomain(domain_store.domains, bd_idx);
            return search(current, domain_store, matching_size_goal, vtx_counts0, vtx_counts1);
        }

    public:
        MCS(Graph & linegraph0, Graph & linegraph1, const std::vector<Edge> & edges0,
                const std::vector<Edge> & edges1,
                const Graph & originalgraph0, const Graph & originalgraph1,
                Params params)
            : linegraph0(linegraph0), linegraph1(linegraph1), edges0(edges0),
              edges1(edges1), originalgraph0(originalgraph0),
              originalgraph1(originalgraph1), params(params)
        { }

        auto run() -> std::pair<vector<Assignment>, McsStats>
        {
            DomainStore domain_store;

            for (int i=0; i<linegraph0.n; i++)
                domain_store.left.push_back(i);
            for (int i=0; i<linegraph1.n; i++)
                domain_store.right.push_back(i);

            domain_store.domains.push_back({0, 0, linegraph0.n, linegraph1.n});

            vector<Assignment> current;

            // vtx_counts0 and vtx_counts1 have the following meaning, for the original pattern
            // graph and original target graph respectively.  The vtx_counts array has n
            // counters.  The i^th counter keeps track of the number of edges adjacent to
            // vertex i of the original graph that are in the current solution.
            //
            // I'm relying on the claim---which I haven't yet properly proven---that we have
            // a K3/claw problem if and only if the two vtx_counts arrays have different numbers
            // of non-zero elements.  TODO: look again at the RASCAL paper, and see how they handle
            // the K3/claw problem.

	    if (params.mcsplit_down) {
		for (unsigned int goal = std::min(linegraph0.n, linegraph1.n) ; goal > 0 ; --goal) {
		    if (incumbent.size() == goal) break;
                    auto domain_store_copy = domain_store;
                    std::vector<int> vtx_counts0(originalgraph0.n);
                    std::vector<int> vtx_counts1(originalgraph1.n);
                    search(current, domain_store_copy, goal, vtx_counts0, vtx_counts1);
		    if (incumbent.size() == goal || abort_due_to_timeout) break;
		    if (!params.quiet) cout << "Upper bound: " << goal-1 << std::endl;
                    cout << stats.nodes << std::endl;
		}
	    } else {
                std::vector<int> vtx_counts0(originalgraph0.n);
                std::vector<int> vtx_counts1(originalgraph1.n);
                search(current, domain_store, std::min(linegraph0.n, linegraph1.n), vtx_counts0, vtx_counts1);
	    }

            return {incumbent, stats};
        }
    };
};

//auto vertices_sorted_by_degree(Graph & linegraph) -> vector<int>
//{
//    auto deg = calculate_degrees(linegraph);
//    vector<int> vv(linegraph.n);
//    std::iota(std::begin(vv), std::end(vv), 0);
//    std::stable_sort(std::begin(vv), std::end(vv), [&](int a, int b) {
//        return deg[a] > deg[b];
//    });
//    return vv;
//}

auto solve_mcs(Graph & linegraph0, Graph & linegraph1,
                const std::vector<Edge> & edges0,
                const std::vector<Edge> & edges1,
                const Graph & originalgraph0,
                const Graph & originalgraph1,
                Params params)
		-> std::pair<vector<Assignment>, McsStats>
{
//    auto vv0 = vertices_sorted_by_degree(linegraph0);
//    auto vv1 = vertices_sorted_by_degree(linegraph1);
//
//    struct Graph linegraph0_sorted = induced_subgraph(linegraph0, vv0);
//    struct Graph linegraph1_sorted = induced_subgraph(linegraph1, vv1);

//    auto solution = MCS(linegraph0_sorted, linegraph1_sorted, params).run();

    auto solution = MCS(linegraph0, linegraph1, edges0, edges1,
            originalgraph0, originalgraph1, params).run();

//    // Convert to indices from original, unsorted graphs
//    for (auto& vtx_pair : solution.first) {
//        vtx_pair.v = vv0[vtx_pair.v];
//        vtx_pair.w = vv1[vtx_pair.w];
//    }

    return solution;
}
