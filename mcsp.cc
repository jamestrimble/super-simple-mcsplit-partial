#include "graph.hh"
#include "solve_mcs.hh"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

#include <argp.h>

using std::cout;
using std::endl;
using std::vector;

static void fail(std::string msg) {
    std::cerr << msg << std::endl;
    exit(1);
}

/*******************************************************************************
                             Command-line arguments
*******************************************************************************/

static char doc[] = "Find a maximum clique in a graph in DIMACS format\vHEURISTIC can be min_max or min_product";
static char args_doc[] = "HEURISTIC FILENAME1 FILENAME2";
static struct argp_option options[] = {
    {"quiet", 'q', 0, 0, "Quiet output"},
    {"verbose", 'v', 0, 0, "Verbose output"},
    {"dimacs", 'd', 0, 0, "Read DIMACS format"},
    {"lad", 'l', 0, 0, "Read LAD format"},
    {"mcsplit-down", 'm', 0, 0, "Use the McSplit-down algorithm rather than branch and bound"},
    {"connected", 'c', 0, 0, "Require the subgraphs to be connected"},
    {"timeout", 't', "timeout", 0, "Specify a timeout (seconds)"},
    { 0 }
};

struct Arguments {
    bool quiet;
    bool verbose;
    bool dimacs;
    bool lad;
    bool mcsplit_down;
    bool connected;
    Heuristic heuristic;
    char *filename1;
    char *filename2;
    int timeout;
    int arg_num;
};

static Arguments arguments;

static error_t parse_opt (int key, char *arg, struct argp_state *state) {
    switch (key) {
        case 'd':
            if (arguments.lad)
                fail("The -d and -l options cannot be used together.\n");
            arguments.dimacs = true;
            break;
        case 'l':
            if (arguments.dimacs)
                fail("The -d and -l options cannot be used together.\n");
            arguments.lad = true;
            break;
        case 'q':
            arguments.quiet = true;
            break;
        case 'v':
            arguments.verbose = true;
            break;
        case 'm':
            arguments.mcsplit_down = true;
            break;
        case 'c':
            arguments.connected = true;
            break;
        case 't':
            arguments.timeout = std::stoi(arg);
            break;
        case ARGP_KEY_ARG:
            if (arguments.arg_num == 0) {
                if (std::string(arg) == "min_max")
                    arguments.heuristic = Heuristic::min_max;
                else if (std::string(arg) == "min_product")
                    arguments.heuristic = Heuristic::min_product;
                else
                    fail("Unknown heuristic (try min_max or min_product)");
            } else if (arguments.arg_num == 1) {
                arguments.filename1 = arg;
            } else if (arguments.arg_num == 2) {
                arguments.filename2 = arg;
            } else {
                argp_usage(state);
            }
            arguments.arg_num++;
            break;
        case ARGP_KEY_END:
            if (arguments.arg_num == 0)
                argp_usage(state);
            break;
        default: return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

/*******************************************************************************
                                     Stats
*******************************************************************************/

unsigned long long nodes{ 0 };

/*******************************************************************************
*******************************************************************************/

vector<Edge> get_edges(const Graph & g)
{
    vector<Edge> edges;
    for (int i=0; i<g.n; i++)
        for (int j=0; j<i; j++)
            if (g.adjmat[i][j])
                edges.push_back({i, j});

    auto deg = calculate_degrees(g);

    std::sort(edges.begin(), edges.end(), [&deg](auto e0, auto e1){
                return deg[e0.v] + deg[e0.w] > deg[e1.v] + deg[e1.w];
            });

    return edges;
}

Graph make_linegraph(const Graph & g, const vector<Edge> & edges)
{
    Graph linegraph(edges.size());
    for (unsigned i=0; i<edges.size(); i++) {
        Edge ei = edges[i];
        for (unsigned j=0; j<edges.size(); j++) {
            Edge ej = edges[j];
            if (i != j && (ei.v == ej.v ||
                           ei.v == ej.w ||
                           ei.w == ej.v ||
                           ei.w == ej.w)) {
                linegraph.adjmat[i][j] = 1;
                linegraph.adjmat[j][i] = 1;
            }
        }
    }

    return linegraph;
}

struct EdgeEdgeMapping
{
    Edge edge0;
    Edge edge1;
};

bool check_sol(const Graph & g0, const Graph & g1,
        const vector<EdgeEdgeMapping> & edge_edge_mappings)
{
    std::vector<int> vtx_vtx_mappings(g0.n, -1);

    for (EdgeEdgeMapping eem : edge_edge_mappings) {
        int v0 = eem.edge0.v;
        int w0 = eem.edge0.w;
        int v1 = eem.edge1.v;
        int w1 = eem.edge1.w;
        if (vtx_vtx_mappings[v0] == -1) {
            vtx_vtx_mappings[v0] = v1;
        } else if (vtx_vtx_mappings[v0] != v1) {
            return false;
        }
        if (vtx_vtx_mappings[w0] == -1) {
            vtx_vtx_mappings[w0] = w1;
        } else if (vtx_vtx_mappings[w0] != w1) {
            return false;
        }
        if (!g0.adjmat[v0][w0])
            return false;
        if (!g1.adjmat[v1][w1])
            return false;
    }
    return true;
}

auto get_edge_edge_mappings(
        const std::vector<Assignment> & solution,
        const std::vector<Edge> & edges0,
        const std::vector<Edge> & edges1,
        const Graph & g0) -> vector<EdgeEdgeMapping>
{
    // Keep track of the vertices in the pattern graph that must
    // be mapped to specific vertices in the target graph
    std::vector<int> vtx_vtx_mappings(g0.n, -1);

    std::vector<EdgeEdgeMapping> edge_edge_mappings;
    for (Assignment a : solution) {
        edge_edge_mappings.push_back({edges0[a.v], edges1[a.w]});
    }

    // Look at each pair of edge mappings.  If the pattern-graph
    // edges have an endpoint in common, then that vertex must be mapped
    // to the vertex that is an endpoint of both the target-graph edges.
    for (unsigned i=0; i<edge_edge_mappings.size(); i++) {
        Edge edge0i = edge_edge_mappings[i].edge0;
        Edge edge1i = edge_edge_mappings[i].edge1;
        for (unsigned j=i+1; j<edge_edge_mappings.size(); j++) {
            Edge edge0j = edge_edge_mappings[j].edge0;
            Edge edge1j = edge_edge_mappings[j].edge1;

            // TODO: refactor this.  Maybe factor out a function?
            int common_pattern_graph_vtx = -1;
            if (edge0i.v == edge0j.v || edge0i.v == edge0j.w) {
                common_pattern_graph_vtx = edge0i.v;
            } else if (edge0i.w == edge0j.v || edge0i.w == edge0j.w) {
                common_pattern_graph_vtx = edge0i.w;
            }
            if (common_pattern_graph_vtx != -1) {
                if (edge1i.v == edge1j.v || edge1i.v == edge1j.w) {
                    vtx_vtx_mappings[common_pattern_graph_vtx] = edge1i.v;
                } else {
                    vtx_vtx_mappings[common_pattern_graph_vtx] = edge1i.w;
                }
            }
        }
    }
    // Orient the edges
    for (EdgeEdgeMapping & eem : edge_edge_mappings) {
        int u = vtx_vtx_mappings[eem.edge0.v];
        if (u != -1) {
            if (eem.edge1.v != u) {
                std::swap(eem.edge1.v, eem.edge1.w);
            }
        } else {
            int u = vtx_vtx_mappings[eem.edge0.w];
            if (u != -1 && eem.edge1.w != u) {
                std::swap(eem.edge1.v, eem.edge1.w);
            }
        }
    }
    return edge_edge_mappings;
}

int main(int argc, char** argv)
{
    argp_parse(&argp, argc, argv, 0, 0, 0);

    char format = arguments.dimacs ? 'D' : arguments.lad ? 'L' : 'B';
    struct Graph g0 = readGraph(arguments.filename1, format);
    struct Graph g1 = readGraph(arguments.filename2, format);

    std::thread timeout_thread;
    std::mutex timeout_mutex;
    std::condition_variable timeout_cv;
    abort_due_to_timeout.store(false);
    bool aborted = false;

    if (0 != arguments.timeout) {
        timeout_thread = std::thread([&] {
                auto abort_time = std::chrono::steady_clock::now() + std::chrono::seconds(arguments.timeout);
                {
                    /* Sleep until either we've reached the time limit,
                     * or we've finished all the work. */
                    std::unique_lock<std::mutex> guard(timeout_mutex);
                    while (! abort_due_to_timeout.load()) {
                        if (std::cv_status::timeout == timeout_cv.wait_until(guard, abort_time)) {
                            /* We've woken up, and it's due to a timeout. */
                            aborted = true;
                            break;
                        }
                    }
                }
                abort_due_to_timeout.store(true);
                });
    }

    auto start = std::chrono::steady_clock::now();

    Params params = {
        arguments.quiet,
        arguments.verbose,
        arguments.mcsplit_down,
        arguments.connected,
        arguments.heuristic,
        start
    };

    //////////////////////////////////////////
    // TODO: move this stuff into the solver file
    auto edges0 = get_edges(g0);
    auto edges1 = get_edges(g1);

    auto linegraph0 = make_linegraph(g0, edges0);
    auto linegraph1 = make_linegraph(g1, edges1);
    //////////////////////////////////////////

    auto result = solve_mcs(linegraph0, linegraph1, edges0, edges1,
            g0, g1, params);
    auto solution = result.first;
    auto stats = result.second;

    auto stop = std::chrono::steady_clock::now();
    auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    /* Clean up the timeout thread */
    if (timeout_thread.joinable()) {
        {
            std::unique_lock<std::mutex> guard(timeout_mutex);
            abort_due_to_timeout.store(true);
            timeout_cv.notify_all();
        }
        timeout_thread.join();
    }

    auto edge_edge_mappings = get_edge_edge_mappings(solution, edges0, edges1, g0);

    if (!check_sol(g0, g1, edge_edge_mappings))
        fail("*** Error: Invalid solution\n");

    for (auto m : edge_edge_mappings) {
        cout << "(" <<
        m.edge0.v << "," << m.edge0.w << " -> " <<
        m.edge1.v << "," << m.edge1.w << ") ";
    }
    cout << std::endl;

    cout << "Solution size " << solution.size() << std::endl;
//    for (int i=0; i<g0.n; i++)
//        for (unsigned int j=0; j<solution.size(); j++)
//            if (solution[j].v == i)
//                cout << "(" << solution[j].v << " -> " << solution[j].w << ") ";
    cout << std::endl;

    cout << "Nodes:                              " << stats.nodes << endl;
    cout << "Time of last incumbent update (ms): " << stats.time_of_last_incumbent_update << endl;
    cout << "CPU time (ms):                      " << time_elapsed << endl;
    if (aborted)
        cout << "TIMEOUT" << endl;

    cout << "Summary:" << std::endl;
    cout << (aborted ? "TIMEOUT" : "COMPLETED") << " " <<
            stats.nodes << " " <<
            stats.time_of_last_incumbent_update << " " <<
            time_elapsed << endl;
}

