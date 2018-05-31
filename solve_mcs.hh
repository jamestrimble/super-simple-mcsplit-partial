#include <atomic>
#include <chrono>

extern std::atomic<bool> abort_due_to_timeout;

enum class Heuristic
{
    min_max,
    min_product
};

struct Params
{
    bool quiet;
    bool verbose;
    bool mcsplit_down;
    Heuristic heuristic;
    std::chrono::time_point<std::chrono::steady_clock> start_time;
};

struct Assignment
{
    int v;
    int w;

    auto operator==(const Assignment & other) const -> bool
    {
        return v==other.v && w==other.w;
    }
};

struct Edge
{
    int v;
    int w;
};

struct McsStats
{
    unsigned long long nodes = 0;
    long long time_of_last_incumbent_update = 0;
};

auto solve_mcs(Graph & linegraph0, Graph & linegraph1,
        const std::vector<Edge> & edges0,
        const std::vector<Edge> & edges1,
        const Graph & originalgraph0,
        const Graph & originalgraph1,
        Params params)
		-> std::pair<std::vector<Assignment>, McsStats>;
