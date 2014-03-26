#include <vector>
#include <map>
#include <set>

typedef std::vector<int> Match;
typedef std::vector<Match> Matches;

struct Node {
  std::map<int, Matches> neighbours;
};
struct Graph {
  std::map<int, Node> nodes;
};

struct Parsed {
  int init;
  std::vector<int> finals;
  Graph graph;
};

Parsed parse_automaton(int k, std::string filename);
std::set<std::vector<int>> DepthFirst(Graph* graph, std::vector<int>& visited, int end, int k);
