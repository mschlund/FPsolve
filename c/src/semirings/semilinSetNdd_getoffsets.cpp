#include <iostream>
#include <cassert>
#include <algorithm>
#include <fstream>
#include "semilinSetNdd_getoffsets.h"

// http://www.technical-recipes.com/2011/a-recursive-algorithm-to-find-all-paths-between-two-given-nodes/

std::vector<std::vector<Match>> appendToAll(std::vector<std::vector<Match>> offsets, Match current) {
  std::vector<std::vector<Match>> result;
  for(auto offset : offsets) {
    auto tmp = offset;
    tmp.push_back(current);
    result.push_back(tmp);
  }
  if(offsets.size() == 0)
    result = {{current}};
  return result;
};

std::vector<std::vector<int>> appendElementToAll(std::vector<std::vector<int>> offsets, int current) {
  std::vector<std::vector<int>> result;
  for(auto offset : offsets) {
    auto tmp = offset;
    tmp.push_back(current);
    result.push_back(tmp);
  }
  if(offsets.size() == 0)
    result = {{current}};
  return result;
};

std::vector<int> addTo(std::vector<int> first, std::vector<int> second) {
  auto result = first;
  assert(first.size() == second.size());
  for(int i = 0; i < first.size(); ++i) {
    result.at(i) += second.at(i);
  }
  return result;
}

std::set<std::vector<int>> DepthFirst(Graph* graph, std::vector<int>& visited, std::vector<bool>& visited_bool, int end, int k) {
  int back = visited.back();
  auto adjNodes = graph->nodes[back].neighbours;

  std::vector<std::vector<Match>> offsets;
  // examine adjacent nodes
  for(auto node : adjNodes) {
    if(visited_bool.at(node.first))
      continue; // we already visited the node.
    if(node.first == end) {
      visited.push_back(node.first);
      visited_bool.at(node.first) = true;
      int hops = (int) visited.size();
      std::vector<std::vector<Match>> tmp; // build up the path of matches in tmp
      // std::cout << visited[0];
      for(int i=1; i<hops; i++) {
        // std::cout << "→" << visited[i];
        // all the possibles matches for this link
        auto matches = graph->nodes[visited[i-1]].neighbours[visited[i]];

        // tmp holds the path until this node
        auto path_before = tmp;
        std::vector<std::vector<Match>> path_after;
        // loop over all possible matches
        for(auto match : matches) {
          // ... and append them to the path until this node
          auto tmp2 = appendToAll(path_before, match);
          // then collect this new path
          path_after.insert(path_after.begin(), tmp2.begin(), tmp2.end());
        }
        // save the newly calculated path until the next node
        tmp = path_after;
      }
      // std::cout << ":\t";
      // tmp now holds one complete path to an end node, collect this path in offsets
      offsets.insert(offsets.end(), tmp.begin(), tmp.end());

      visited_bool.at(visited.back()) = false;
      visited.pop_back(); // delete last element → backtracking
      break;
    }

  }

  std::set<std::vector<int>> result;
  int i = 0;
  for(auto offset : offsets) {
    std::vector<int> sum(k);
    int pow = 0;
    for(auto match : offset) {
      int i = 0;
      for(auto m : match) {
        sum.at(i++) += m*(1<<pow); // + m*2^pow
      }
      pow++;
    }
    for(auto s : sum)
    result.insert(sum);
  }

  // in breadt-first, recursion needs to come after visiting adjacent nodes
  for(auto node : adjNodes) {
    if(visited_bool.at(node.first) || node.first == end)
      continue; // we already visited the node.
    visited.push_back(node.first);
    visited_bool.at(node.first) = true;
    auto ret = DepthFirst(graph, visited, visited_bool, end, k);
    result.insert(ret.begin(), ret.end());

    visited_bool.at(visited.back()) = false;
    visited.pop_back(); // delete last element → backtracking
  }

  return result;
}

Parsed parse_automaton(int k, std::string filename) {
  std::vector<int> finals;
  std::ifstream file(filename);
  int init;
  int number;
  int node_count;
  Graph graph;
  while(file >> number)
  {
    finals.push_back(number);
    // std::cout << number << ",";
    if(file.peek() == '\n')
      break;
  }
  // std::cout << std::endl;
  file >> init;
  // std::cout << init << std::endl;
  int first, second;
  while(file >> first >> second) {
    // std::cout << first << "→" << second << ":";

    char tmp;
    std::vector<int> input; // accumulate input before using it
    while(file >> tmp) {
      switch(tmp) {
        case '0': // std::cout << 0;
                  input.push_back(0);
                  break;
        case '1': // std::cout << 1;
                  input.push_back(1);
                  break;
        case 'X': // std::cout << 'X';
                  input.push_back(2); // 2 is 0 or 1
                  break;
        default:
                  assert(false); // this should not happen
      }

      if(file.peek() == '\n')
        break;
    }

    // process input. we have k variables and input.size inputs
    // => n = input.size()/k possible matches
    int n = input.size()/k;
    Matches matches(n);
    for(int i = 0; i < input.size(); ++i) {
      matches.at(i%n).push_back(input.at(i));
    }

    // postprocess parsed matches for Xs
    std::vector<std::vector<int>> accu;
    std::vector<std::vector<int>> matches2;
    for(int i = 0; i < matches.size(); ++i) {
      assert(matches.at(i).size() == k);
      accu = {};
      for(int j = 0; j < matches.at(i).size(); ++j) {
        if(matches.at(i).at(j) == 2) {
          auto tmp0 = appendElementToAll(accu, 0);
          auto tmp1 = appendElementToAll(accu, 1);
          accu = tmp0;
          accu.insert(accu.end(), tmp1.begin(), tmp1.end());
        } else {
          accu = appendElementToAll(accu, matches.at(i).at(j));
        }
      }
      matches2.insert(matches2.end(), accu.begin(), accu.end());
    }
    graph.nodes[first].neighbours[second] = matches2;
    // std::cout << std::endl;
    node_count = std::max(std::max(first, second), node_count);
  }
  return {init, finals, graph, node_count};
}
