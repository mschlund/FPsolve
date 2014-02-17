#include <algorithm>

#include "hash.h"

#include "free-structure.h"


/*
 * NodeVisitor
 */

void NodeVisitor::Visit(const Addition &a) {
  a.GetLhs()->Accept(*this);
  a.GetRhs()->Accept(*this);
}

void NodeVisitor::Visit(const Multiplication &m) {
  m.GetLhs()->Accept(*this);
  m.GetRhs()->Accept(*this);
}

void NodeVisitor::Visit(const Star &s) {
  s.GetNode()->Accept(*this);
}

void NodeVisitor::Visit(const Element &e) {}
void NodeVisitor::Visit(const Epsilon &e) {}
void NodeVisitor::Visit(const Empty &e) {}

/*
 * StringPrinter
 */

class StringPrinter : public NodeVisitor {
  public:
    StringPrinter(std::ostream &out) : out_(out) {}
    ~StringPrinter() = default;

    void Visit(const Addition &a) {
      out_ << "(";
      a.GetLhs()->Accept(*this);
      out_ << " + ";
      a.GetRhs()->Accept(*this);
      out_ << ")";
    }

    void Visit(const Multiplication &m) {
      out_ << "(";
      m.GetLhs()->Accept(*this);
      out_ << " * ";
      m.GetRhs()->Accept(*this);
      out_ << ")";
    }

    void Visit(const Star &s) {
      out_ << "(";
      s.GetNode()->Accept(*this);
      out_ << ")*";
    }

    void Visit(const Element &e) {
      out_ << e.GetVar();
    }

    void Visit(const Epsilon &e) {
      out_ << "_";
    }

    void Visit(const Empty &e) {}

  private:
    std::ostream &out_;
};

class WordsPrinter : public NodeVisitor {
  public:
    WordsPrinter() = default;
    ~WordsPrinter() = default;

    void Visit(const Addition &a) {
      a.GetLhs()->Accept(*this);
      auto tmp = std::move(words_);
      words_ = {};
      a.GetRhs()->Accept(*this);
      std::move(tmp.begin(), tmp.end(), std::back_inserter(words_));
    }

    void Visit(const Multiplication &m) {
      m.GetLhs()->Accept(*this);
      auto tmp = std::move(words_);
      words_ = {};
      m.GetRhs()->Accept(*this);
      std::vector<std::string> result;
      for (auto &lhs_word : tmp) {
        for (auto &rhs_word : words_) {
          result.emplace_back(lhs_word + rhs_word);
        }
      }
      words_ = std::move(result);
    }

    void Visit(const Star &s) {
      s.GetNode()->Accept(*this);
      std::sort(words_.begin(), words_.end());
      std::stringstream ss;
      ss << "(";
      bool first = true;
      for (const auto &word : words_) {
        if (!first) {
          ss << "+";
        } else {
          first = false;
        }
        ss << word;
      }
      ss << ")*";
      words_ = { ss.str() };
    }

    void Visit(const Element &e) {
      words_.emplace_back(Var::GetVar(e.GetVar()).string());
    }

    void Visit(const Epsilon &e) {
      words_.emplace_back("1");
    }

    void Visit(const Empty &e) {}

    const std::vector<std::string>& GetWords() const {
      return words_;
    }

    std::string GetString() {
      std::sort(words_.begin(), words_.end());
      std::stringstream ss;
      bool first = true;
      ss << "(";
      for (const auto &word : words_) {
        if (!first) {
          ss << "+";
        } else {
          first = false;
        }
        ss << word;
      }
      ss << ")";
      return ss.str();
    }

  private:
    std::vector<std::string> words_;
};

std::string NodeToString(const Node &node) {
  WordsPrinter printer;
  node.Accept(printer);
  return printer.GetString();
}

std::string NodeToRawString(const Node &node) {
  std::stringstream ss;
  StringPrinter printer{ss};
  node.Accept(printer);
  return ss.str();
}


/*
 * NodeFactory
 *
 * The following is a bit repetitive, but it's short enough that it's probably
 * not worth the effort to abstract away the common parts...
 */

NodePtr NodeFactory::NewAddition(NodePtr lhs, NodePtr rhs) {
  assert(lhs);
  assert(rhs);

  if (lhs == empty_) {
    return rhs;
  }
  if (rhs == empty_) {
    return lhs;
  }

  /* Since + is commutative, use pointers to order the arguments.
   * This makes it possible to have that:
   *   NewAddition(a, b) == NewAddition(b, a) */
  if (lhs > rhs) {
    std::swap(lhs, rhs);
  }

  auto iter = additions_.find({lhs, rhs});
  if (iter != additions_.end()) {
    return iter->second;
  }
  NodePtr node_ptr{new Addition(lhs, rhs)};
  additions_.insert({ {lhs, rhs}, node_ptr });
  return node_ptr;
}

NodePtr NodeFactory::NewMultiplication(NodePtr lhs, NodePtr rhs) {
  assert(lhs);
  assert(rhs);

  if (lhs == epsilon_) {
    return rhs;
  }
  if (rhs == epsilon_) {
    return lhs;
  }

  if (lhs == empty_ || rhs == empty_) {
    return empty_;
  }

  auto iter = multiplications_.find({lhs, rhs});
  if (iter != multiplications_.end()) {
    return iter->second;
  }
  NodePtr node_ptr{new Multiplication(lhs, rhs)};
  multiplications_.insert({ {lhs, rhs}, node_ptr });
  return node_ptr;
}

NodePtr NodeFactory::NewStar(NodePtr node) {
  assert(node);

  if (node == empty_) {
    return epsilon_;
  }

  auto iter = stars_.find(node);
  if (iter != stars_.end()) {
    return iter->second;
  }
  NodePtr node_ptr{new Star(node)};
  stars_.insert({node, node_ptr});
  return node_ptr;
}

NodePtr NodeFactory::NewElement(VarId var) {
  auto iter = elems_.find(var);
  if (iter != elems_.end()) {
    return iter->second;
  }
  NodePtr node_ptr{new Element(var)};
  elems_.insert({var, node_ptr});
  return node_ptr;
}

void NodeFactory::GC() {
  assert(false);
}

void NodeFactory::PrintStats(std::ostream &out) {
  std::cout << "Size (free-struct): "
    << additions_.size() + multiplications_.size()+ stars_.size()
    << std::endl;
  std::cout << "Add (free-struct): " << additions_.size() << std::endl;
  std::cout << "Mult (free-struct): " << multiplications_.size() << std::endl;
  std::cout << "Stars (free-struct): " << stars_.size() << std::endl;
  std::cout << "Elems (free-struct): " << elems_.size() << std::endl;
}



void NodeFactory::PrintDot(std::ostream &out) {
  out << "digraph {" << std::endl;

  struct TypePrinter : public NodeVisitor {
    TypePrinter(std::ostream &out) : out_(out) {};
    ~TypePrinter() = default;

    void Visit(const Addition &a) { out_ << "+"; }
    void Visit(const Multiplication &m) { out_ << "*"; }
    void Visit(const Star &s) { out_ << "(-)*"; }
    void Visit(const Element &e) { out_ << e.GetVar(); }
    void Visit(const Epsilon &e) { out_ << "epsilon"; }
    void Visit(const Empty &e) { out_ << "empty"; }
    std::ostream &out_;
  };

  auto print_edge = [this,&out](NodePtr a, NodePtr b) {
    out << "\"" << a << "\""
        << " -> "
        << "\"" << b << "\""
        << std::endl;
  };

  TypePrinter type_printer{out};
  auto print_node = [&out,&type_printer](NodePtr node) {
    out << "\"" << node << "\"" << " [label=\"";
    node->Accept(type_printer);
    out << "\"]" << std::endl;
  };

  for (auto &children_parent : additions_) {
    print_node(children_parent.second);
    print_edge(children_parent.second, children_parent.first.first);
    print_edge(children_parent.second, children_parent.first.second);
  }

  for (auto &children_parent : multiplications_) {
    print_node(children_parent.second);
    print_edge(children_parent.second, children_parent.first.first);
    print_edge(children_parent.second, children_parent.first.second);
  }

  for (auto &child_parent : stars_) {
    print_node(child_parent.second);
    print_edge(child_parent.second, child_parent.first);
  }

  for (auto &child_parent : elems_) {
    out << "\"" << child_parent.second << "\""
        << " [label=\"" << child_parent.first << "\"]"
        << std::endl;
  }

  print_node(empty_);
  print_node(epsilon_);

  out << "}" << std::endl;
}

