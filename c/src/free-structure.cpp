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

std::ostream& operator<<(std::ostream &out, const Node &node) {
  StringPrinter printer{out};
  node.Accept(printer);
  return out;
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
