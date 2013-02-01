#include <cassert>

#include "free-structure.h"


/*
 * NodeVisitor
 */

void NodeVisitor::Visit(Addition &a) {
  a.GetLhs()->Accept(*this);
  a.GetRhs()->Accept(*this);
}

void NodeVisitor::Visit(Multiplication &m) {
  m.GetLhs()->Accept(*this);
  m.GetRhs()->Accept(*this);
}

void NodeVisitor::Visit(Star &s) {
  s.GetNode()->Accept(*this);
}

void NodeVisitor::Visit(Element &e) {}
void NodeVisitor::Visit(Epsilon &e) {}
void NodeVisitor::Visit(Empty &e) {}

/*
 * NodeFactory
 *
 * The following is a bit repetitive, but it's short enough that it's probably
 * not worth the effort to abstract away the common parts...
 */

NodePtr NodeFactory::NewAddition(NodePtr lhs, NodePtr rhs) {
  if (lhs == empty_) {
    return rhs;
  }
  if (rhs == empty_) {
    return lhs;
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
  if (lhs == epsilon_) {
    return rhs;
  }
  if (rhs == epsilon_) {
    return lhs;
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
  if (node == empty_ || node == epsilon_) {
    return node;
  }

  auto iter = stars_.find(node);
  if (iter != stars_.end()) {
    return iter->second;
  }
  NodePtr node_ptr{new Star(node)};
  stars_.insert({node, node_ptr});
  return node_ptr;
}

NodePtr NodeFactory::NewElement(VarPtr var) {
  auto iter = elems_.find(var);
  if (iter != elems_.end()) {
    return iter->second;
  }
  NodePtr node_ptr{new Element(var)};
  elems_.insert({var, node_ptr});
  return node_ptr;
}
