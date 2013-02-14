#pragma once

#include <memory>
#include <unordered_map>

#include "var.h"
#include "hash.h"


class Node;

class Addition;
class Multiplication;
class Star;
class Element;
class Epsilon;
class Empty;

/* See [Note: Garbage] */
typedef Node* NodePtr;

class NodeVisitor;

class NodeFactory;

class Node {
  public:
    virtual ~Node() {}
    virtual void Accept(NodeVisitor &visitor) const = 0;
};

class StringPrinter;

std::ostream& operator<<(std::ostream &out, const Node &node);


/* We're providing the basic default implementations of Visit() functions, so
 * that the visitors inheriting from NodeVisitor can implement only some of
 * them, e.g., if you want to print all the variables used in a Node tree, you
 * only need to override the Visit(Element &). */
class NodeVisitor {
  public:
    virtual ~NodeVisitor() {}
    virtual void Visit(const Addition &a);
    virtual void Visit(const Multiplication &m);
    virtual void Visit(const Star &s);
    virtual void Visit(const Element &e);
    virtual void Visit(const Epsilon &e);
    virtual void Visit(const Empty &e);
};


class Addition : public Node {
  public:
    ~Addition() = default;

    void Accept(NodeVisitor &visitor) const { visitor.Visit(*this); }

    NodePtr GetLhs() const { return lhs; }
    NodePtr GetRhs() const { return rhs; }

  private:
    Addition(NodePtr l, NodePtr r) : lhs(l), rhs(r) {}
    const NodePtr lhs;
    const NodePtr rhs;
    friend class NodeFactory;
};

class Multiplication : public Node {
  public:
    ~Multiplication() = default;

    void Accept(NodeVisitor &visitor) const { visitor.Visit(*this); }

    NodePtr GetLhs() const { return lhs; }
    NodePtr GetRhs() const { return rhs; }

  private:
    Multiplication(NodePtr l, NodePtr r) : lhs(l), rhs(r) {}
    const NodePtr lhs;
    const NodePtr rhs;
    friend class NodeFactory;
};

class Star : public Node {
  public:
    ~Star() = default;

    void Accept(NodeVisitor &visitor) const { visitor.Visit(*this); }

    NodePtr GetNode() const { return node; }

  private:
    Star(NodePtr n) : node(n) {}
    const NodePtr node;
    friend class NodeFactory;
};

class Element : public Node {
  public:
    ~Element() = default;

    void Accept(NodeVisitor &visitor) const { visitor.Visit(*this); }

    VarPtr GetVar() const { return var; }

  private:
    Element(VarPtr v) : var(v) {}
    const VarPtr var;
    friend class NodeFactory;
};

class Empty : public Node {
  public:
    ~Empty() = default;
    void Accept(NodeVisitor &visitor) const { visitor.Visit(*this); }
  private:
    Empty() = default;
    friend class NodeFactory;
};

class Epsilon : public Node {
  public:
    ~Epsilon() = default;
    void Accept(NodeVisitor &visitor) const { visitor.Visit(*this); }
  private:
    Epsilon() = default;
    friend class NodeFactory;
};

class NodeFactory {
  public:
    NodeFactory() : empty_(new Empty), epsilon_(new Epsilon) {}
    virtual ~NodeFactory() {
      for (auto &pair : additions_) { delete pair.second; }
      for (auto &pair : multiplications_) { delete pair.second; }
      for (auto &pair : stars_) { delete pair.second; }
      for (auto &pair : elems_) { delete pair.second; }
      delete empty_;
      delete epsilon_;
    }

    virtual NodePtr NewAddition(NodePtr lhs, NodePtr rhs);
    virtual NodePtr NewMultiplication(NodePtr lhs, NodePtr rhs);
    virtual NodePtr NewStar(NodePtr node);
    virtual NodePtr NewElement(VarPtr var);
    virtual NodePtr GetEmpty() const { return empty_; }
    virtual NodePtr GetEpsilon() const { return epsilon_; }

    virtual void PrintDot(std::ostream &out);
    virtual void GC();

  private:
    std::unordered_map< std::pair<NodePtr, NodePtr>, NodePtr > additions_;
    std::unordered_map< std::pair<NodePtr, NodePtr>, NodePtr > multiplications_;
    std::unordered_map< NodePtr, NodePtr > stars_;
    std::unordered_map< VarPtr, NodePtr > elems_;
    NodePtr empty_;
    NodePtr epsilon_;
};

/*
 * [Note: Garbage]
 *
 * The current version of NodeFactory will never deallocate Nodes since we
 * always keep at least one reference in the std::unordered_maps.  We should
 * probably introduce external reference count (external as in not counting the
 * references kept by NodeFactory itself).
 *
 * But we have to be careful with this because we can easily end up destroying
 * and creating the same object repeatedly during the fixed-point computation.
 * One way out of this would be to perform something like a garbage collection.
 * So instead of deallocating (and removing the entry from unordered_map) as
 * soon as the external reference count is 0, do that from time to time.
 *
 */
