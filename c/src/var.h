#pragma once

#include <cassert>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <limits>

#include "debug_output.h"

class Var;

class VarId {
  public:
    /* This is a temporary hack to get CommutativeRExp to compile...  We
     * initialize it with max() to indicate that it's really invalid (has no
     * entry in Var).  */
    VarId() : id_(std::numeric_limits<std::uint_fast32_t>::max()) {}

    VarId(std::uint_fast32_t i) : id_(i) {}

    VarId(const VarId &v) = default;
    VarId(VarId &&v) = default;

    VarId& operator=(const VarId &v) = default;
    VarId& operator=(VarId &&v) = default;

    inline bool operator==(const VarId rhs) const { return id_ == rhs.id_; }
    inline bool operator!=(const VarId rhs) const { return id_ != rhs.id_; }
    inline bool operator<(const VarId rhs) const { return id_ < rhs.id_; }
    inline bool operator>(const VarId rhs) const { return id_ > rhs.id_; }

    VarId operator++() { ++id_; return *this; }

    // VarId operator++(int i) { return VarId{id_++}; }

    // VarId operator+(const VarId rhs) const {
    //   return VarId(id_ + rhs.id_);
    // }

    inline std::size_t Hash() const {
      std::hash<std::uint_fast32_t> h;
      return h(id_);
    }

    // friend std::ostream& operator<<(std::ostream &out, const VarId &vid) {
    //   return out << vid.id_;
    // }

    std::uint_fast32_t GetRawId() const { return id_; }

  private:
    std::uint_fast32_t id_;
};

std::ostream& operator<<(std::ostream &out, const VarId &vid);


namespace std {

template <>
struct hash<VarId> {
  inline std::size_t operator()(const VarId vid) const {
    return vid.Hash();
  }

};

}  /* namespace std */



typedef VarId VarPtr;

class Var {
  public:

    static VarId GetVarId() {
      std::stringstream ss;
      /* Prefix auto-generated variables with underscore. */
      ss << "_" << next_id_.GetRawId();
      /* GetVarId(std::string) will use next_id_ and create the new mappings. */
      return GetVarId(ss.str());
    }

    static VarId GetVarId(const std::string &name) {
      auto iter = name_to_id_.find(name);
      /* Var exists, return reference to it. */
      if (iter != name_to_id_.end()) {
        return iter->second;
      }

      /* Var doesn't exists, create a new one. */
      std::unique_ptr<Var> var{new Var(next_id_, name)};
      ++next_id_;
      auto iter_inserted = id_to_var_.emplace(var->GetId(), std::move(var));
      assert(iter_inserted.second);
      auto &inserted_var = iter_inserted.first->second;
      name_to_id_.emplace(inserted_var->GetName(), inserted_var->GetId());
      return inserted_var->GetId();
    }

    static const Var& GetVar(const VarId vid) {
      auto iter = id_to_var_.find(vid);
      assert(iter != id_to_var_.end());
      return *iter->second;
    }

    inline bool operator==(const Var &rhs) const {
      return id_ == rhs.id_;
    }

    inline bool operator<(const Var &rhs) const {
      return id_ < rhs.id_;
    }

    inline VarId GetId() const { return id_; }

    inline std::string string() const { return name_; }

    friend std::ostream& operator<<(std::ostream &out, const Var &var) {
      return out << var.string();
    }

  private:

    Var(const VarId i, const std::string &n) : id_(i), name_(n) {}
    Var(const VarId i, std::string &&n) : id_(i), name_(std::move(n)) {}

    VarId id_;
    std::string name_;

    static VarId next_id_;
    static std::unordered_map<std::string, VarId> name_to_id_;
    static std::unordered_map<VarId, std::unique_ptr<Var> > id_to_var_;

    std::string GetName() const { return name_; }

};

std::ostream& operator<<(std::ostream &out, const std::vector<VarId> vids);
