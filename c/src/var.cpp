#include "var.h"

VarId Var::next_id_ = 0;

std::unordered_map<std::string, VarId> Var::name_to_id_;
std::unordered_map<VarId, std::unique_ptr<Var> > Var::id_to_var_;

std::ostream& operator<<(std::ostream &out, const VarId &vid) {
  return out << Var::GetVar(vid);
}

std::ostream& operator<<(std::ostream &out, const std::vector<VarId> &vids) {
  out << "{";
  for (auto iter = vids.begin(); iter != vids.end(); ++iter) {
    if (iter != vids.begin()) {
      out << ",";
    }
    out << Var::GetVar(*iter);
  }
  out << "}";
  return out;
}
