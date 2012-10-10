#ifndef VAR_H
#define VAR_H

#include <string>
#include <map>
#include <sstream>
#include <set>
#include <vector>
#include <memory>

class Var;
typedef std::shared_ptr<Var> VarPtr;

class Var
{
private:
	int id;
	std::string name;
	static int max_id;
	static std::map<std::string, VarPtr> vars; // name → Var*
	Var();
	Var(std::string name);
	std::string getName();
public:
	static VarPtr getVar();
	static VarPtr getVar(std::string name);
	static VarPtr getVar(VarPtr var);
	friend bool operator==(const VarPtr& l, const VarPtr& r);
	friend bool operator!=(const VarPtr& l, const VarPtr& r);
	friend bool operator<(const VarPtr& l, const VarPtr& r);
	std::string string() const;
};

std::ostream& operator<<(std::ostream& os, const VarPtr var);
std::ostream& operator<<(std::ostream& os, const std::multiset<VarPtr> vars);
std::ostream& operator<<(std::ostream& os, const std::vector<VarPtr> vars);

#endif
