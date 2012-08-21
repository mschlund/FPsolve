#ifndef VAR_H
#define VAR_H

#include <string>
#include <map>
#include <sstream>
#include <set>
#include <vector>

class Var
{
private:
	int id;
	std::string name;
	static int max_id;
	static std::map<std::string, int>* vars; // name → id
public:
	Var();
	// if there is already a variable with this name, return a reference to this variable
	Var(std::string name);
	bool operator==(const Var& var) const;
	bool operator!=(const Var& var) const;
	bool operator<(const Var& var) const;
	std::string string() const;
};

std::ostream& operator<<(std::ostream& os, const Var var);
std::ostream& operator<<(std::ostream& os, const std::multiset<Var> vars);
std::ostream& operator<<(std::ostream& os, const std::vector<Var> vars);

#endif
