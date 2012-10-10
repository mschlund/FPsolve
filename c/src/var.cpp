#include <assert.h>
#include "var.h"

Var::Var()
{
	this->id = max_id;
	std::stringstream ss;
	ss << "_"; // prefix for auto-generated variables
	ss << max_id;
	this->name = ss.str();
	Var::max_id++;
}

Var::Var(std::string name)
{
	assert(name[0] != '_'); // name must not begin with underscore
	this->id = max_id++;
	this->name = name;
}

std::string Var::getName()
{
	return this->name;
}

// generates a new unnamed var
VarPtr Var::getVar()
{
	VarPtr var(new Var());
	Var::vars.insert(Var::vars.begin(), std::pair<std::string,VarPtr>(var->getName(), var));
	return var;
}

// returns a reference to 
VarPtr Var::getVar(std::string name)
{
	auto v_it = Var::vars.find(name);
	if(v_it != Var::vars.end()) // var exists, return reference to it
		return v_it->second;
	else // create a new one
	{
		VarPtr var(new Var(name));
		Var::vars.insert(Var::vars.begin(), std::pair<std::string,VarPtr>(name, var));
		return var;
	}
}

VarPtr Var::getVar(VarPtr var)
{
	return var;
}

bool operator==(const VarPtr& l, const VarPtr& r)
{
	return l->id == r->id;
}

bool operator!=(const VarPtr& l, const VarPtr& r)
{
	return l->id != r->id;
}

bool operator<(const VarPtr& l, const VarPtr& r)
{
	return l->id < r->id;
}

std::string Var::string() const
{
	std::stringstream ss;
	//ss << "Var(" << this->name << ")";
	ss << this->name;
	return ss.str();
}

int Var::max_id = 0;
std::map<std::string, VarPtr> Var::vars;

std::ostream& operator<<(std::ostream& os, const VarPtr var)
{
	return os << var->string();
}

std::ostream& operator<<(std::ostream& os, const std::multiset<VarPtr> vars)
{
	std::stringstream ss;
	ss << "{";
	for(std::multiset<VarPtr>::const_iterator var = vars.begin(); var != vars.end(); ++var)
	{
		if(var != vars.begin() )
		{
			ss << ",";
		}
		ss << (*var)->string();
	}
	ss << "}";
	return os << ss.str();
}

std::ostream& operator<<(std::ostream& os, const std::vector<VarPtr>& vars)
{
	std::stringstream ss;
	ss << "{";
	for(std::vector<VarPtr>::const_iterator var = vars.begin(); var != vars.end(); ++var)
	{
		if(var != vars.begin())
		{
			ss << ",";
		}
		ss << (*var)->string();
	}
	ss << "}";
	return os << ss.str();
}
