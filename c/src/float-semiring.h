#include <string>

class FloatSemiring
{
private:
	float val;
public:
	FloatSemiring(const float val);
	~FloatSemiring();
	FloatSemiring operator + (const FloatSemiring elem);
	FloatSemiring operator * (const FloatSemiring elem);
	std::string getString();
};
