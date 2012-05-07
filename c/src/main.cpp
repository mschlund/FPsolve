#include <iostream>
#include "float-semiring.h"

int main(int argc, char* argv[])
{
	std::cout << "hello world" << std::endl;
	FloatSemiring foo(5);
	FloatSemiring bar(6);
	FloatSemiring baz = foo + bar;
	std::cout << foo.getString() << " + " << bar.getString() << " = " << baz.getString() << std::endl;
	return 0;
}
