#pragma once
#ifdef OPCOUNT
#undef OPCOUNT
#include <iostream>
#define OPADD std::cout << "Addition on line: " << __LINE__ << " in file " << __FILE__ << std::endl
#define OPMULT std::cout << "Multiplication on line: " << __LINE__ << " in file " << __FILE__ << std::endl
#define OPSTAR std::cout << "Star on line: " << __LINE__ << " in file " << __FILE__ << std::endl
//#define OPADD __count_add++
//#define OPMULT __count_mult++
//#define SHOWOPS std::cout << std::endl << "Statistics:" << std::endl << "\tAdditions:\t" << __count_add << std::endl << "\tMultiplications:\t" << __count_mult << std::endl
#else
#define OPADD
#define OPMULT
#define OPSTAR
//#define SHOWOPS
#endif // OPCOUNT
