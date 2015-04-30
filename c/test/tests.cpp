#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#ifdef USE_GENEPI
#include "test-semilinSetNdd.h"
#endif

int main(int argc, char* argv[])
{
#ifdef USE_GENEPI
  SemilinSetNdd::genepi_init();
#endif

  // Get the top level suite from the registry
  CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();

  // Adds the test to the list of test to run
  CppUnit::TextUi::TestRunner runner;
  runner.addTest( suite );

  // Change the default outputter to a compiler error format outputter
  runner.setOutputter( new CppUnit::CompilerOutputter( &runner.result(), std::cerr ) );

  // Run the tests.
  bool wasSucessful = runner.run();

  //SemilinSetNdd::genepi_dealloc();
  // Return error code 1 if the one of test failed.
  return wasSucessful ? 0 : 1;
}

