// T. M. Kelley (c) 2011 LANS LLC

#include "copyright.hh"
#include "gtest/gtest.h"
#include "test_aux.hh"
#include "types.hh"

TEST(species_name,nu_e)
{
  bool passed(true);

  std::cout << nut::copyright() << std::endl;

  std::string n = nut::species_name(nut::nu_e);

  passed = n == "nu_e";

  EXPECT_TRUE(passed);
  return;
}

// End of file
