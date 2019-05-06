// constants.hh
// T. M. Kelley
// Jan 11, 2011
// Header for constants
// (c) Copyright 2011 LANSLLC all rights reserved.

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>

namespace nut {

double const pi = 3.141592653589793;

double const huge = 1e300;

double const tiny = 1e-100;

double const c = 2.99792458e10;  // cm/sec

double const k_B = 8.61734e-11;  // MeV/K

double const pmg = 1.67262e-24;  // proton mass in grams

std::string const time_unit = "second";       // 1e-8 sec
std::string const space_unit = "centimeter";  // 1e8 Angstrom

}  // namespace nut

#endif

// version
// $Id$

// End of file
