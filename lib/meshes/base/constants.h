/**
 * Header file containing physical constants
 *
 * @author: Hans R. Hammer
 */

#pragma once

#include <cmath>
#include <limits>

namespace murmeln {
namespace constants {
/// infinity
const double INF = std::numeric_limits<double>::infinity();

/// maximum finite positive double
constexpr double Max_Dbl = std::numeric_limits<double>::max();

constexpr double pi = 3.14159265358979323846;

/// not a number
//    const double NAN = std::numeric_limits<double>::nan();
/// Pi, circle constant [-]
constexpr long double PI =
    3.1415926535897932384626433832795028841971693993751058;
/// Speed of light c [m/s]
//    const long double SPEED_LIGHT = 299792458.0;
/// Speed of light c [cm/s]
const double SPEED_LIGHT = 2.99792458e10;

/// Natural charge Q [C]
const long double CHARGE = 1.6021766208e-19;

/// Gas Constant R [J K^-1 mol^-1]
//    const long double GAS_CONSTANT = 8.3144598;
/// Gas Constant R [erg eV^-1 mol^-1]
const long double GAS_CONSTANT = 8.3144598 * 11604.505 * 1.0e7;

/// Avogadro constant [mol^-1]
const long double AVOGADRO = 6.022140857e23;

/// Planck constant [J s]
//    const long double PLANK = 6.626070040e-34;
/// Planck constant [eVs]
//    const long double PLANK = 4.135667662e-15;
/// Planck constant [erg s]
const long double PLANK = 6.626070040e-34 * 1.0e7;

/// Boltzmann constant [J / K]
//    const long double BOLTZMANN = 1.38064852e-23;
/// Boltzmann constant [eV/K]
//    const long double BOLTZMANN = 8.6173303e-5;
const long double BOLTZMANN = GAS_CONSTANT / AVOGADRO;

/// Stefan-Boltzmann constant sigma [J s^-1 m^-2 K^-4]
//    const long double STEFAN_BOLTZMANN = 5.670367e-8;
const long double STEFAN_BOLTZMANN =
    2.0 * std::pow(PI, 5) * std::pow(BOLTZMANN, 4) /
    (15.0 * std::pow(PLANK, 3) * std::pow(SPEED_LIGHT, 2));

/// Radiation constant a [J/(m^3K^4)]
//    const long double RADIATION = 7.5657e-16;
const long double RADIATION = 4.0 * STEFAN_BOLTZMANN / SPEED_LIGHT;

/**
 * @brief Planck spectrum
 *
 * The Planck spectrum dependent on the temperature and frequency
 *
 * @param _temperature Temperature [eV]
 * @param _frequency Frequency [s^-1]
 * @return
 */
inline long double planckFct(const long double &_temperaure,
                             const long double &_frequency) {
  return 2.0 * PLANK / std::pow(SPEED_LIGHT, 2) * std::pow(_frequency, 3) *
         1.0 /
         (std::exp(-PLANK * _frequency / (BOLTZMANN * _temperaure)) - 1.0);
}

} // namespace constants
} // namespace murmeln
