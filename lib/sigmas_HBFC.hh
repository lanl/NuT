// sigmas_HBFC.hh
// T. M. Kelley
// Jan 12, 2011
// (c) Copyright 2011 LANSLLC, all rights reserved

/*!\file Neutrino process cross sections as used in Herant, Benz, Fryer, and
 * Colgate, Ap. J. 435:339-361 (1994). */

namespace nut {
namespace sigmas {
/** cross sections in cm^2, energies in MeV */

/** neutrino--baryon events */

// nucleon absorption/emission
template <typename fp_t>
fp_t
nu_N_abs(fp_t const e)
{
  return 9e-44 * e * e;
}

template <typename fp_t>
fp_t
nu_N_elastic(fp_t const e)
{
  return 1.7e-44 * e * e;
}

template <typename fp_t>
fp_t
nu_N_Z_elastic(fp_t const e, fp_t const N, fp_t const Z)
{
  return 1.7e-44 * 1. / 6. * (N - 0.08 * Z) * (N - 0.08 * Z) * e * e;
}

template <typename fp_t>
fp_t
nu_A_elastic(fp_t const e, fp_t const A)
{
  return 1.7e-44 * 1. / 6. * A * A * e * e;
}

// nuetrino-lepton scattering
// nu_e--e^-  and nu_e_bar--e^+
template <typename fp_t>
fp_t
nu_e_e_minus(fp_t const e_nu, fp_t const e_e)
{
  return 9.2e-45 * e_nu * e_e;
}

template <typename fp_t>
fp_t
nu_e_bar_e_plus(fp_t const e_nu, fp_t const e_e)
{
  return 9.2e-45 * e_nu * e_e;
}

// nu_bar_e--e^- & nu_e--e^+ scattering
template <typename fp_t>
fp_t
nu_e_e_plus(fp_t const e_nu, fp_t const e_e)
{
  return 3.9e-45 * e_nu * e_e;
}

template <typename fp_t>
fp_t
nu_e_bar_e_minus(fp_t const e_nu, fp_t const e_e)
{
  return 3.9e-45 * e_nu * e_e;
}

// nu_x--e^- & nu_bar_x--e^+ scattering
template <typename fp_t>
fp_t
nu_x_e_minus(fp_t const e_nu, fp_t const e_e)
{
  return 1.8e-45 * e_nu * e_e;
}

template <typename fp_t>
fp_t
nu_x_bar_e_plus(fp_t const e_nu, fp_t const e_e)
{
  return 1.8e-45 * e_nu * e_e;
}

// nu_bar_x--e^- & nu_x--e^+ scattering
template <typename fp_t>
fp_t
nu_x_e_plus(fp_t const e_nu, fp_t const e_e)
{
  return 1.3e-45 * e_nu * e_e;
}

template <typename fp_t>
fp_t
nu_x_bar_e_minus(fp_t const e_nu, fp_t const e_e)
{
  return 1.3e-45 * e_nu * e_e;
}

// nu--nu_bar annhilation
template <typename fp_t>
fp_t
nu_e_nu_bar_e(fp_t const e_nu, fp_t const e_nu_bar)
{
  return 6.6e-45 * e_nu * e_nu_bar;
}

template <typename fp_t>
fp_t
nu_x_nu_bar_x(fp_t const e_nu, fp_t const e_nu_bar)
{
  return 1.4e-45 * e_nu * e_nu_bar;
}

/** Total cross sections */
template <typename fp_t>
fp_t
nu_N_total(fp_t const e)
{
  return nu_N_abs(e) + nu_N_elastic(e);
}

template <typename fp_t>
fp_t
nu_nuclei_total(fp_t const e, fp_t const A)
{
  return nu_A_elastic(e, A);
}

// electron (anti)neutrino--electron/positron scattering
template <typename fp_t>
fp_t
nu_e_e_total(fp_t const e_nu, fp_t const e_e)
{
  return nu_e_e_minus(e_nu, e_e) + nu_e_e_plus(e_nu, e_e);
}

// other flavor (anti)neutrino--electron/positron scattering
template <typename fp_t>
fp_t
nu_x_e_total(fp_t const e_nu, fp_t const e_e)
{
  return nu_x_e_minus(e_nu, e_e) + nu_x_e_plus(e_nu, e_e);
}

// initial implementation for total collision cross section will
// only use neutrino--nucleon scatterings (evading question of how
// to sample antineutrino/electron/positron energies).
template <typename fp_t>
fp_t
collide(fp_t const e_nu)
{
  return nu_N_total(e_nu);
}

}  // namespace sigmas

}  // namespace nut

// version
// $Id$

// End of file
