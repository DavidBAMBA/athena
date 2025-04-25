//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file harris_laminar.cpp
//! \brief Problem generator for laminar reconnection in a Harris current sheet.
//!
//! Implements a 2D resistive MHD simulation of a Harris sheet:
//!   Bx(z) = B0 * tanh(z / a)
//!   rho(z) = rho0 * sech^2(z / a) + rho_bg
//!   P = constant
//! A small perturbation in vy is added to trigger the reconnection.
//! This configuration is designed to replicate the Sweet-Parker reconnection regime
//! without turbulent forcing, as described in Loureiro et al. (2009).
//========================================================================================

#include <cmath>
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../coordinates/coordinates.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // --- Physics Parameters ---
  Real B0     = pin->GetReal("problem", "b0");
  Real a      = pin->GetReal("problem", "a");
  Real rho0   = pin->GetReal("problem", "rho0");
  Real rho_bg = pin->GetReal("problem", "rho_bg");
  Real P0     = pin->GetReal("problem", "p0");
  Real dvamp  = pin->GetReal("problem", "dvamp");
  Real sigma  = pin->GetReal("problem", "sigma");

  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        Real z = pcoord->x2v(j);
        Real x = pcoord->x1v(i);

        Real tanh_za = std::tanh(z/a);
        Real sech2_za = 1.0 / (std::cosh(z/a) * std::cosh(z/a));

        Real rho = rho0 * sech2_za + rho_bg;
        Real vx = 0.0;
        Real vy = dvamp * std::sin(2.0 * PI * x / (pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min)) *
                  std::exp(- (z*z) / (sigma*sigma));
        Real vz = 0.0;

        phydro->u(IDN, k, j, i) = rho;
        phydro->u(IM1, k, j, i) = rho * vx;
        phydro->u(IM2, k, j, i) = rho * vy;
        phydro->u(IM3, k, j, i) = rho * vz;

        if (NON_BAROTROPIC_EOS) {
          Real kinetic = 0.5 * rho * (vx*vx + vy*vy + vz*vz);
          phydro->u(IEN, k, j, i) = P0 / (pin->GetReal("hydro", "gamma") - 1.0) + kinetic;
        }

        // Bx = B0 * tanh(z/a), By = 0, Bz = 0
        pfield->b.x1f(k,j,i)   = B0 * tanh_za;
        if (i == ie) pfield->b.x1f(k,j,i+1) = B0 * tanh_za;
        pfield->b.x2f(k,j,i)   = 0.0;
        pfield->b.x3f(k,j,i)   = 0.0;

        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) += 0.5 * SQR(B0 * tanh_za);
        }
      }
    }
  }

  return;
}
