// Local Lax-Friedrichs Riemann solver for relativistic hydrodynamics

// Primary header
#include "../../hydro_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../hydro.hpp"                       // Hydro
#include "../../../eos/eos.hpp"                     // HydroEqnOfState
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../../mesh.hpp"                     // MeshBlock
#include "../../../../coordinates/coordinates.hpp"  // Coordinates

// Riemann solver
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   bb: 3D array of normal magnetic fields (not used)
//   prim_l, prim_r: left and right primitive states
// Outputs:
//   flux: fluxes across interface
// Notes:
//   prim_l, prim_r overwritten
//   implements LLF algorithm similar to that of fluxcalc() in step_ch.c in Harm
//   references Mignone & Bodo 2005, MNRAS 364 126 (MB)
void HydroIntegrator::RiemannSolver(const int k, const int j, const int il,
    const int iu, const int ivx, const AthenaArray<Real> &bb, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &flux)
{
  // Transform primitives to locally flat coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_hydro->pmy_block->pcoord->PrimToLocal1(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
      case IVY:
        pmy_hydro->pmy_block->pcoord->PrimToLocal2(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
      case IVZ:
        pmy_hydro->pmy_block->pcoord->PrimToLocal3(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
    }

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_hydro->peos->GetGamma();

  // Go through each interface
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract left primitives
    const Real &rho_l = prim_l(IDN,i);
    const Real &pgas_l = prim_l(IEN,i);
    Real u_l[4];
    if (GENERAL_RELATIVITY)
    {
      u_l[1] = prim_l(ivx,i);
      u_l[2] = prim_l(ivy,i);
      u_l[3] = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 + SQR(u_l[1]) + SQR(u_l[2]) + SQR(u_l[3]));
    }
    else  // SR
    {
      const Real &vx_l = prim_l(ivx,i);
      const Real &vy_l = prim_l(ivy,i);
      const Real &vz_l = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 / (1.0 - SQR(vx_l) - SQR(vy_l) - SQR(vz_l)));
      u_l[1] = u_l[0] * vx_l;
      u_l[2] = u_l[0] * vy_l;
      u_l[3] = u_l[0] * vz_l;
    }

    // Extract right primitives
    const Real &rho_r = prim_r(IDN,i);
    const Real &pgas_r = prim_r(IEN,i);
    Real u_r[4];
    if (GENERAL_RELATIVITY)
    {
      u_r[1] = prim_r(ivx,i);
      u_r[2] = prim_r(ivy,i);
      u_r[3] = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 + SQR(u_r[1]) + SQR(u_r[2]) + SQR(u_r[3]));
    }
    else  // special relativity
    {
      const Real &vx_r = prim_r(ivx,i);
      const Real &vy_r = prim_r(ivy,i);
      const Real &vz_r = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 / (1.0 - SQR(vx_r) - SQR(vy_r) - SQR(vz_r)));
      u_r[1] = u_r[0] * vx_r;
      u_r[2] = u_r[0] * vy_r;
      u_r[3] = u_r[0] * vz_r;
    }

    // Calculate wavespeeds in left state (MB 23)
    Real lambda_p_l, lambda_m_l;
    Real wgas_l = rho_l + gamma_adi/(gamma_adi-1.0) * pgas_l;
    pmy_hydro->peos->SoundSpeedsSR(wgas_l, pgas_l, u_l[1]/u_l[0], SQR(u_l[0]),
        &lambda_p_l, &lambda_m_l);

    // Calculate wavespeeds in right state (MB 23)
    Real lambda_p_r, lambda_m_r;
    Real wgas_r = rho_r + gamma_adi/(gamma_adi-1.0) * pgas_r;
    pmy_hydro->peos->SoundSpeedsSR(wgas_r, pgas_r, u_r[1]/u_r[0], SQR(u_r[0]),
        &lambda_p_r, &lambda_m_r);

    // Calculate extremal wavespeed
    Real lambda_l = std::min(lambda_m_l, lambda_m_r);
    Real lambda_r = std::max(lambda_p_l, lambda_p_r);
    Real lambda = std::max(lambda_r, -lambda_l);

    // Calculate conserved quantities in L region (MB 3)
    Real cons_l[NWAVE];
    cons_l[IDN] = rho_l * u_l[0];
    cons_l[IEN] = wgas_l * u_l[0] * u_l[0] - pgas_l;
    cons_l[ivx] = wgas_l * u_l[1] * u_l[0];
    cons_l[ivy] = wgas_l * u_l[2] * u_l[0];
    cons_l[ivz] = wgas_l * u_l[3] * u_l[0];

    // Calculate fluxes in L region (MB 2,3)
    Real flux_l[NWAVE];
    flux_l[IDN] = rho_l * u_l[1];
    flux_l[IEN] = wgas_l * u_l[0] * u_l[1];
    flux_l[ivx] = wgas_l * u_l[1] * u_l[1] + pgas_l;
    flux_l[ivy] = wgas_l * u_l[2] * u_l[1];
    flux_l[ivz] = wgas_l * u_l[3] * u_l[1];

    // Calculate conserved quantities in R region (MB 3)
    Real cons_r[NWAVE];
    cons_r[IDN] = rho_r * u_r[0];
    cons_r[IEN] = wgas_r * u_r[0] * u_r[0] - pgas_r;
    cons_r[ivx] = wgas_r * u_r[1] * u_r[0];
    cons_r[ivy] = wgas_r * u_r[2] * u_r[0];
    cons_r[ivz] = wgas_r * u_r[3] * u_r[0];

    // Calculate fluxes in R region (MB 2,3)
    Real flux_r[NWAVE];
    flux_r[IDN] = rho_r * u_r[1];
    flux_r[IEN] = wgas_r * u_r[0] * u_r[1];
    flux_r[ivx] = wgas_r * u_r[1] * u_r[1] + pgas_r;
    flux_r[ivy] = wgas_r * u_r[2] * u_r[1];
    flux_r[ivz] = wgas_r * u_r[3] * u_r[1];

    // Set fluxes
    for (int n = 0; n < NWAVE; ++n)
      flux(n,i) = 0.5 * (flux_l[n] + flux_r[n] - lambda * (cons_r[n] - cons_l[n]));

    // Set conserved quantities in GR
    if (GENERAL_RELATIVITY)
      for (int n = 0; n < NWAVE; ++n)
        cons_(n,i) = 0.5 * (cons_r[n] + cons_l[n] + (flux_l[n] - flux_r[n]) / lambda);
  }

  // Transform fluxes to global coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_hydro->pmy_block->pcoord->FluxToGlobal1(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
      case IVY:
        pmy_hydro->pmy_block->pcoord->FluxToGlobal2(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
      case IVZ:
        pmy_hydro->pmy_block->pcoord->FluxToGlobal3(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
    }
  return;
}