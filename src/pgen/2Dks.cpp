//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file rt.c
//  \brief Problem generator for RT instabilty.
//
// Note the gravitational acceleration is hardwired to be 0.1. Density difference is
// hardwired to be 2.0 in 2D, and is set by the input parameter <problem>/rhoh in 3D
// (default value is 3.0). This reproduces 2D results of Liska & Wendroff, 3D results of
// Dimonte et al.
//
// FOR 2D HYDRO:
// Problem domain should be -1/6 < x < 1/6; -0.5 < y < 0.5 with gamma=1.4 to match Liska
// & Wendroff. Interface is at y=0; perturbation added to Vy. Gravity acts in y-dirn.
// Special reflecting boundary conditions added in x2 to improve hydrostatic eqm
// (prevents launching of weak waves) Atwood number A=(d2-d1)/(d2+d1)=1/3. Options:
//     iprob = 1  -- Perturb V2 using single mode
//     iprob != 1 -- Perturb V2 using multiple mode
//
// FOR 3D:
// Problem domain should be -.05 < x < .05; -.05 < y < .05, -.1 < z < .1, gamma=5/3 to
// match Dimonte et al.  Interface is at z=0; perturbation added to Vz. Gravity acts in
// z-dirn. Special reflecting boundary conditions added in x3.  A=1/2.  Options:
//     iprob = 1 -- Perturb V3 using single mode
//     iprob = 2 -- Perturb V3 using multiple mode
//     iprob = 3 -- B rotated by "angle" at interface, multimode perturbation
//
// REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
//========================================================================================

#include <cmath>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);

// made global to share with BC functions
static Real grav_acc;

int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive==true)
    EnrollUserRefinementCondition(RefinementCondition);
  if (mesh_size.nx3 == 1) {  // 2D problem
    // Enroll special BCs
    EnrollUserBoundaryFunction(INNER_X2, ProjectPressureInnerX2);
    EnrollUserBoundaryFunction(OUTER_X2, ProjectPressureOuterX2);
  } else { // 3D problem
    // Enroll special BCs
    EnrollUserBoundaryFunction(INNER_X3, ProjectPressureInnerX3);
    EnrollUserBoundaryFunction(OUTER_X3, ProjectPressureOuterX3);
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the Rayleigh-Taylor instability test
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  int64_t iseed = -1;
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;
  grav_acc = phydro->psrc->GetG2();

  Real kx = 2.0*(PI)/(pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real ky = 2.0*(PI)/(pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  Real kz = 2.0*(PI)/(pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min);

  // Read perturbation amplitude, problem switch, density ratio
  Real amp = pin->GetReal("problem","amp");
  Real sigma = pin->GetReal("problem","sigma");
  Real delta = pin->GetReal("problem","delta");
  int nx = pin->GetReal("problem","nx");
  int ny = pin->GetReal("problem","ny");

  Real drat = pin->GetOrAddReal("problem","drat",3.0);
  
  Real d0 = 1.0;			//Density in the cold magnetized layer
  Real b0 = sqrt(sigma*d0);		//b0 = B0/sqrt(4*PI)
  Real p0 = 0.5*b0*b0;			//d0 and p0 chosen that sound speed = 1
  Real h0 = drat*d0 + gamma*p0/gm1;	//Enthalpy density in the hot layer
  Real v1y;


// 2D PROBLEM ---------------------------------------------------------------

  if (block_size.nx3 == 1) {
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {

        //Velocity perturbation (* ZERO NET MOMENTUM *) -----------------------//
        v1y = amp*0.5*sin(kx*nx*pcoord->x1v(i))*(1+cos(ky*ny*pcoord->x2v(j)));
        //---------------------------------------------------------------------//

        phydro->u(IDN,k,j,i) = drat*d0;
        phydro->u(IEN,k,j,i) = drat*d0+p0/gm1;
        phydro->u(IM1,k,j,i) = 0.0;
        // Velocity perturbation in x and y direction
        phydro->u(IM2,k,j,i) = h0*v1y;
        //------------------------------------------//
        phydro->u(IM3,k,j,i) = 0.0;
        if (MAGNETIC_FIELDS_ENABLED) {
          pfield->b.x1f(k,j,i) = 0.0;
          pfield->b.x2f(k,j,i) = 0.0;
          pfield->b.x3f(k,j,i) = 0.0;
        }
        // Lower B-field Zone ----------------------//
        if (pcoord->x2v(j) <= -0.5*delta) {
          phydro->u(IDN,k,j,i) = d0;
          phydro->u(IEN,k,j,i) = d0+p0*(1+v1y*v1y);
          phydro->u(IM2,k,j,i) = (d0+2*p0)*v1y;
          if (MAGNETIC_FIELDS_ENABLED) pfield->b.x3f(k,j,i) = b0;      
        }
        // Upper B-field Zone ----------------------//
        if (pcoord->x2v(j) >= 0.5*delta) {
          phydro->u(IDN,k,j,i) = d0;
          phydro->u(IEN,k,j,i) = d0+p0*(1+v1y*v1y);
          phydro->u(IM2,k,j,i) = (d0+2*p0)*v1y;
          if (MAGNETIC_FIELDS_ENABLED) pfield->b.x3f(k,j,i) = -b0; 
        }
      }
    }}

  return;
}


//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureInnerX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      if (n==(IVY)) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(IVY,k,js-j,i) = -prim(IVY,k,js+j-1,i);  // reflect 2-velocity
        }
      } else if (n==(IPR)) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(IPR,k,js-j,i) = prim(IPR,k,js+j-1,i)
             - prim(IDN,k,js+j-1,i)*grav_acc*(2*j-1)*pco->dx2f(j);
        }
      } else {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,js-j,i) = prim(n,k,js+j-1,i);
        }
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(js-j),i) =  b.x1f(k,(js+j-1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(js-j),i) = -b.x2f(k,(js+j  ),i);  // reflect 2-field
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(js-j),i) =  b.x3f(k,(js+j-1),i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void ProjectPressureOuterX2()
//  \brief  Pressure is integated into ghost cells to improve hydrostatic eqm

void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      if (n==(IVY)) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(IVY,k,je+j,i) = -prim(IVY,k,je-j+1,i);  // reflect 2-velocity
        }
      } else if (n==(IPR)) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(IPR,k,je+j,i) = prim(IPR,k,je-j+1,i)
             + prim(IDN,k,je-j+1,i)*grav_acc*(2*j-1)*pco->dx2f(j);
        }
      } else {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,je+j,i) = prim(n,k,je-j+1,i);
        }
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones, reflecting b2
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(je+j  ),i) =  b.x1f(k,(je-j+1),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+j+1),i) = -b.x2f(k,(je-j+1),i);  // reflect 2-field
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(je+j  ),i) =  b.x3f(k,(je-j+1),i);
      }
    }}
  }

  return;
}

//----------------------------------------------------------------------------------------

// refinement condition: velocity gradient
int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real vgmax=0.0;
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real vgy=std::fabs(w(IVY,k,j,i+1)-w(IVY,k,j,i-1))*0.5;
        Real vgx=std::fabs(w(IVX,k,j+1,i)-w(IVX,k,j-1,i))*0.5;
        if (vgy > vgmax) vgmax=vgy;
        if (vgx > vgmax) vgmax=vgx;
      }
    }
  }
  if (vgmax > 0.01) return 1;
  return -1;
}
