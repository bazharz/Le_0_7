#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

// -----------------------------------------------------------
// Search for the closest index in an array to a given value
// using the bisection technique.
// INPUTS/OUTPUTS:
// xtable(0:n-1) => array to search in (ascending order)
// n             => array size
// x             => x location
// idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
// -----------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void 
pelelmex_initdata(
  int i,
  int j,
  int k,
  int is_incompressible,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& aux,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const * pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx      = geomdata.CellSize();

  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
               const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
	       const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);
  AMREX_D_TERM(
               const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
               const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
               const amrex::Real Lz = prob_hi[2] - prob_lo[2];)


  // Fill in the velocities


  int inSize = prob_parm.input_resolution;
  int Nx = prob_parm.input_nx;
  int Ny = prob_parm.input_ny;
  int Nz = prob_parm.input_nz;


  //
  // Fill Velocity
  //
  state(i,j,k,VELX) = 0.0;  //0.0;
  state(i,j,k,VELY) = 0.0;
  state(i,j,k,VELZ) = 0.0;

  if ( is_incompressible != 0) {
     return;
  }

  state(i,j,k,TEMP) = 300.0;  //prob_parm.d_Tinput[i + Nx * j + Nx * Ny * k];
  amrex::Print() << "i =" << i << "j =" << j << "k =" << k << "state(i,j,k,TEMP) =" << state(i,j,k,TEMP) << "\n";

  amrex::Real Yt[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  /*Yt[H_ID] = prob_parm.d_Y_H[i + Nx * j + Nx * Ny * k];  
  Yt[H2_ID] = prob_parm.d_Y_H2[i + Nx * j + Nx * Ny * k];
  Yt[O_ID] = prob_parm.d_Y_O[i + Nx * j + Nx * Ny * k];
  Yt[O2_ID] = prob_parm.d_Y_O2[i + Nx * j + Nx * Ny * k];
  Yt[OH_ID] = prob_parm.d_Y_OH[i + Nx * j + Nx * Ny * k];
  Yt[H2O_ID] = prob_parm.d_Y_H2O[i + Nx * j + Nx * Ny * k];
  Yt[HO2_ID] = prob_parm.d_Y_HO2[i + Nx * j + Nx * Ny * k];
  Yt[H2O2_ID] = prob_parm.d_Y_H2O2[i + Nx * j + Nx * Ny * k];
  Yt[N2_ID] = prob_parm.d_Y_N2[i + Nx * j + Nx * Ny * k];*/
  
  Yt[O2_ID] = 0.226;
  Yt[N2_ID] = 0.745;
  Yt[H2_ID] = 0.029;
  amrex::Print() << "Size of Yt =" << sizeof(Yt) << "\n";
  amrex::Print() << "Yt[H2] = " << Yt[H2_ID] << "\n";
 
  auto eos = pele::physics::PhysicsType::eos();   
  eos.X2Y(Yt, massfrac);
  amrex::Real rho_cgs, P_cgs;
  P_cgs = prob_parm.P_mean * 10.0;

  eos.PYT2R(P_cgs, Yt, state(i,j,k,TEMP), rho_cgs);  //Yt, state(i,j,k,TEMP), rho_cgs);
  state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

  eos.TY2H(state(i,j,k,TEMP), Yt, state(i,j,k,RHOH));  //Yt, state(i,j,k,RHOH));
  state(i,j,k,RHOH) *= 1.0e-4 * state(i,j,k,DENSITY);   // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i,j,k,FIRSTSPEC+n) = Yt[n] * state(i,j,k,DENSITY);  //Yt[n] * state(i,j,k,DENSITY);
  }
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}
#endif
