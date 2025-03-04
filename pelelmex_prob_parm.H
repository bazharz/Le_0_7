#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean = 101325.0;
   amrex::Real T_mean = 300.0;

   int input_resolution = 0;
   int input_nx = 0;
   int input_ny = 0;
   int input_nz = 0;
   amrex::Real urms0    = 1.0;
   amrex::Real uin_norm = 1.0;

   amrex::Real* d_Tinput = nullptr;
   amrex::Real* d_Y_H = nullptr;
   amrex::Real* d_Y_H2 = nullptr;
   amrex::Real* d_Y_O = nullptr;
   amrex::Real* d_Y_O2 = nullptr;
   amrex::Real* d_Y_OH = nullptr;
   amrex::Real* d_Y_H2O = nullptr;
   amrex::Real* d_Y_HO2 = nullptr;
   amrex::Real* d_Y_H2O2 = nullptr;
   amrex::Real* d_Y_N2 = nullptr;
};
#endif
