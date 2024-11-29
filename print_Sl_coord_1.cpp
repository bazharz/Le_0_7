#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BLFort.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_FabConv.H>

using namespace amrex;

// new function for print MultiFab

#include <vector>
#include <utility> // for std::pair
#include <algorithm> // for std::sort

// Structure de données pour stocker les coordonnées et les valeurs de chaque point
struct PointData {
    std::array<Real, AMREX_SPACEDIM> coords;
    Real value;
};


void printMultiFab(const std::vector<MultiFab>& mf_vec, const std::vector<Geometry>& geom_vec, const std::string& label)
{
    //Print() << label << " info:\n";

    // Vector pour stocker les données de chaque point de tous les niveaux de maillage
    std::vector<PointData> pointData;

    // Collecte de toutes les données de tous les niveaux de maillage
    for (int lev = 0; lev < mf_vec.size(); ++lev) {
        const MultiFab& mf = mf_vec[lev];
        const Geometry& geom = geom_vec[lev];

        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();

            // Récupérer les coordonnées de la boîte
            RealBox real_box(geom.ProbLo(), geom.ProbHi());
            std::vector<std::pair<Real, Real>> coord_limits;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                coord_limits.push_back(std::make_pair(real_box.lo(idim), real_box.hi(idim)));
            }

            // Parcourir tous les points de la boîte
            for (BoxIterator bit(bx); bit.ok(); ++bit) {
                const IntVect& iv = bit();
                PointData pd;

                // Calculer les coordonnées du point
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    Real dx = geom.CellSize(idim);
                    pd.coords[idim] = geom.ProbLo(idim) + (iv[idim] + 0.5) * dx;
                }

                // Accéder à la valeur de C au point
                pd.value = mf[mfi](iv, 0); // Supposant qu'il n'y a qu'une composante dans mf

                // Ajouter les données du point à la structure de données
                pointData.push_back(pd);
            }
        }
    }

    // Trier les données de tous les points
    std::sort(pointData.begin(), pointData.end(), [](const PointData& a, const PointData& b) {
        // Trier par coordonnées en x pour une valeur de y donnée
        if (a.coords[1] != b.coords[1]) {
            return a.coords[1] < b.coords[1]; // Tri par y d'abord
        } else {
            return a.coords[0] < b.coords[0]; // Si les coordonnées en y sont identiques, tri par x
        }
    });

    // Imprimer les données triées
    for (const auto& pd : pointData) {
        // Afficher les coordonnées et la valeur de C
        //Print() << "Coords: ";
        for (const auto& coord : pd.coords) {
            Print() << coord << " ";
        }
        Print() << pd.value << "\n";  //", C: " << pd.value << "\n";
    }
}

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> \n\tOptions:\n\tis_per=<L M N> T0=<T0> T1=<T1> \n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {
    if (argc < 2) {
      print_usage(argc,argv);
    }

    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    std::string infile        = "";  
    int finestLevel           = 1000;
    int nAuxVar               = 0;

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    if (pp.contains("help")) {
      print_usage(argc,argv);
    }

    pp.get("infile",infile);
    pp.query("finestLevel",finestLevel);

    // Initialize DataService
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Plotfile global infos
    finestLevel = std::min(finestLevel,amrData.FinestLevel());
    int Nlev = finestLevel + 1;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    RealBox rb(&(amrData.ProbLo()[0]), 
               &(amrData.ProbHi()[0]));

    // Gradient variable
    // Commented out because it's not used in this script
    // int idC = -1;

    // Auxiliary variables
    nAuxVar = pp.countval("Aux_Variables");
    Vector<std::string> AuxVar(nAuxVar);
    for(int ivar = 0; ivar < nAuxVar; ++ivar) { 
         pp.get("Aux_Variables", AuxVar[ivar],ivar);
    }

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    const int idCst = 0;
    int nCompIn = idCst + 1;
    Vector<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idCst]; // Using idCst instead of idC

    if (nAuxVar>0)
    {
        inVarNames.resize(nCompIn+nAuxVar);
        for (int ivar=0; ivar<nAuxVar; ++ivar) {
            if ( amrData.StateNumber(AuxVar[ivar]) < 0 ) {
               amrex::Abort("Unknown auxiliary variable name: "+AuxVar[ivar]);
            }
            inVarNames[nCompIn] = AuxVar[ivar];
            nCompIn++;
        } 
    }

    Vector<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i) {
      destFillComps[i] = i;
    }

    const int idGr = nCompIn;
    const int nCompOut = idGr + AMREX_SPACEDIM + 18 ; // 1 component stores the ||gradC||
    //Print() << "nCompOut =" << nCompOut << "\n";

    // Check symmetry/periodicity in given coordinate direction
    Vector<int> sym_dir(AMREX_SPACEDIM,0);
    pp.queryarr("sym_dir",sym_dir,0,AMREX_SPACEDIM);  

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    //Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        //Print() << is_per[idim] << " ";
    }
    //Print() << "\n";
    BCRec gradVarBC;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        gradVarBC.setLo(idim,BCType::foextrap);
        gradVarBC.setHi(idim,BCType::foextrap);
        if ( is_per[idim] ) {
            gradVarBC.setLo(idim, BCType::int_dir);
            gradVarBC.setHi(idim, BCType::int_dir);
        }
    }

    int coord = 0;

    // ---------------------------------------------------------------------
    // Let's start the real work
    // ---------------------------------------------------------------------
    Vector<MultiFab> state(Nlev);
    Vector<MultiFab> state1(Nlev);
    Vector<MultiFab> state2(Nlev);
    Vector<MultiFab> state3(Nlev);
    Vector<MultiFab> omegaC(Nlev);
    Vector<MultiFab> minus(Nlev);
    Vector<MultiFab> YC(Nlev);
    Vector<MultiFab> YC_min(Nlev);
    Vector<MultiFab> YC_max(Nlev);
    Vector<MultiFab> YC_eq(Nlev);
    Vector<MultiFab> YC_numerateur(Nlev);
    Vector<MultiFab> YC_denominateur(Nlev);
    Vector<MultiFab> C(Nlev);
    Vector<MultiFab> nC(Nlev);
    Vector<MultiFab> Sl_x(Nlev);
    Vector<MultiFab> Sl_y(Nlev);
    Vector<MultiFab> Sl_x_2(Nlev);
    Vector<MultiFab> Sl_y_2(Nlev);
    Vector<MultiFab> Sl(Nlev);
    Vector<MultiFab> D_N2(Nlev);
    Vector<MultiFab> D_H2(Nlev);
    Vector<MultiFab> D_O2(Nlev);
    Vector<MultiFab> D_H2O(Nlev);
    Vector<MultiFab> rho(Nlev);
    Vector<MultiFab> rho_zero(Nlev);
    Vector<MultiFab> VC(Nlev);
    Vector<MultiFab> rhoDgradYC_x(Nlev);
    Vector<MultiFab> rhoDgradYC_y(Nlev);
    Vector<MultiFab> d_x_rhoDgradYC_x(Nlev);
    Vector<MultiFab> d_y_rhoDgradYC_y(Nlev);
    Vector<MultiFab> rho_zerogradYC_x(Nlev);
    Vector<MultiFab> rho_zerogradYC_y(Nlev);
    Vector<MultiFab> gradYC_x(Nlev);
    Vector<MultiFab> gradYC_y(Nlev);
    Vector<MultiFab> maggradYC(Nlev);
    Vector<MultiFab> gradstate2_x(Nlev);
    Vector<MultiFab> gradstate2_y(Nlev);
    Vector<MultiFab> gradstate3_x(Nlev);
    Vector<MultiFab> gradstate3_y(Nlev);
    Vector<MultiFab> Temp(Nlev);
    Vector<MultiFab> Y_H2(Nlev);
    Vector<MultiFab> Y_O2(Nlev);
    Vector<MultiFab> Y_H2O(Nlev);
    Vector<MultiFab> Y_H2_eq(Nlev);
    Vector<MultiFab> Y_O2_eq(Nlev);
    Vector<MultiFab> Y_H2O_eq(Nlev);
    Vector<MultiFab> omegaH2(Nlev);
    Vector<MultiFab> omegaO2(Nlev);
    Vector<MultiFab> omegaH2O(Nlev);
    Vector<MultiFab> tableau_1(Nlev);
    Vector<MultiFab> omegaT(Nlev);
    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    const int nGrow = 1;

    // Read data on all the levels
    for (int lev=0; lev<Nlev; ++lev) {

      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));//,nullptr);
      state[lev].define(grids[lev], dmap[lev], nCompOut, nGrow); 
      state1[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      state2[lev].define(grids[lev], dmap[lev], nCompOut, nGrow); // rho * D * gradC_x --> compute gradient of this quantity
      state3[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      omegaC[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      minus[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      nC[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Sl_x[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Sl_y[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Sl_x_2[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Sl_y_2[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Sl[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Sl[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      D_N2[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      D_H2[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      D_O2[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      D_H2O[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      VC[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      rhoDgradYC_x[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      rhoDgradYC_y[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      rho_zerogradYC_x[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      rho_zerogradYC_y[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      d_x_rhoDgradYC_x[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      d_y_rhoDgradYC_y[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Temp[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Y_H2[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Y_O2[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Y_H2O[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Y_H2_eq[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Y_O2_eq[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      Y_H2O_eq[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      YC[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      YC_min[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      YC_max[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      YC_eq[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      YC_numerateur[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      YC_denominateur[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      C[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      omegaH2[lev].define(grids[lev], dmap[lev], 1, nGrow);
      omegaO2[lev].define(grids[lev], dmap[lev], 1, nGrow);
      omegaH2O[lev].define(grids[lev], dmap[lev], 1, nGrow);
      omegaT[lev].define(grids[lev], dmap[lev], 1, nGrow);
      rho[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      rho_zero[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      gradYC_x[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      gradYC_y[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      gradstate2_x[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      gradstate2_y[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      gradstate3_x[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      gradstate3_y[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      maggradYC[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      tableau_1[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);

      //Print() << "Reading data for level: " << lev << std::endl; 
      amrData.FillVar(Temp[lev], lev, "temp", 0); 
      amrData.FillVar(Y_H2[lev], lev, "Y(H2)", 0);
      amrData.FillVar(Y_O2[lev], lev, "Y(O2)", 0);
      amrData.FillVar(Y_H2O[lev], lev, "Y(H2O)", 0);
      amrData.FillVar(omegaH2[lev], lev, "I_R(H2)", 0);
      amrData.FillVar(omegaO2[lev], lev, "I_R(O2)", 0);
      amrData.FillVar(omegaH2O[lev], lev, "I_R(H2O)", 0);
      amrData.FillVar(omegaT[lev], lev, "HeatRelease", 0);
      
      // Call rho
      amrData.FillVar(rho[lev], lev, "density", 0);

      // Array -1.0
      minus[lev].setVal(-1.0, 0, 1, 1); 
      
      // Find the maximum of rho[lev]
      Real max_rho = rho[lev].max(0, 0, false);  // Max of component 0, no ghost cells
      // Aray rho_zero
      rho_zero[lev].setVal(max_rho, 0, 1, 1); 


      //compute Yc = YH2O - YH2 - YO2
      MultiFab::Multiply(Y_O2[lev], minus[lev], 0, 0, 1, 1);
      MultiFab::Multiply(Y_H2[lev], minus[lev], 0, 0, 1, 1);
      MultiFab::Add(Y_H2O[lev], Y_O2[lev], 0, 0, 1, 1);
      MultiFab::Add(Y_H2O[lev], Y_H2[lev], 0, 0, 1, 1);
      MultiFab::Copy(YC[lev], Y_H2O[lev], 0, 0, 1, 0);
      MultiFab::Copy(state[lev], YC[lev], 0, 0, 1, 0); // copy YC into state --> for computing the gradient of YC
      
      // Find the minimum of YC[lev]
      Real min_value = YC[lev].min(0, 0, false);  // Min of component 0, no ghost cells

      // Print the minimum value
      //Print() << "Minimum value of YC[lev] on level " << lev << ": " << min_value << std::endl;

      // Fill YC_min array by -min_value
      YC_min[lev].setVal(-min_value, 0, 1, 1);

      // Find the maximum of YC[lev]
      Real max_value = YC[lev].max(0, 0, false);  // Max of component 0, no ghost cells

      // Print the minimum value
      //Print() << "Maximum value of YC[lev] on level " << lev << ": " << max_value << std::endl;

      // Fill YC_min array by -min_value
      YC_max[lev].setVal(max_value, 0, 1, 1);
      
      // Find the minimum of Y_H2[lev]
      Real min_Y_H2 = Y_H2[lev].min(0, 0, false);  // Min of component 0, no ghost cells
      // Find the minimum of Y_O2[lev]
      Real min_Y_O2 = Y_O2[lev].min(0, 0, false);  // Min of component 0, no ghost cells
      // Find the maximum of Y_H2O[lev]
      Real max_Y_H2O = Y_H2O[lev].max(0, 0, false);  // Min of component 0, no ghost cells

      // compute Yc_eq = YH2O_eq - YH2_eq - YO2_eq

      Y_H2_eq[lev].setVal(-min_Y_H2, 0, 1, 1); 
      Y_O2_eq[lev].setVal(-min_Y_O2, 0, 1, 1); 
      Y_H2O_eq[lev].setVal(max_Y_H2O, 0, 1, 1);
      MultiFab::Add(Y_H2O_eq[lev], Y_O2_eq[lev], 0, 0, 1, 1);
      MultiFab::Add(Y_H2O_eq[lev], Y_H2_eq[lev], 0, 0, 1, 1);
      MultiFab::Copy(YC_eq[lev], Y_H2O_eq[lev], 0, 0, 1, 0);

      //compute C = (Yc - Yc,min)/(Yc,max - Yc,min)
      MultiFab::Copy(YC_denominateur[lev], YC_max[lev], 0, 0, 1, 0); // Yc_denominateur
      MultiFab::Add(YC_denominateur[lev], YC_min[lev], 0, 0, 1, 1); // Yc_denominateur = Yc,max - Yc,min
      MultiFab::Copy(YC_numerateur[lev], YC[lev], 0, 0, 1, 0); // Yc_numerateur
      MultiFab::Add(YC_numerateur[lev], YC_min[lev], 0, 0, 1, 1); // Yc_numerateur = Yc - Yc,min
      MultiFab::Divide(YC_numerateur[lev], YC_denominateur[lev], 0, 0, 1, 1);
      MultiFab::Copy(C[lev], YC_numerateur[lev], 0, 0, 1, 0);

      // Find the maximum of C[lev]
      Real max_C = C[lev].max(0, 0, false);  // Max of component 0, no ghost cells

      // Print the minimum value
      //Print() << "Maximum value of C[lev] on level " << lev << ": " << max_C << std::endl;

      // Find the minimum of C[lev]
      Real min_C = C[lev].min(0, 0, false);  // Min of component 0, no ghost cells

      // Print the minimum value
      //Print() << "Minimum value of C[lev] on level " << lev << ": " << min_C << std::endl;

      //printMultiFab(C[lev], "C[lev]", lev);

      // compute omegaC = omegaH2O - omegaH2 - omegaO2
      MultiFab::Multiply(omegaO2[lev], minus[lev], 0, 0, 1, 1);
      MultiFab::Multiply(omegaH2[lev], minus[lev], 0, 0, 1, 1);
      MultiFab::Add(omegaH2O[lev], omegaH2[lev], 0, 0, 1, 1);
      MultiFab::Add(omegaH2O[lev], omegaO2[lev], 0, 0, 1, 1);
      MultiFab::Copy(omegaC[lev], omegaH2O[lev], 0, 0, 1, 0);
      //printMultiFab(omegaC[lev], "omegaC[lev]", lev);

      // Diffusion coefficient : N2
      amrData.FillVar(D_N2[lev], lev, "D_N2", 0);
      amrData.FillVar(D_H2[lev], lev, "D_H2", 0);
      amrData.FillVar(D_O2[lev], lev, "D_O2", 0);
      amrData.FillVar(D_H2O[lev], lev, "D_H2O", 0);
      //printMultiFab(D_N2[lev], "D_N2[lev]", lev);
      
      // Defintion of rhoDgradYC_x & rhoDgradYC_y
      MultiFab::Copy(rhoDgradYC_x[lev],D_N2[lev], 0, 0, 1, 0);
      MultiFab::Multiply(rhoDgradYC_x[lev], rho[lev], 0, 0, 1, 1);
      MultiFab::Copy(rhoDgradYC_y[lev],rhoDgradYC_x[lev], 0, 0, 1, 0);


      state[lev].FillBoundary(idCst,1,geoms[lev].periodicity());
    }

    // Get face-centered gradients from MLMG
    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);
    MLPoisson poisson({geoms}, {grids}, {dmap}, info);
    poisson.setMaxOrder(4);
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
    for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
       if (is_per[idim] == 1) {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
       } else {
          if (sym_dir[idim] == 1) {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::reflect_odd;
          } else {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
          }
       }
    }
    poisson.setDomainBC(lo_bc, hi_bc);

    // Need to apply the operator to ensure CF consistency with composite solve
    int nGrowGrad = 0;                   // No need for ghost face on gradient
    Vector<Array<MultiFab,AMREX_SPACEDIM>> gradYC(Nlev);
    Vector<std::unique_ptr<MultiFab>> phi;
    Vector<MultiFab> laps;
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         const auto& ba = grids[lev];
         gradYC[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                dmap[lev], 1, nGrowGrad);
      }    
      phi.push_back(std::make_unique<MultiFab> (state[lev],amrex::make_alias,idCst,1));
      poisson.setLevelBC(lev, phi[lev].get());
      laps.emplace_back(grids[lev], dmap[lev], 1, 1);
    }

    MLMG mlmg(poisson);
    mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    mlmg.getFluxes(GetVecOfArrOfPtrs(gradYC), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);

    for (int lev = 0; lev < Nlev; ++lev) {
        // Convert to cell avg gradient
        MultiFab gradAlias(state[lev], amrex::make_alias, idGr, AMREX_SPACEDIM);
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(gradYC[lev]));
        gradAlias.mult(-1.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(state[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {    
           const Box& bx = mfi.tilebox();
           auto const& grad_a   = gradAlias.const_array(mfi);
           auto const& gradMag  = state[lev].array(mfi,idGr+AMREX_SPACEDIM);

           amrex::ParallelFor(bx, [=]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {    
              gradMag(i,j,k) = std::sqrt(AMREX_D_TERM(grad_a(i,j,k,0) * grad_a(i,j,k,0),
                                                      + grad_a(i,j,k,1) * grad_a(i,j,k,1),
                                                      + grad_a(i,j,k,2) * grad_a(i,j,k,2)));
           });
	   const Box& gradMagBox = mfi.validbox();
        }

	//Copy of gradient of C and its magnitude in independant arrays
	const MultiFab& gradC_x_in_state = gradYC[lev][0];
	const MultiFab& gradC_y_in_state = gradYC[lev][1];
	const MultiFab& maggradC_in_state = state[lev];

	MultiFab::Copy(gradYC_x[lev], gradC_x_in_state, 0, 0, 1, nGrow);
	MultiFab::Copy(gradYC_y[lev], gradC_y_in_state, 0, 0, 1, nGrow);
	MultiFab::Copy(maggradYC[lev], maggradC_in_state, 0, 0, 1, nGrow);
	// Computation of nc = gradYC_x / ||gradC||
	MultiFab::Copy(nC[lev], gradYC_x[lev], 0, 0, 1, nGrow);
	MultiFab::Divide(nC[lev], maggradYC[lev], 0, 0, 1, 1);
    }

    // -----------------------------------------------------------------
    // Compute rhoDgradYC_x & rho_zero_gradYC_x
    // -----------------------------------------------------------------

    for (int lev =0; lev < Nlev; ++lev) {
      
      // Computation of rhoDgradYC_x = rho * D_N2 * gradYC_x, already rhoDgradYC_x = rho * D_N2
      MultiFab::Multiply(rhoDgradYC_x[lev], gradYC_x[lev], 0, 0, 1, 1);
      MultiFab::Copy(state2[lev], rhoDgradYC_x[lev], 0, 0, 1, 0); // copy into state2 for computing gradient

      // Computing rho_zero * gradYC_x
      MultiFab::Copy(rho_zerogradYC_x[lev], rho_zero[lev], 0, 0, 1, 0); // rho_zerogradYC = rho_zero
      MultiFab::Multiply(rho_zerogradYC_x[lev], gradYC_x[lev], 0, 0, 1, 1); // rho_zerogradYC = rho_zero * gradYC_x

      state2[lev].FillBoundary(idCst,1,geoms[lev].periodicity()); // Fill boundary periodicity
      
    }


    // -----------------------------------------------------------------
    // Compute the gradient(rho * D * grad(YC)_x)
    // -----------------------------------------------------------------

    // Get face-centered gradients from MLMG
    LPInfo info1;
    info1.setAgglomeration(1);
    info1.setConsolidation(1);
    info1.setMetricTerm(false);
    info1.setMaxCoarseningLevel(0);
    MLPoisson poisson1({geoms}, {grids}, {dmap}, info1);
    poisson1.setMaxOrder(4);
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc1;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc1;
    for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
       if (is_per[idim] == 1) {
          lo_bc1[idim] = hi_bc1[idim] = LinOpBCType::Periodic;
       } else {
          if (sym_dir[idim] == 1) {
             lo_bc1[idim] = hi_bc1[idim] = LinOpBCType::reflect_odd;
          } else {
             lo_bc1[idim] = hi_bc1[idim] = LinOpBCType::Neumann;
          }
       }
    }
    poisson1.setDomainBC(lo_bc1, hi_bc1);

    // Need to apply the operator to ensure CF consistency with composite solve
    int nGrowGrad1 = 0;                   // No need for ghost face on gradient
    Vector<Array<MultiFab,AMREX_SPACEDIM>> gradstate2(Nlev);
    Vector<std::unique_ptr<MultiFab>> phi1;
    Vector<MultiFab> laps1;
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         const auto& ba = grids[lev];
         gradstate2[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                dmap[lev], 1, nGrowGrad1);
      }    
      phi1.push_back(std::make_unique<MultiFab> (state2[lev],amrex::make_alias,idCst,1));
      poisson1.setLevelBC(lev, phi1[lev].get());
      laps1.emplace_back(grids[lev], dmap[lev], 1, 1);
    }

    MLMG mlmg1(poisson1);
    mlmg1.apply(GetVecOfPtrs(laps1), GetVecOfPtrs(phi1));
    mlmg1.getFluxes(GetVecOfArrOfPtrs(gradstate2), GetVecOfPtrs(phi1), MLMG::Location::FaceCenter);

    for (int lev = 0; lev < Nlev; ++lev) {
        // Convert to cell avg gradient
        MultiFab gradAlias1(state2[lev], amrex::make_alias, idGr, AMREX_SPACEDIM);
        average_face_to_cellcenter(gradAlias1, 0, GetArrOfConstPtrs(gradstate2[lev]));
        gradAlias1.mult(-1.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(state2[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {    
           const Box& bx = mfi.tilebox();
           auto const& grad_a1   = gradAlias1.const_array(mfi);
           auto const& gradMag1  = state2[lev].array(mfi,idGr+AMREX_SPACEDIM);

           amrex::ParallelFor(bx, [=]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {    
              gradMag1(i,j,k) = std::sqrt(AMREX_D_TERM(grad_a1(i,j,k,0) * grad_a1(i,j,k,0),
                                                      + grad_a1(i,j,k,1) * grad_a1(i,j,k,1),
                                                      + grad_a1(i,j,k,2) * grad_a1(i,j,k,2)));
           });
	   const Box& gradMagBox1 = mfi.validbox();
        }

	//Copy of gradient terms of C and its magnitude in independant arrays
	const MultiFab& gradstate2_x_in_state2 = gradstate2[lev][0];
	const MultiFab& gradstate2_y_in_state2 = gradstate2[lev][1];
	const MultiFab& maggradstate2_in_state2 = state2[lev];

	MultiFab::Copy(gradstate2_x[lev], gradstate2_x_in_state2, 0, 0, 1, nGrow);
	MultiFab::Copy(gradstate2_y[lev], gradstate2_y_in_state2, 0, 0, 1, nGrow);
    }
    

    // ---------------------------------------------------------------------
    //  Compute Sl_x
    // ---------------------------------------------------------------------
    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab::Copy(Sl_x[lev], gradstate2_x[lev], 0, 0, 1, 0); //  copy of d_x_rhoDgradYC_x
      MultiFab::Multiply(Sl_x[lev], minus[lev], 0, 0, 1, 1); // Sl = - grad_x(rho * D * grad(YC))
      MultiFab::Add(Sl_x[lev], omegaC[lev], 0, 0, 1, 1); // Sl = -grad_x(rho * D * grad(YC)) + omegaC
      MultiFab::Divide(Sl_x[lev], rho_zerogradYC_x[lev], 0, 0, 1, 1); // Sl = (-grad_x(rho * D * grad(YC)) + omegaC) / (rho_zero * grad_x(Yc))
    }

    // Après le remplissage du tableau Sl_x, ajoutez la condition suivante :
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int i = 0; i < Sl_x[lev].boxArray().size(); ++i) {
        const Box& bx = Sl_x[lev][i].box();
        const auto& sl_arr = Sl_x[lev][i].array();
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (std::isinf(sl_arr(i, j, k, n))) {
               sl_arr(i, j, k, n) = 0.0;  // Replace -inf by 0
               }
            });
        }
      }
    }


    // ---------------------------------------------------------------------
    //  print Sl_x
    // ---------------------------------------------------------------------
    //for (int lev = 0; lev < Nlev; ++lev) {
      //printMultiFab(Sl_x[lev], "Sl_x[lev]", lev);
      //printMultiFab(C[lev], "C[lev]", lev);
      //printMultiFab(C[lev], geoms[lev], "C[lev]", lev);
    //}
    //printMultiFab(C, geoms, "C");
    //printMultiFab(Sl_x, geoms, "Sl_x");


    // -----------------------------------------------------------------
    // Compute rhoDgradYC_y & rho_zero_gradYC_y
    // -----------------------------------------------------------------

    for (int lev =0; lev < Nlev; ++lev) {
      
      // Computation of rhoDgradYC_y = rho * D_N2 * gradYC_y, already rhoDgradYC_y = rho * D_N2
      MultiFab::Multiply(rhoDgradYC_y[lev], gradYC_y[lev], 0, 0, 1, 1);
      MultiFab::Copy(state3[lev], rhoDgradYC_y[lev], 0, 0, 1, 0); // copy into state3 for computing gradient

      // Computing rho_zero * gradYC_y
      MultiFab::Copy(rho_zerogradYC_y[lev], rho_zero[lev], 0, 0, 1, 0); // rho_zerogradYC = rho_zero
      MultiFab::Multiply(rho_zerogradYC_y[lev], gradYC_y[lev], 0, 0, 1, 1); // rho_zerogradYC = rho_zero * gradYC_y

      state3[lev].FillBoundary(idCst,1,geoms[lev].periodicity()); // Fill boundary periodicity
      
    }


    // -----------------------------------------------------------------
    // Compute the gradient(rho * D * grad(YC)_y)
    // -----------------------------------------------------------------

    // Get face-centered gradients from MLMG
    LPInfo info2;
    info2.setAgglomeration(1);
    info2.setConsolidation(1);
    info2.setMetricTerm(false);
    info2.setMaxCoarseningLevel(0);
    MLPoisson poisson2({geoms}, {grids}, {dmap}, info2);
    poisson2.setMaxOrder(4);
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc2;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc2;
    for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
       if (is_per[idim] == 1) {
          lo_bc2[idim] = hi_bc2[idim] = LinOpBCType::Periodic;
       } else {
          if (sym_dir[idim] == 1) {
             lo_bc2[idim] = hi_bc2[idim] = LinOpBCType::reflect_odd;
          } else {
             lo_bc2[idim] = hi_bc2[idim] = LinOpBCType::Neumann;
          }
       }
    }
    poisson2.setDomainBC(lo_bc2, hi_bc2);

    // Need to apply the operator to ensure CF consistency with composite solve
    int nGrowGrad2 = 0;                   // No need for ghost face on gradient
    Vector<Array<MultiFab,AMREX_SPACEDIM>> gradstate3(Nlev);
    Vector<std::unique_ptr<MultiFab>> phi2;
    Vector<MultiFab> laps2;
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         const auto& ba = grids[lev];
         gradstate3[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                dmap[lev], 1, nGrowGrad2);
      }    
      phi2.push_back(std::make_unique<MultiFab> (state3[lev],amrex::make_alias,idCst,1));
      poisson2.setLevelBC(lev, phi2[lev].get());
      laps2.emplace_back(grids[lev], dmap[lev], 1, 1);
    }

    MLMG mlmg2(poisson2);
    mlmg2.apply(GetVecOfPtrs(laps2), GetVecOfPtrs(phi2));
    mlmg2.getFluxes(GetVecOfArrOfPtrs(gradstate3), GetVecOfPtrs(phi2), MLMG::Location::FaceCenter);

    for (int lev = 0; lev < Nlev; ++lev) {
        // Convert to cell avg gradient
        MultiFab gradAlias2(state3[lev], amrex::make_alias, idGr, AMREX_SPACEDIM);
        average_face_to_cellcenter(gradAlias2, 0, GetArrOfConstPtrs(gradstate3[lev]));
        gradAlias2.mult(-1.0);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(state3[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {    
           const Box& bx = mfi.tilebox();
           auto const& grad_a2   = gradAlias2.const_array(mfi);
           auto const& gradMag2  = state3[lev].array(mfi,idGr+AMREX_SPACEDIM);

           amrex::ParallelFor(bx, [=]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {    
              gradMag2(i,j,k) = std::sqrt(AMREX_D_TERM(grad_a2(i,j,k,0) * grad_a2(i,j,k,0),
                                                      + grad_a2(i,j,k,1) * grad_a2(i,j,k,1),
                                                      + grad_a2(i,j,k,2) * grad_a2(i,j,k,2)));
           });
	   const Box& gradMagBox2 = mfi.validbox();
        }

	//Copy of gradient terms of C and its magnitude in independant arrays
	const MultiFab& gradstate3_x_in_state3 = gradstate3[lev][0];
	const MultiFab& gradstate3_y_in_state3 = gradstate3[lev][1];
	const MultiFab& maggradstate3_in_state3 = state3[lev];

	MultiFab::Copy(gradstate3_x[lev], gradstate3_x_in_state3, 0, 0, 1, nGrow);
	MultiFab::Copy(gradstate3_y[lev], gradstate3_y_in_state3, 0, 0, 1, nGrow);
    }
    

    // ---------------------------------------------------------------------
    //  Compute Sl_y
    // ---------------------------------------------------------------------
    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab::Copy(Sl_y[lev], gradstate3_y[lev], 0, 0, 1, 0); //  copy of d_y_rhoDgradYC_y
      MultiFab::Multiply(Sl_y[lev], minus[lev], 0, 0, 1, 1); // Sl_y = - grad_y(rho * D * grad(YC))
      MultiFab::Add(Sl_y[lev], omegaC[lev], 0, 0, 1, 1); // Sl_y = -grad_xyrho * D * grad(YC)) + omegaC
      MultiFab::Divide(Sl_y[lev], rho_zerogradYC_y[lev], 0, 0, 1, 1); // Sl_y = (-grad_y(rho * D * grad(YC)) + omegaC) / (rho_zero * grad_y(Yc))
    }

    // Après le remplissage du tableau Sl_x, ajoutez la condition suivante :
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int i = 0; i < Sl_y[lev].boxArray().size(); ++i) {
        const Box& bx = Sl_y[lev][i].box();
        const auto& sl_arry = Sl_y[lev][i].array();
        for (int n = 0; n < AMREX_SPACEDIM; ++n) {
          amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (std::isinf(sl_arry(i, j, k, n))) {
               sl_arry(i, j, k, n) = 0.0;  // Replace -inf by 0
               }
            });
        }
      }
    }

    // Print Sl_y
    //printMultiFab(Sl_y, geoms, "Sl_x");



    // ---------------------------------------------------------------------
    // Compute Sl = sqrt_root(Sl_x**2 + Sl_y**2)
    // ---------------------------------------------------------------------
    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab::Copy(Sl_x_2[lev], Sl_x[lev], 0, 0, 1, 0); //  copy of Sl_x
      MultiFab::Copy(Sl_y_2[lev], Sl_y[lev], 0, 0, 1, 0); //  copy of Sl_y
      MultiFab::Multiply(Sl_x_2[lev], Sl_x[lev], 0, 0, 1, 1); // Sl_x_2 = Sl_x * Sl_x
      MultiFab::Multiply(Sl_y_2[lev], Sl_y[lev], 0, 0, 1, 1); // Sl_y_2 = Sl_y * Sl_y
      MultiFab::Copy(Sl[lev], Sl_x_2[lev], 0, 0, 1, 0); //  copy of Sl_x_2
      MultiFab::Add(Sl[lev], Sl_y_2[lev], 0, 0, 1, 1); // Sl**2 = Sl_x_2 + Sl_y_2
      // Compute the square root
      for (MFIter mfi(Sl[lev]); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        auto const& Sl_array = Sl[lev].array(mfi);

        amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Sl_array(i, j, k) = std::sqrt(Sl_array(i, j, k));
        });
      }

    }



    int indexC_05 = -1;  // Initialisation à une valeur non valide
    int lev_with_index = -1;  // Niveau où l'index a été trouvé
    for (int lev = 0; lev < Nlev; ++lev) {
      // Boucle sur les données pour trouver l'index où C est égale à 0.5
      for (int i = 0; i < C[lev].size(); ++i) {
	//double C_value = C[lev][i];
	//const MultiFab& current_C = C[lev][i];
	const FArrayBox& current_C = C[lev][i];
	if (std::abs(current_C.dataPtr()[0] - 0.5) < 8e-1) {
        //if (std::abs(current_C[0] - 0.5) < 1e-3) {  // tol est la tolérance pour la comparaison
            indexC_05 = i;
	    lev_with_index = lev;
	    //std::cout << "indexC_05 = " << indexC_05 << "lev_with_index =" << lev_with_index << std::endl;
            break;
        }
      }
      // Si l'index a été trouvé, sortir de la boucle externe
      //if (indexC_05 != -1) {
      //  break;
      //}
      if (indexC_05 != -1) {
        // Utiliser l'index trouvé pour obtenir la valeur de Sl à cet endroit
	double Sl_at_C_05 = Sl[lev_with_index][0].dataPtr()[indexC_05];  // Utilisez dataPtr() pour accéder aux données brutes
        //double Sl_at_C_05 = Sl[lev_with_index][indexC_05];
	//std::cout << "Sl[C=0.5] = " << Sl_at_C_05 << std::endl;
	//Print() << "Sl_at_C_05 =" << Sl_at_C_05 << "\n";
        // Maintenant, vous avez la valeur de Sl lorsque C est égale à 0.5
	} else {
             // Gérer le cas où l'index n'a pas été trouvé
             std::cerr << "Erreur : Aucun index trouvé pour C = 0.5." << std::endl;
	}
    }



    // ---------------------------------------------------------------------
    //  print Sl
    // ---------------------------------------------------------------------
    //for (int lev = 0; lev < Nlev; ++lev) {
    //  printMultiFab(Sl[lev], "Sl[lev]", lev);
    //}


    for (int lev = 0; lev < Nlev; ++lev) {

      MultiFab::Copy(state[lev], YC[lev], 0, idGr + AMREX_SPACEDIM + 1, 1, 0); // Copy YC into state
      MultiFab::Copy(state[lev], YC_eq[lev], 0, idGr + AMREX_SPACEDIM + 2, 1, 0); // Copy YC_eq into state
      MultiFab::Copy(state[lev], omegaC[lev], 0, idGr + AMREX_SPACEDIM + 3, 1, 0); // Copy omegaC into state
      MultiFab::Copy(state[lev], VC[lev], 0, idGr + AMREX_SPACEDIM + 4, 1, 0); // Copy VC into state
      MultiFab::Copy(state[lev], rho[lev], 0, idGr + AMREX_SPACEDIM + 5, 1, 0); // Copy rho into state
      MultiFab::Copy(state[lev], Temp[lev], 0, idGr + AMREX_SPACEDIM + 6, 1, 0); // Copy temp into state
      MultiFab::Copy(state[lev], D_N2[lev], 0, idGr + AMREX_SPACEDIM + 7, 1, 0); // Copy D_N2 into state
      MultiFab::Copy(state[lev], nC[lev], 0, idGr + AMREX_SPACEDIM + 8, 1, 0); // Copy nC into state
      MultiFab::Copy(state[lev], gradYC_x[lev], 0, idGr + AMREX_SPACEDIM + 9, 1, 0); // Copy gradYC into state
      MultiFab::Copy(state[lev], rho_zero[lev], 0, idGr + AMREX_SPACEDIM + 10, 1, 0); // Copy rho_zero into state
      MultiFab::Copy(state[lev], rho_zerogradYC_x[lev], 0, idGr + AMREX_SPACEDIM + 11, 1, 0); // Copy rho_zerogradYC into state
      MultiFab::Copy(state[lev], gradstate2_x[lev], 0, idGr + AMREX_SPACEDIM + 12, 1, 0); // Copy gradstate2_x into state
      MultiFab::Copy(state[lev], Sl_x[lev], 0, idGr + AMREX_SPACEDIM + 13, 1, 0); // Copy Sl_x into state
      MultiFab::Copy(state[lev], Sl_y[lev], 0, idGr + AMREX_SPACEDIM + 14, 1, 0); // Copy Sl_y into state
      MultiFab::Copy(state[lev], Sl[lev], 0, idGr + AMREX_SPACEDIM + 15, 1, 0); // Copy Sl_y into state
      MultiFab::Copy(state[lev], C[lev], 0, idGr + AMREX_SPACEDIM + 16, 1, 0); // Copy C into state
      MultiFab::Copy(state[lev], omegaT[lev], 0, idGr + AMREX_SPACEDIM + 17, 1, 0); // Copy C into state
    }


    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------
    Vector<std::string> nnames(nCompOut);
    for (int i=0; i<nCompIn; ++i) {
      nnames[i] = inVarNames[i];
    }
    nnames[idCst] = "YC";
    nnames[idGr+0] = "YC_gx";
    nnames[idGr+1] = "YC_gy";
#if AMREX_SPACEDIM==3
    nnames[idGr+2] = "YC_gz";
#endif
    nnames[idGr+AMREX_SPACEDIM] = "||gradYC||";

    // Additional variables
    nnames[idGr + AMREX_SPACEDIM + 1] = "YC";
    nnames[idGr + AMREX_SPACEDIM + 2] = "YC_eq";
    nnames[idGr + AMREX_SPACEDIM + 3] = "omegaC";
    nnames[idGr + AMREX_SPACEDIM + 4] = "VC";
    nnames[idGr + AMREX_SPACEDIM + 5] = "rho";
    nnames[idGr + AMREX_SPACEDIM + 6] = "Temp";
    nnames[idGr + AMREX_SPACEDIM + 7] = "D_N2";
    nnames[idGr + AMREX_SPACEDIM + 8] = "nC";
    nnames[idGr + AMREX_SPACEDIM + 9] = "gradYC_x";
    nnames[idGr + AMREX_SPACEDIM + 10] = "rho_zero";
    nnames[idGr + AMREX_SPACEDIM + 11] = "rho_zerogradYC";
    nnames[idGr + AMREX_SPACEDIM + 12] = "gradstate2_x";
    nnames[idGr + AMREX_SPACEDIM + 13] = "Sl_x";
    nnames[idGr + AMREX_SPACEDIM + 14] = "Sl_y";
    nnames[idGr + AMREX_SPACEDIM + 15] = "Sl";
    nnames[idGr + AMREX_SPACEDIM + 16] = "C";
    nnames[idGr + AMREX_SPACEDIM + 17] = "omegaT";

    std::string outfile(getFileRoot(infile)); // + "Sl_mag"); //pp.query("outfile",outfile);

    //Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(state), nnames,
                                   geoms, 0.0, isteps, refRatios);

  }
  amrex::Finalize();
  return 0;
}

