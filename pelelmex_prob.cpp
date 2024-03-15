#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

// -----------------------------------------------------------
// Read a binary file
// INPUTS/OUTPUTS:
// iname => filename
// nx    => input resolution
// ny    => input resolution
// nz    => input resolution
// data  <= output data
// -----------------------------------------------------------
void
read_binary(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  const size_t ncol,
  amrex::Vector<double>& data /*needs to be double*/)
{
  std::ifstream infile(iname, std::ios::in | std::ios::binary);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }

  for (size_t i = 0; i < nx * ny * nz * ncol; i++) {
    infile.read(reinterpret_cast<char*>(&data[i]), sizeof(data[i]));
  }
  infile.close();
}

AMREX_FORCE_INLINE
std::string
read_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

// -----------------------------------------------------------
// Read a csv file
// INPUTS/OUTPUTS:
// iname => filename
// nx    => input resolution
// ny    => input resolution
// nz    => input resolution
// data  <= output data
// -----------------------------------------------------------
void
read_csv(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<amrex::Real>& data)
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile);
  if (not infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  std::istringstream iss(memfile);

  // Read the file
  size_t nlines = 0;
  std::string firstline;
  std::string line;
  std::getline(iss, firstline); // skip header
  while (getline(iss, line)) {
    ++nlines;
  }

  // Quick sanity check
  if (nlines != nx * ny * nz) {
    amrex::Abort(
      "Number of lines in the input file (= " + std::to_string(nlines) +
      ") does not match the input resolution (=" + std::to_string(nx) + ")");
  }

  // Read the data from the file
  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline); // skip header
  int cnt = 0;
  while (std::getline(iss, line)) {
    std::istringstream linestream(line);
    std::string value;
    while (getline(linestream, value, ',')) {
      std::istringstream sinput(value);
      sinput >> data[cnt];
      cnt++;
    }
  }
}

void PeleLM::readProbParm()
{
    amrex::ParmParse pp("prob");

  ProbParm local_prob_parm;

  // Read HIT user-defined options from the input file
  pp.query("input_resolution", local_prob_parm.input_resolution);
    pp.query("input_nx",local_prob_parm.input_nx);
    pp.query("input_ny",local_prob_parm.input_ny);
    pp.query("input_nz",local_prob_parm.input_nz);
    if (local_prob_parm.input_resolution == 0 & local_prob_parm.input_nx == 0 & local_prob_parm.input_ny == 0 & local_prob_parm.input_nz == 0){
    amrex::Abort("for HIT, the input resolution cannot be 0 !");
   }

  std::string datafile;
  pp.query("input_name", datafile);
  int binfmt = 0; // Default is ASCII format
  pp.query("input_binaryformat", binfmt);
  pp.query("urms0", local_prob_parm.urms0);

  amrex::Real lambda0 = 0.5;
  amrex::Real tau = lambda0 / local_prob_parm.urms0;

  // Output initial conditions
  std::ofstream ofs("initialConditions.txt", std::ofstream::out);
  amrex::Print(ofs) << "lambda0, urms0, tau \n";
    amrex::Print(ofs) << "Hello \n";
  amrex::Print(ofs).SetPrecision(17) << lambda0 << "," << local_prob_parm.urms0 << "," << tau << std::endl;
  ofs.close();

  // Read initial velocity field
  const size_t nx = local_prob_parm.input_nx;  //prob_parm.input_resolution;
  const size_t ny = local_prob_parm.input_ny;  //prob_parm.input_resolution;
  const size_t nz = local_prob_parm.input_nz;  //prob_parm.input_resolution;
  amrex::Print() << "nx =" << nx << "ny =" << ny << "nz =" << nz << "\n";
  amrex::Vector<amrex::Real> data(nx * ny * nz * 10); /* this needs to be double */
  if (binfmt) {
    read_binary(datafile, nx, ny, nz, 10, data);
  } else {
    read_csv(datafile, nx, ny, nz, data);
  }

  // Extract position and velocities
  amrex::Vector<amrex::Real> Tinput(nx * ny * nz);
  amrex::Vector<amrex::Real> Y_Hinput(nx * ny * nz);
  amrex::Vector<amrex::Real> Y_H2input(nx * ny * nz);
  amrex::Vector<amrex::Real> Y_Oinput(nx * ny * nz);
  amrex::Vector<amrex::Real> Y_O2input(nx * ny * nz);
  amrex::Vector<amrex::Real> Y_OHinput(nx * ny * nz);
  amrex::Vector<amrex::Real> Y_H2Oinput(nx * ny * nz);
  amrex::Vector<amrex::Real> Y_HO2input(nx * ny * nz);
  amrex::Vector<amrex::Real> Y_H2O2input(nx * ny * nz);
  amrex::Vector<amrex::Real> Y_N2input(nx * ny * nz);

  //xinput.resize(nx * ny * nz);
  /*Tinput.resize(nx * ny * nz);
  Y_Hinput.resize(nx * ny * nz);
  Y_H2input.resize(nx * ny * nz);
  Y_Oinput.resize(nx * ny * nz);
  Y_O2input.resize(nx * ny * nz);
  Y_OHinput.resize(nx * ny * nz);
  Y_H2Oinput.resize(nx * ny * nz);
  Y_HO2input.resize(nx * ny * nz);
  Y_H2O2input.resize(nx * ny * nz);
  Y_N2input.resize(nx * ny * nz);*/

  for (long i = 0; i < nx * ny * nz; i++) {
    //xinput[i] = data[0 + i * 4];
    Tinput[i] = data [0 + i * 10];
    Y_Hinput[i] = data[1 + i * 10];
    Y_H2input[i] = data[2 + i * 10];
    Y_Oinput[i] = data[3 + i * 10];
    Y_O2input[i] = data[4 + i * 10];
    Y_OHinput[i] = data[5 + i * 10];
    Y_H2Oinput[i] = data[6 + i * 10];
    Y_HO2input[i] = data[7 + i * 10];
    Y_H2O2input[i] = data[8 + i * 10];
    Y_N2input[i] = data[9 + i * 10];
    amrex::Print() << "i =" << i << "Tinput =" << Tinput[i] << "H =" << Y_Hinput[i] <<"H2 =" << Y_H2input[i] << "N2 =" << Y_N2input[i] << "\n";
  }
    
    local_prob_parm.d_Tinput = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
    local_prob_parm.d_Y_H = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
    local_prob_parm.d_Y_H2 = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
    local_prob_parm.d_Y_O = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
    local_prob_parm.d_Y_O2 = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
    local_prob_parm.d_Y_OH = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
    local_prob_parm.d_Y_H2O = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
    local_prob_parm.d_Y_HO2 = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
    local_prob_parm.d_Y_H2O2 = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));
    local_prob_parm.d_Y_N2 = (amrex::Real*) amrex::The_Arena()->alloc(nx*ny*nz*sizeof(amrex::Real));

    for (int i = 0; i < nx*ny*nz; i++) {

       local_prob_parm.d_Tinput[i] = Tinput[i];
       local_prob_parm.d_Y_H[i] = Y_Hinput[i];
       local_prob_parm.d_Y_H2[i] = Y_H2input[i];
       local_prob_parm.d_Y_O[i] = Y_Oinput[i];
       local_prob_parm.d_Y_O2[i] = Y_O2input[i];
       local_prob_parm.d_Y_OH[i] = Y_OHinput[i];
       local_prob_parm.d_Y_H2O[i] = Y_H2Oinput[i];
       local_prob_parm.d_Y_HO2[i] = Y_HO2input[i];
       local_prob_parm.d_Y_H2O2[i] = Y_H2O2input[i];
       local_prob_parm.d_Y_N2[i] = Y_N2input[i];
    }

  // Initialize PeleLM::prob_parm container
  PeleLM::prob_parm->d_Tinput = (amrex::Real*)amrex::The_Arena()->alloc(nx * ny * nz * sizeof(amrex::Real));

  // Copy into PeleLM::prob_parm: CPU only for now
  std::memcpy(&PeleLM::prob_parm->d_Tinput, &local_prob_parm.d_Tinput,sizeof(local_prob_parm.d_Tinput));
  std::memcpy(&PeleLM::prob_parm->d_Y_H, &local_prob_parm.d_Y_H,sizeof(local_prob_parm.d_Y_H));
  std::memcpy(&PeleLM::prob_parm->d_Y_H2, &local_prob_parm.d_Y_H2,sizeof(local_prob_parm.d_Y_H2));
  std::memcpy(&PeleLM::prob_parm->d_Y_O, &local_prob_parm.d_Y_O,sizeof(local_prob_parm.d_Y_O));
  std::memcpy(&PeleLM::prob_parm->d_Y_O2, &local_prob_parm.d_Y_O2,sizeof(local_prob_parm.d_Y_O2));
  std::memcpy(&PeleLM::prob_parm->d_Y_OH, &local_prob_parm.d_Y_OH,sizeof(local_prob_parm.d_Y_OH));
  std::memcpy(&PeleLM::prob_parm->d_Y_H2O, &local_prob_parm.d_Y_H2O,sizeof(local_prob_parm.d_Y_H2O));
  std::memcpy(&PeleLM::prob_parm->d_Y_HO2, &local_prob_parm.d_Y_HO2,sizeof(local_prob_parm.d_Y_HO2));
  std::memcpy(&PeleLM::prob_parm->d_Y_H2O2, &local_prob_parm.d_Y_H2O2,sizeof(local_prob_parm.d_Y_H2O2));
  std::memcpy(&PeleLM::prob_parm->d_Y_N2, &local_prob_parm.d_Y_N2,sizeof(local_prob_parm.d_Y_N2));
  PeleLM::prob_parm->input_resolution = local_prob_parm.input_resolution;
  PeleLM::prob_parm->input_nx = local_prob_parm.input_nx;
  PeleLM::prob_parm->input_ny = local_prob_parm.input_ny;
  PeleLM::prob_parm->input_nz = local_prob_parm.input_nz;
  PeleLM::prob_parm->urms0 = local_prob_parm.urms0;
  PeleLM::prob_parm->uin_norm = local_prob_parm.uin_norm;
}
