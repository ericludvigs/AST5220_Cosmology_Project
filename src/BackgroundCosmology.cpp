#include "BackgroundCosmology.h"

//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double Omega_B0,
    double Omega_CDM0,
    double Omega_K0,
    double N_eff,
    double T_CMB0) :
  h(h),
  Omega_B0(Omega_B0),
  Omega_CDM0(Omega_CDM0),
  Omega_K0(Omega_K0),
  N_eff(N_eff),
  T_CMB0(T_CMB0)
{

  H0 = h*Constants.H0_over_h;

  double gamma_term1 = 2.0 * (pow(M_PI,2))/30.0;
  double gamma_term2 = (pow(Constants.k_b * T_CMB0,4))/
    (pow(Constants.hbar,3)* pow(Constants.c,5));
  double gamma_term3 = (8*M_PI*Constants.G)/(3*pow(H0,2));
  Omega_gamma0 = gamma_term1 * gamma_term2 * gamma_term3;


  //=============================================================================
  // TODO: Compute Omega_gamma0, Omega_Nu0, OmegaLambda, H0, ...
  //=============================================================================
  //...
  //...
  //...
  //...
}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    
  //=============================================================================
  // TODO: Set the range of x and the number of points for the splines
  // For this Utils::linspace(x_start, x_end, npts) is useful
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, num_x_points);
  Vector a_array = exp(x_array);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    //=============================================================================
    // TODO: Set the rhs of the detadx ODE
    //=============================================================================
    //...
    //...

    detadx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set the initial condition, set up the ODE system, solve and make
  // the spline eta_of_x_spline 
  //=============================================================================
  // ...
  // ...
  // ...
  // ...

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...
  double a = exp(x);
  double H = 0.0;

  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  double a = exp(x);
  double H = H_of_x(x);
  double Hp = a*H;

  return Hp;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::ddHpddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omega_B(double x) const{
  if(x == 0.0) return Omega_B0;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omega_gamma(double x) const{
  if(x == 0.0) return Omega_gamma0;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omega_Nu(double x) const{
  if(x == 0.0) return Omega_Nu0;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omega_CDM(double x) const{
  if(x == 0.0) return Omega_CDM0;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omega_Lambda(double x) const{
  if(x == 0.0) return Omega_Lambda0;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::get_Omega_K(double x) const{
  if(x == 0.0) return Omega_K0;

  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  //=============================================================================
  // TODO: Implement...
  //=============================================================================
  //...
  //...

  return 0.0;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_N_eff() const{
  return N_eff;
}

double BackgroundCosmology::get_T_CMB(double x) const{
  if(x == 0.0) return T_CMB0;
  return T_CMB0 * exp(-x);
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "Omega_B0:      " << Omega_B0      << "\n";
  std::cout << "Omega_CDM0:    " << Omega_CDM0    << "\n";
  std::cout << "Omega_Lambda+: " << Omega_Lambda0 << "\n";
  std::cout << "Omega_K0:      " << Omega_K0      << "\n";
  std::cout << "Omega_Nu0:     " << Omega_Nu0     << "\n";
  std::cout << "Omega_gamma0:  " << Omega_gamma0  << "\n";
  std::cout << "N_eff:         " << N_eff          << "\n";
  std::cout << "h:             " << h             << "\n";
  std::cout << "T_CMB+:        " << T_CMB0        << "\n";
  std::cout << "H0:            " << H0            << "\n";
  std::cout << "Omega_gamma0:  " << Omega_gamma0  << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                   << " ";
    fp << eta_of_x(x)         << " ";
    fp << Hp_of_x(x)          << " ";
    fp << dHpdx_of_x(x)       << " ";
    fp << get_Omega_B(x)      << " ";
    fp << get_Omega_CDM(x)    << " ";
    fp << get_Omega_Lambda(x) << " ";
    fp << get_Omega_gamma(x)  << " ";
    fp << get_Omega_Nu(x)     << " ";
    fp << get_Omega_K(x)      << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

