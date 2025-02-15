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

  Omega_Nu0 = N_eff * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0) * Omega_gamma0;

  Omega_Lambda0 = 1 - (Omega_K0 + Omega_B0 + Omega_CDM0 + Omega_gamma0 + Omega_Nu0);

}

//====================================================
// Do all the solving. Compute eta(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    

  Vector x_array = Utils::linspace(x_start, x_end, num_x_points);
  //Vector a_array = exp(x_array);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    detadx[0] = Constants.c/Hp_of_x(x);

    return GSL_SUCCESS;
  };

  double eta_ini = (Constants.c*exp(x_start))/(H0*(Omega_gamma0 + Omega_Nu0));
  Vector y_ic{eta_ini};

  ODESolver ode;
  ode.solve(detadx, x_array, y_ic);
  auto eta_array = ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array, "eta of x Spline");

  Utils::EndTiming("Eta");
}

//====================================================
// Get methods
//====================================================

double BackgroundCosmology::H_of_x(double x) const{
  // Hubble parameter
  double a = exp(x);
  double matter_term = (Omega_B0 + Omega_CDM0)*pow(a,-3);
  double radiation_term = (Omega_gamma0 + Omega_Nu0)*pow(a,-4);
  double curvature_term = Omega_K0*pow(a,-2);
  double dark_energy_term = Omega_Lambda0;
  double H = H0*sqrt(matter_term + radiation_term + curvature_term + dark_energy_term);

  return H;
}

double BackgroundCosmology::Hp_of_x(double x) const{
  // H prime, modified Hubble parameter
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
  double a = exp(x);
  double Omega_B = ( Omega_B0*pow(H0,2) )/( pow(a,3)*pow(H_of_x(x),2) );

  return Omega_B;
}

double BackgroundCosmology::get_Omega_gamma(double x) const{
  if(x == 0.0) return Omega_gamma0;
  double a = exp(x);
  double Omega_gamma = ( Omega_gamma0*pow(H0,2) )/( pow(a,4)*pow(H_of_x(x),2) );

  return Omega_gamma;
}

double BackgroundCosmology::get_Omega_Nu(double x) const{
  if(x == 0.0) return Omega_Nu0;
  double a = exp(x);
  double Omega_Nu = ( Omega_Nu0*pow(H0,2) )/( pow(a,4)*pow(H_of_x(x),2) );

  return Omega_Nu;
}

double BackgroundCosmology::get_Omega_CDM(double x) const{
  if(x == 0.0) return Omega_CDM0;
  double a = exp(x);
  double Omega_CDM = ( Omega_CDM0*pow(H0,2) )/( pow(a,3)*pow(H_of_x(x),2) );

  return Omega_CDM;
}

double BackgroundCosmology::get_Omega_Lambda(double x) const{
  if(x == 0.0) return Omega_Lambda0;
  double a = exp(x);
  double Omega_Lambda = ( Omega_Lambda0*pow(H0,2) )/( pow(H_of_x(x),2) );

  return Omega_Lambda;
}

double BackgroundCosmology::get_Omega_K(double x) const{
  if(x == 0.0) return Omega_K0;
  double a = exp(x);
  double Omega_K = ( Omega_K0*pow(H0,2) )/( pow(a,2)*pow(H_of_x(x),2) );

  return Omega_K;
}

double BackgroundCosmology::get_sum_DensityParams(double x) const{
  if (x == 0.0) {
    return Omega_B0 + Omega_CDM0 + Omega_gamma0 + Omega_Nu0 + Omega_K0 + Omega_Lambda0;
  }
  // mass parameter
  double Omega_M = get_Omega_B(x) + get_Omega_CDM(x);
  // radiation parameter
  double Omega_R = get_Omega_gamma(x) + get_Omega_Nu(x);
  // sum
  return Omega_M + Omega_R + get_Omega_K(x) + get_Omega_Lambda(x);
}

double BackgroundCosmology::get_r_distance_of_x(double x) const {
  double Chi = get_comoving_distance_of_x(x);
  if (Omega_K0 == 0.0) return Chi;
  else {
    double A = sqrt(abs(Omega_K0))*H0*(Chi/Constants.c);

    if (Omega_K0 < 0.0) return Chi*(sin(A)/A);
    if (Omega_K0 > 0.0) return Chi*(sinh(A)/A);
  }
}


double BackgroundCosmology::get_angular_diameter_distance_of_x(double x) const{

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

  double Chi = eta_of_x(x_end) - eta_of_x(x);

  return Chi;
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
  std::cout << "Sum of dens. params.:  " << get_sum_DensityParams(0.0) << "\n";
  std::cout << "Omega_B0:      " << Omega_B0      << "\n";
  std::cout << "Omega_CDM0:    " << Omega_CDM0    << "\n";
  std::cout << "Omega_Lambda+: " << Omega_Lambda0 << "\n";
  std::cout << "Omega_K0:      " << Omega_K0      << "\n";
  std::cout << "Omega_Nu0:     " << Omega_Nu0     << "\n";
  std::cout << "Omega_gamma0:  " << Omega_gamma0  << "\n";
  std::cout << "N_eff:         " << N_eff         << "\n";
  std::cout << "h:             " << h             << "\n";
  std::cout << "T_CMB+:        " << T_CMB0        << "\n";
  std::cout << "H0:            " << H0            << "\n";
  std::cout << "H early radiation era:     " << Hp_of_x(Constants.x_start) << "\n";
  std::cout << "H little bit before current era:     " << Hp_of_x(0.00001) << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string& filename) const{
  const double x_min = -10.0;
  const double x_max =  0.0;
  const int    n_pts =  100;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  fp << "x | eta | Hp | dHp/dx | Sum_Dens_Params | Omega_B | Omega_CDM | Omega_Lambda | Omega_gamma | Omega_Nu | Omega_K |\n";
  auto print_data = [&] (const double x) {
    fp << x                        << " | ";
    fp << eta_of_x(x)              << " | ";
    fp << Hp_of_x(x)               << " | ";
    fp << dHpdx_of_x(x)            << " | ";
    fp << get_sum_DensityParams(x) << " | ";
    fp << get_Omega_B(x)           << " | ";
    fp << get_Omega_CDM(x)         << " | ";
    fp << get_Omega_Lambda(x)      << " | ";
    fp << get_Omega_gamma(x)       << " | ";
    fp << get_Omega_Nu(x)          << " | ";
    fp << get_Omega_K(x)           << " | ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

