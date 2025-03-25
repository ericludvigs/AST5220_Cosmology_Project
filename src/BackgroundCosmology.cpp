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
// Do all the solving. Compute eta(x), t(x)
//====================================================

// Solve the background
void BackgroundCosmology::solve(){
  //Utils::StartTiming("Eta");

  Vector x_array = Utils::linspace(x_start, x_end, num_x_points);
  //Vector a_array = exp(x_array);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    detadx[0] = Constants.c/Hp_of_x(x);

    return GSL_SUCCESS;
  };

  double eta_ini = (Constants.c*exp(x_start))/(H0*(Omega_gamma0 + Omega_Nu0));
  Vector y_ic{eta_ini};

  ODESolver eta_ode;
  eta_ode.solve(detadx, x_array, y_ic);
  auto eta_array = eta_ode.get_data_by_component(0);

  eta_of_x_spline.create(x_array, eta_array, "eta of x Spline");

  //Utils::EndTiming("Eta");

  //Utils::StartTiming("t(x)");
  // The ODE for dt/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    dtdx[0] = 1.0/H_of_x(x);

    return GSL_SUCCESS;
  };

  double t_ini = 1.0/(2*(H_of_x(x_start)));
  Vector t_ic{t_ini};

  ODESolver t_ode;
  t_ode.solve(dtdx, x_array, t_ic);
  auto t_array = t_ode.get_data_by_component(0);

  t_of_x_spline.create(x_array, t_array, "t of x Spline");

  //Utils::EndTiming("t(x)");
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
  const double a = exp(x);
  const double Omega_B = ( Omega_B0*pow(H0,2) )/( pow(a,3)*pow(H_of_x(x),2) );

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
  // compiler throws a fit about not having a return here
  throw std::invalid_argument("Omega_K0 is somehow neither 0, nor smaller or greater than 0");
}


double BackgroundCosmology::get_angular_diameter_distance_of_x(double x) const{
  // d_A = a*r
  return exp(x)*get_r_distance_of_x(x);
}
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  // d_L = r/a
  return get_r_distance_of_x(x)/exp(x);
}
double BackgroundCosmology::get_comoving_distance_of_x(double x) const{
  // Chi = eta_0 - eta
  double Chi = eta_of_x(x_end) - eta_of_x(x);

  return Chi;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
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
void BackgroundCosmology::info(const std::string& filename = "") const{
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "Sum of dens. params.:  " << std::setw(5) << get_sum_DensityParams(0.0) << "\n";
  std::cout << "Omega_B0:      " << std::setw(15) << Omega_B0      << "\n";
  std::cout << "Omega_CDM0:    " << std::setw(15) << Omega_CDM0    << "\n";
  std::cout << "Omega_Lambda0: " << std::setw(15) << Omega_Lambda0 << "\n";
  std::cout << "Omega_K0:      " << std::setw(15) << Omega_K0      << "\n";
  std::cout << "Omega_Nu0:     " << std::setw(15) << Omega_Nu0     << "\n";
  std::cout << "Omega_gamma0:  " << std::setw(15) << Omega_gamma0  << "\n";
  std::cout << "N_eff:         " << std::setw(15) << N_eff         << "\n";
  std::cout << "h:             " << std::setw(15) << h             << "\n";
  std::cout << "T_CMB0:        " << std::setw(15) << T_CMB0        << "\n";
  std::cout << "H0:            " << std::setw(15) << H0            << "\n";
  std::cout << "H early radiation era:           " << std::setw(15) << Hp_of_x(Constants.x_start) << "\n";
  std::cout << "H little bit before current era: " << std::setw(15) << Hp_of_x(0.00001) << "\n";
  std::cout << std::endl;

  if (!filename.empty()) {
    std::ofstream fp(filename.c_str());
    fp << "{" << "\n";
    fp << "    \"Omega_B0\": " << Omega_B0 << ",\n";
    fp << "    \"Omega_CDM0\": " << Omega_CDM0 << ",\n";
    fp << "    \"Omega_Lambda0\": " << Omega_Lambda0 << ",\n";
    fp << "    \"Omega_K0\": " << Omega_K0 << ",\n";
    fp << "    \"Omega_Nu0\": " << Omega_Nu0 << ",\n";
    fp << "    \"Omega_gamma0\": " << Omega_gamma0 << ",\n";
    fp << "    \"N_eff\": " << N_eff << ",\n";
    fp << "    \"T_CMB0\": " << T_CMB0 << ",\n";
    fp << "    \"h\": " << h << ",\n";
    fp << "    \"H0\": " << H0 << "\n";
    fp << "}" << "\n";
  }
}

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string& filename) const{
  const double x_min = -12.0;
  const double x_max =  0.0;
  const int    n_pts =  1000;
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  fp << "      x      |      eta     |      t       |      Hp      |    dHp/dx    |   ddHp/ddx   |      Chi     |      d_A     |      d_L     | Sum_Dens_Params |    Omega_B   |   Omega_CDM  | Omega_Lambda |  Omega_gamma |   Omega_Nu   |   Omega_K    \n";
  auto print_data = [&] (const double x) {
    fp << std::setw(12) << x                                     << " | ";
    fp << std::setw(12) << eta_of_x(x)                           << " | ";
    fp << std::setw(12) << t_of_x(x)                             << " | ";
    fp << std::setw(12) << Hp_of_x(x)                            << " | ";
    fp << std::setw(12) << dHpdx_of_x(x)                         << " | ";
    fp << std::setw(12) << ddHpddx_of_x(x)                       << " | ";
    fp << std::setw(12) << get_comoving_distance_of_x(x)         << " | ";
    fp << std::setw(12) << get_angular_diameter_distance_of_x(x) << " | ";
    fp << std::setw(12) << get_luminosity_distance_of_x(x)       << " | ";
    fp << std::setw(15) << get_sum_DensityParams(x)              << " | ";
    fp << std::setw(12) << get_Omega_B(x)                        << " | ";
    fp << std::setw(12) << get_Omega_CDM(x)                      << " | ";
    fp << std::setw(12) << get_Omega_Lambda(x)                   << " | ";
    fp << std::setw(12) << get_Omega_gamma(x)                    << " | ";
    fp << std::setw(12) << get_Omega_Nu(x)                       << " | ";
    fp << std::setw(12) << get_Omega_K(x)                        << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

