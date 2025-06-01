#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector log_Xe_arr  = Vector(npts_rec_arrays);
  Vector log_ne_arr  = Vector(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){
    //==============================================================
    // Get X_e from solving the Saha equation
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    //std::cout << "iteration i: "<< i << ", Xe=" << Xe_current << ", ne=" << ne_current << "\n";

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit) {
      saha_regime = false;
      saha_regime_index = i;
      saha_regime_exit_x = x_array[i];
      std::cout << "exited saha regime at i=" << saha_regime_index << "\n";
      std::cout << "current electron fraction Xe = " << Xe_current << "\n";
      std::cout << std::endl;
    }

    if(saha_regime){
      
      //=============================================================================
      // Store the result we got from the Saha equation
      //=============================================================================
      log_Xe_arr[i] = log(Xe_current);
      log_ne_arr[i] = log(ne_current);

    } else {

      //==============================================================
      // Compute X_e from current time til today by solving
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================

      // copy of our x array from point i to the end
      // this should contain all points after saha approx stops applying
      Vector peebles_x_array = Vector(x_array.begin()+i, x_array.end());
      //std::cout << "peebles_x_array: " << peebles_x_array.size() << "\n";

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // Set up IC, solve the ODE and fetch the result
      //=============================================================================
      double Xe_ini = Xe_current;
      Vector y_ic{Xe_ini};

      peebles_Xe_ode.solve(dXedx, peebles_x_array, y_ic);
      auto Xe_solution = peebles_Xe_ode.get_data_by_component(0);
      //std::cout << "Xe_solution: " << Xe_solution.size() << "\n";
      Vector log_Xe_solution = log(Xe_solution);
      //std::cout << "log_Xe_solution: " << log_Xe_solution.size() << "\n";

      // copy over to Xe array
      // copy from `begin` to `end` into `destination` starting at `position`
      std::copy(log_Xe_solution.begin(), log_Xe_solution.end(), log_Xe_arr.begin()+i);
      //std::cout << "log_Xe_arr: " << log_Xe_arr.size() << "\n";

      // make array of same length as xe solution (not as long as the full xe array)
      auto ne_arr_solution = Vector(Xe_solution.begin(), Xe_solution.end());
      // apply the electron abundance function to every point of xe
      // count up the rest of the points that remain in loop
      for(int j = 0; j < ne_arr_solution.size(); j++) {
        ne_arr_solution[j] = get_electron_abundance(peebles_x_array[j], Xe_solution[j]);
      }
      // take the log for splining
      auto log_ne_arr_solution = log(ne_arr_solution);
      // overwrite our vector with new info
      std::copy(log_ne_arr_solution.begin(), log_ne_arr_solution.end(), log_ne_arr.begin()+i);

      // end looping, we did everything with Peebles in one big batch
      break;
    
    }
  }

  //=============================================================================
  // Spline the result. Implement and make sure the Xe_of_x, ne_of_x
  // functions are working
  //=============================================================================
  //std::cout << "x_array: " << x_array.size() << "\n";
  //std::cout << "log_Xe_arr: " << log_Xe_arr.size() << "\n";
  //std::cout << "log_ne_arr: " << log_ne_arr.size() << "\n";
  log_Xe_of_x_spline.create(x_array, log_Xe_arr, "Xe of x Spline");
  log_ne_of_x_spline.create(x_array, log_ne_arr, "ne of x Spline");

  Utils::EndTiming("Xe");
}

//====================================================
// Function to specifically get n_e
//====================================================
double RecombinationHistory::get_electron_abundance(const double x, const double Xe) const {
  const double a           = exp(x);
  // Physical constants
  const double G           = Constants.G;
  const double m_H         = Constants.m_H;

  // Fetch cosmological parameters
  const double H0 = cosmo->get_H0();
  const double OmegaB0 = cosmo->get_Omega_B(0.0);

  // number density of baryons - assume no heavier elements,
  // so this is approx. baryon density over number hydrogen mass
  const double rho_c0 = (3.0)/(8.0*M_PI*G)*pow(H0, 2.0);
  const double rho_b = (rho_c0 * OmegaB0)/pow(a, 3);
  const double n_b = rho_b/m_H;
  // approximation made
  const double n_H = (1-Yp)*n_b;

  // from definition
  const double ne = Xe * n_H;
  return ne;
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;

  // Fetch cosmological parameters
  const double H0 = cosmo->get_H0();
  const double OmegaB0 = cosmo->get_Omega_B(0.0);
  const double T_CMB = cosmo->get_T_CMB(x);
  // assumption we made in text
  const double T_b = T_CMB/a;

  // Electron fraction and number density
  double Xe;
  double ne;

  // number density of baryons - assume no heavier elements,
  // so this is approx. baryon density over number hydrogen mass
  const double rho_c0 = (3.0)/(8.0*M_PI*G)*pow(H0, 2.0);
  const double rho_b = (rho_c0 * OmegaB0)/pow(a, 3);
  const double n_b = rho_b/m_H;
  //=============================================================================
  // Compute Xe and ne from the Saha equation
  //=============================================================================
  //std::cout << "T_b = " << T_b << "\n";
  double rhs_1 = (1/n_b);
  //std::cout << "rhs_1 = " << rhs_1 << "\n";
  double rhs_2 = pow((m_e*T_b*k_b)/(2*M_PI*pow(hbar,2)), 3.0/2.0);
  //std::cout << "rhs_2 = " << rhs_2 << "\n";
  double rhs_3 = exp(-epsilon_0/(k_b*T_b));
  //std::cout << "rhs_3 = " << rhs_3 << "\n";
  double S_RHS = rhs_1*rhs_2*rhs_3;
  //std::cout << "S_RHS = " << S_RHS << "\n";
  // if R ~ 1
  if (abs(S_RHS) < 4.0/(1e-4)) {
    //std::cout << "S_RHS = " << S_RHS << " < " << 4.0/(1e-4) << "\n";
    // abc formula solution, rewritten somewhat for numerics
    Xe = (S_RHS/2.0)*(-1 + sqrt(1 + 4.0/S_RHS));
    //std::cout << "Xe = " << Xe << "\n";
  }
  else {
    //std::cout << "S_RHS = " << S_RHS << " > " << 4.0/(1e-4) << "\n";
    // taylor expansion approximation
    Xe = (S_RHS/2.0)*(-1 + 1 + (1.0/2.0)*(4.0/S_RHS));
  }

  ne = get_electron_abundance(x, Xe);

  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;

  // Cosmological parameters
  const double H0 = cosmo->get_H0();
  const double H = cosmo->H_of_x(x);
  const double OmegaB0 = cosmo->get_Omega_B(0.0);
  const double T_CMB = cosmo->get_T_CMB(x);
  // assumption we made in text
  const double T_b = T_CMB/a;

  //=============================================================================
  // Write the expression for dXedx
  //=============================================================================
  double epsilon_0_kb_tb = (epsilon_0)/(k_b*T_b);
  double phi_2_tb = 0.448 * log(epsilon_0_kb_tb);
  double alpha_2_tb = (8*c*sigma_T)/(sqrt(3*M_PI)) * sqrt(epsilon_0_kb_tb) * phi_2_tb;

  double beta_parenthesis_factor = (m_e*k_b*T_b)/(2*M_PI*pow(hbar, 2));
  double beta_tb_no_exp = alpha_2_tb * pow(beta_parenthesis_factor, 3.0/2.0);
  double beta_tb = beta_tb_no_exp * exp(-epsilon_0_kb_tb);
  double beta_2_tb = beta_tb_no_exp * exp(-(1.0/4.0)*epsilon_0_kb_tb); // condensed exponential from the 2 betas

  double n_b = (3*pow(H0, 2)*OmegaB0)/(8*M_PI*G*m_H*pow(a, 3));
  double n_H = (1-Yp) * n_b;
  double n_1s = (1-X_e) * n_H;

  double lambda_alpha = H * (pow(3*epsilon_0, 3))/(pow(8*M_PI, 2)*pow(c, 3)*pow(hbar, 3)*n_1s);

  double Cr_tb = (lambda_2s1s+lambda_alpha)/(lambda_2s1s+lambda_alpha+beta_2_tb);

  double rhs_1 = Cr_tb/H;
  double rhs_2 = beta_tb*(1-X_e) - n_H*alpha_2_tb*pow(X_e, 2);
  
  dXedx[0] = rhs_1*rhs_2;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================
void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  // redshift given to happen around z = 10 ish
  // which is x = -2.34 ish
  const int npts = 1000;
  const double x_just_before_reion = -3;
  const double x_just_after_reion = -2;
  // small addon to prevent final array from having the same point twice
  const double x_epsilon = 1e-4;
  Vector x_before_reion_array = Utils::linspace(x_start, x_just_before_reion-x_epsilon, npts);
  Vector x_during_reion_array = Utils::linspace(x_just_before_reion, x_just_after_reion, npts);
  Vector x_after_reion_array = Utils::linspace(x_just_after_reion+x_epsilon, x_end, npts);
  // concatenate
  Vector x_array;
  x_array.reserve( x_before_reion_array.size() + x_during_reion_array.size() + x_after_reion_array.size() ); // preallocate memory
  x_array.insert( x_array.end(), x_before_reion_array.begin(), x_before_reion_array.end() );
  x_array.insert( x_array.end(), x_during_reion_array.begin(), x_during_reion_array.end() );
  x_array.insert( x_array.end(), x_after_reion_array.begin(), x_after_reion_array.end() );

  // make a reversed x array
  // std::cout << "x array" << "\n";
  // for (auto i: x_array)
  //   std::cout << i << ' ';
  // std::cout << "\n";

  Vector reverse_x_array = Vector(npts*3);
  std::reverse_copy(x_array.begin(), x_array.end(), reverse_x_array.begin());

  // std::cout << "reversed x array" << "\n";
  // for (auto i: reverse_x_array)
  //   std::cout << i << ' ';
  // std::cout << "\n";


  // all in one go
  //Vector x_array = Utils::linspace(x_start, x_end, npts);
  //Vector reverse_x_array = Utils::linspace(x_end, x_start, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // the expression for dtaudx
    //=============================================================================
    const double H = cosmo->H_of_x(x);
    // Set the derivative for photon optical depth
    dtaudx[0] = -(Constants.c * ne_of_x(x) * Constants.sigma_T)/(H);

    return GSL_SUCCESS;
  };

  //=============================================================================
  // Set up and solve the ODE
  //=============================================================================
  // by starting at x = 0 (today), we know that tau(0) = 0
  double tau_ini = 0.0;
  Vector y_ic{tau_ini};

  ODESolver tau_ode;
  tau_ode.solve(dtaudx, reverse_x_array, y_ic);
  auto tau_array = tau_ode.get_data_by_component(0);
  // reverse array to match normal x
  std::reverse(tau_array.begin(), tau_array.end());

  //=============================================================================
  // Set up and solve sound horizon ODE
  //=============================================================================
  // lambda function for sound wave speed
  auto c_s = [&](const double x) {
    const double omega_gamma0 = cosmo->get_Omega_gamma(0.0);
    const double omega_b0 = cosmo->get_Omega_B(0.0);
    const double a = exp(x);
    const double R = (4*omega_gamma0)/(3*omega_b0*a);
    return Constants.c * sqrt((R)/(3*(1+R)));
  };
  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dsdx = [&](double x, const double *s, double *dsdx){

    //=============================================================================
    // the expression for dsdx
    //=============================================================================
    const double Hp = cosmo->Hp_of_x(x);
    // Set the derivative for photon optical depth
    dsdx[0] = (c_s(x))/(Hp);

    return GSL_SUCCESS;
  };

  double s_ini = c_s(x_start)/(cosmo->Hp_of_x(x_start));
  Vector s_ic{s_ini};

  ODESolver s_ode;
  s_ode.solve(dsdx, x_array, s_ic);
  auto s_array = s_ode.get_data_by_component(0);

  //=============================================================================
  // Compute visibility functions and spline everything
  //=============================================================================
  tau_of_x_spline.create(x_array, tau_array, "tau of x spline");

  // use analytical derivation
  auto dtau_dx_array = Vector(3*npts);
  auto dtau_func = [&](const double x){ return -(Constants.c * ne_of_x(x) * Constants.sigma_T)/(cosmo->H_of_x(x)); };
  std::transform(x_array.begin(), x_array.end(), dtau_dx_array.begin(), dtau_func);
  dtau_of_dx_spline.create(x_array, dtau_dx_array, "dtau of dx spline");

  // g_tilde with a transform applied to all points with tau as input variable
  auto gtilde_array = Vector(3*npts);
  auto gtilde_func = [&](const double x){ return -dtau_func(x)*exp(-tau_of_x_spline(x)); };
  std::transform(x_array.begin(), x_array.end(), gtilde_array.begin(), gtilde_func);

  g_tilde_of_x_spline.create(x_array, gtilde_array, "g tilde of x spline");

  sound_horizon_of_x_spline.create(x_array, s_array, "s of x spline");

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // Either from the tau-spline tau_of_x_spline.deriv_(x) or
  // from a separate spline if you choose to do this
  //=============================================================================

  //return tau_of_x_spline.deriv_x(x);
  return dtau_of_dx_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{
  //return tau_of_x_spline.deriv_xx(x);
  return dtau_of_dx_spline.deriv_x(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return exp(log_Xe_of_x_spline(x));
}

double RecombinationHistory::ne_of_x(double x) const{
  return exp(log_ne_of_x_spline(x));
}

double RecombinationHistory::sound_horizon_of_x(double x) const{
  return sound_horizon_of_x_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << "Saha approx ended at iteration i=" << saha_regime_index << "\n";
  std::cout << "With x = " << saha_regime_exit_x << "\n";
  std::cout << "Tau today: " << tau_of_x(0) << "\n";
  std::cout << std::endl;
}

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string& filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  fp << "x"              << "; ";
  fp << "Xe"             << "; ";
  fp << "ne"             << "; ";
  fp << "tau"            << "; ";
  fp << "dtau_dx"        << "; ";
  fp << "ddtau_ddx"      << "; ";
  fp << "g_tilde"        << "; ";
  fp << "dg_tilde_dx"    << "; ";
  fp << "ddg_tilde_ddx" << "; ";
  fp << "sound_horizon" << "; ";
  fp << "t(x)";
  fp << "\n";
  auto print_data = [&] (const double x) {
    fp << x                    << "; ";
    fp << Xe_of_x(x)           << "; ";
    fp << ne_of_x(x)           << "; ";
    fp << tau_of_x(x)          << "; ";
    fp << dtaudx_of_x(x)       << "; ";
    fp << ddtauddx_of_x(x)     << "; ";
    fp << g_tilde_of_x(x)      << "; ";
    fp << dgdx_tilde_of_x(x)   << "; ";
    fp << ddgddx_tilde_of_x(x) << "; ";
    fp << sound_horizon_of_x(x) << "; ";
    fp << cosmo->t_of_x(x);
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

