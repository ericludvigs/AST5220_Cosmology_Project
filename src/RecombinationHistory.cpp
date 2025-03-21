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
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr  = Vector(npts_rec_arrays);
  Vector ne_arr  = Vector(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit) {
      saha_regime = false;
      std::cout << "exited saha regime at i=" << i;
    }

    if(saha_regime){
      
      //=============================================================================
      // Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;
      ne_arr[i] = ne_current;

    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);
      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...
      double Xe_ini = (Constants.c*exp(x_start)); // TODO:
      Vector y_ic{Xe_ini};

      peebles_Xe_ode.solve(dXedx, x_array, y_ic);
      auto Xe_solution = peebles_Xe_ode.get_data_by_component(0);

      Xe_arr[i] = 1.0;
      double n_H = 1.0;
      ne_arr[i] = Xe_arr[i]*n_H;

      break;
    
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...
  log_Xe_of_x_spline.create(x_array, Xe_arr, "Xe of x Spline");

  Utils::EndTiming("Xe");
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
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double H0 = cosmo->get_H0();
  const double OmegaB0 = cosmo->get_Omega_B(0.0);
  const double T_CMB = cosmo->get_T_CMB(x);
  // assumption we made in text
  const double T_b = T_CMB/a;

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;

  // number density of baryons - assume no heavier elements,
  // so this is approx. baryon density over number hydrogen mass
  const double rho_c0 = (3.0)/(8.0*M_PI*G)*pow(H0, 2.0);
  const double rho_b = rho_c0 * OmegaB0;
  const double n_b = rho_b/m_H;
  // approximation made
  const double n_H = n_b;
  //=============================================================================
  // Compute Xe and ne from the Saha equation
  //=============================================================================
  double S_RHS = 1/n_b*pow((m_e*T_b)/(2*M_PI), 3/2)*exp(-epsilon_0/(k_b*T_b));
  // if R ~ 1
  if (abs(S_RHS) > 4.0/1e-4) {
    // abc formula solution, rewritten somewhat for numerics
    Xe = (S_RHS/2.0)*(-1 + sqrt(1 + 4.0/S_RHS));
  }
  else {
    // taylor expansion approximation
    Xe = (S_RHS/2.0)*(-1 + 1 + (1/2.0)*(4.0/S_RHS));
  }

  // from definition
  ne = Xe*n_H;

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
  // const double OmegaB      = cosmo->get_OmegaB();
  // ...
  // ...

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  //...
  //...
  
  dXedx[0] = 1.0;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================
void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth
    dtaudx[0] = 1.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  double tau_ini = (Constants.c*exp(x_start)); // TODO: this
  Vector y_ic{tau_ini};

  ODESolver tau_ode;
  tau_ode.solve(dtaudx, x_array, y_ic);
  auto tau_array = tau_ode.get_data_by_component(0);

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  tau_of_x_spline.create(x_array, tau_array, "tau of x spline");

  Vector gtilde_array = Vector(npts);

  auto gtilde_func = [&](double x){ return -tau_of_x_spline.deriv_x(x)*exp(-tau_of_x_spline(x)); };
  std::transform(x_array.begin(), x_array.end(), gtilde_array.begin(), gtilde_func);

  g_tilde_of_x_spline.create(x_array, gtilde_array, "g tilde of x spline");

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
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================

  return tau_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================

  return tau_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================

  return g_tilde_of_x_spline.deriv_x(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================

  return g_tilde_of_x_spline.deriv_xx(x);
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ne_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
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
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << "; ";
    fp << Xe_of_x(x)           << "; ";
    fp << ne_of_x(x)           << "; ";
    fp << tau_of_x(x)          << "; ";
    fp << dtaudx_of_x(x)       << "; ";
    fp << ddtauddx_of_x(x)     << "; ";
    fp << g_tilde_of_x(x)      << "; ";
    fp << dgdx_tilde_of_x(x)   << "; ";
    fp << ddgddx_tilde_of_x(x) << " ";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

