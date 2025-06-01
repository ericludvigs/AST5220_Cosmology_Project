#include"Perturbations.h"

//====================================================
// Constructors
//====================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//====================================================
// Do all the solving
//====================================================

void Perturbations::solve(){

  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();

  // Compute source functions and spline the result
  compute_source_functions();
}

//====================================================
// The main work: integrate all the perturbations
// and spline the results
//====================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //===================================================================
  // Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //===================================================================
  //Vector k_array(n_k);
  // log of endpoints to avoid underflow to 0
  // then exp to get actual values
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));

  /*
  std::cout << "k array:" << "\n[";
  for (auto i: k_array)
    std::cout << i << ",";
  std::cout << "]\n";
  */

  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // allocate containers for ODE results
  Vector delta_cdm_array(n_x * n_k);
  Vector delta_b_array(n_x * n_k);
  Vector v_cdm_array(n_x * n_k);
  Vector v_b_array(n_x * n_k);
  Vector Phi_array(n_x * n_k);
  Vector Pi_array(n_x * n_k);
  Vector Psi_array(n_x * n_k);
  Vector2D Theta_ells_arrays(Constants.n_ell_theta, Vector((n_x * n_k)));
  // each entry for [l] is akin to Vector Theta_array(n_x * n_k)

  // Loop over all wavenumbers
  for(int ik = 0; ik < n_k; ik++){

    bool progress_bar_print_active = false;
    // Progress bar...
    if ( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << "- " << (100*ik+100)/n_k << "% " << "-\n" << std::flush;
      progress_bar_print_active = true;
      if(ik == n_k-1) std::cout << std::endl;
    }
    else { progress_bar_print_active = false; }

    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    std::pair<double, int> result = get_tight_coupling_time(k);
    double x_end_tight = result.first;
    int x_end_index = result.second;

    if (progress_bar_print_active) {
      std::cout << "tight decoupling for k=" << k << " ends at x = "<< x_end_tight << "\n";
    }

    // allocate sub-arrays of x
    Vector x_tc_array(x_end_index+1); // from 0 to x_end_index
    Vector x_after_tc_array(n_x-x_end_index); // from x_end_index to n_x

    // copies exactly `count` values from the range beginning at `first` to the range beginning at `result`
    // for each integer `i` in `[0, count)`, performs `*(result + i) = *(first + i)`
    std::copy_n(x_array.begin(), x_end_index, x_tc_array.begin());

    // second array
    std::copy_n(x_array.begin()+x_end_index, n_x-x_end_index, x_after_tc_array.begin());

    //===================================================================
    // Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // Integrate from x_start -> x_end_tight
    //Vector x_tc_array = Utils::linspace(x_start, x_end_tight, n_x);

    ODESolver tight_coupling_ode;
    tight_coupling_ode.solve(dydx_tight_coupling, x_tc_array, y_tight_coupling_ini);
    auto y_tight_coupling = tight_coupling_ode.get_final_data();
    //auto y_tight_coupling_all_data = tight_coupling_ode.get_data();

    Vector delta_cdm_tc =  tight_coupling_ode.get_data_by_component(Constants.ind_deltacdm_tc);
    Vector delta_b_tc   =  tight_coupling_ode.get_data_by_component(Constants.ind_deltab_tc);
    Vector v_cdm_tc     =  tight_coupling_ode.get_data_by_component(Constants.ind_vcdm_tc);
    Vector v_b_tc       =  tight_coupling_ode.get_data_by_component(Constants.ind_vb_tc);
    Vector Phi_tc       =  tight_coupling_ode.get_data_by_component(Constants.ind_Phi_tc);
    std::vector<Vector> Theta_tc(Constants.n_ell_theta_tc);
    for (int l = 0; l < Constants.n_ell_theta_tc; l++) {
      Theta_tc[l] = tight_coupling_ode.get_data_by_component(Constants.ind_start_theta_tc + l);
    }

    int num_theta_elements = 0;
    std::cout << "Theta_0 during tight coupling:\n";
    for (double theta : Theta_tc[0]) {
      std::cout << theta << "\n";
      num_theta_elements++;
    }

    std::cout << "Num of elements in theta_0: " << num_theta_elements << "\n";
    std::cout << "Num of x values: " << x_tc_array.size() << "\n";

    //===================================================================
    // Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===================================================================

    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k);

    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
    };

    // Integrate from x_end_tight -> x_end
    //Vector x_after_tc_array = Utils::linspace(x_end_tight, x_end, n_x);

    ODESolver full_ode;
    full_ode.solve(dydx_full, x_after_tc_array, y_full_ini);
    auto y_after_tc_end = full_ode.get_final_data();
    //auto y_after_tc_all_data = full_ode.get_data();

    Vector delta_cdm_after_tc =  full_ode.get_data_by_component(Constants.ind_deltacdm);
    Vector delta_b_after_tc   =  full_ode.get_data_by_component(Constants.ind_deltab);
    Vector v_cdm_after_tc     =  full_ode.get_data_by_component(Constants.ind_vcdm);
    Vector v_b_after_tc       =  full_ode.get_data_by_component(Constants.ind_vb);
    Vector Phi_after_tc       =  full_ode.get_data_by_component(Constants.ind_Phi);
    std::vector<Vector> Theta_after_tc(Constants.n_ell_theta);
    for (int l = 0; l < Constants.n_ell_theta; l++) {
      Theta_after_tc[l] = full_ode.get_data_by_component(Constants.ind_start_theta + l);
    }

    std::cout << "Theta_1 after tight coupling:\n";
    for (double theta : Theta_after_tc[1]) {
      std::cout << theta << "\n";
    }

    //===================================================================
    // store the data found from integrating so we can
    // spline it below
    //
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================

    // loop over the x values, better to do a collected loop than individual loops with std::copy
    double H0 = cosmo->get_H0();
    const double Omega_gamma0 = cosmo->get_Omega_gamma(0.0);
    Vector Theta_tc_2(n_x);
    for(int ix = 0; ix < n_x; ix++) {
      // Cosmological parameters and variables
      const double a = exp(x_array[ix]);
      double curr_dtau_dx = rec->dtaudx_of_x(x_array[ix]);
      double curr_H_prime = cosmo->Hp_of_x(x_array[ix]);

      if (ix < x_end_index) {
        for (int l = 0; l < Constants.n_ell_theta_tc; l++) {
          Theta_ells_arrays[l][ix + n_x * ik] = Theta_tc[l][ix];
          // TODO: include the extra theta values in early times too from the fixed initial condition expressions
        }
        // having to manually set a value here is horrible
        // cant be done with ODE solution because it wants the derivative
        Theta_tc_2[ix] = -((20*Constants.c*k)/(45*curr_H_prime*curr_dtau_dx))*Theta_tc[1][ix];
        Theta_ells_arrays[2][ix + n_x * ik] = Theta_tc_2[ix];
        delta_cdm_array[ix + n_x * ik] = delta_cdm_tc[ix];
        delta_b_array[ix + n_x * ik] = delta_b_tc[ix];
        v_cdm_array[ix + n_x * ik] = v_cdm_tc[ix];
        v_b_array[ix + n_x * ik] = v_b_tc[ix];
        Phi_array[ix + n_x * ik] = Phi_tc[ix];
        Pi_array[ix + n_x * ik] = Theta_tc_2[ix]; // ignore polarization terms
        double psi = -Phi_tc[ix] - ((12.0*pow(H0,2))/(pow(Constants.c, 2)*pow(k,2)*pow(a,2))) * (Omega_gamma0*Theta_tc_2[ix]);
        Psi_array[ix + n_x * ik] = psi;
      }
      else { // ix >= x_end_index
        for (int l = 0; l < Constants.n_ell_theta; l++) {
          Theta_ells_arrays[l][ix + n_x * ik] = Theta_after_tc[l][ix - x_end_index];
        }
        delta_cdm_array[ix + n_x * ik] = delta_cdm_after_tc[ix - x_end_index];
        delta_b_array[ix + n_x * ik] = delta_b_after_tc[ix - x_end_index];
        v_cdm_array[ix + n_x * ik] = v_cdm_after_tc[ix - x_end_index];
        v_b_array[ix + n_x * ik] = v_b_after_tc[ix - x_end_index];
        Phi_array[ix + n_x * ik] = Phi_after_tc[ix - x_end_index];
        Pi_array[ix + n_x * ik] = Theta_after_tc[2][ix - x_end_index]; // ignore polarization terms
        double psi = -Phi_after_tc[ix - x_end_index] - ((12.0*pow(H0,2))/(pow(Constants.c, 2)*pow(k,2)*pow(a,2))) * (Omega_gamma0*Theta_after_tc[2][ix - x_end_index]);
        Psi_array[ix + n_x * ik] = psi;
      }
    }

  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
  delta_cdm_spline.create(x_array, k_array, delta_cdm_array);
  delta_b_spline.create(x_array, k_array, delta_b_array);
  v_cdm_spline.create(x_array, k_array, v_cdm_array);
  v_b_spline.create(x_array, k_array, v_b_array);
  Phi_spline.create(x_array, k_array, Phi_array);
  Pi_spline.create(x_array, k_array, Pi_array);
  Psi_spline.create(x_array, k_array, Psi_array);
  Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
  for (int l = 0; l < Constants.n_ell_theta; l++) {
    Theta_spline[l].create(x_array, k_array, Theta_ells_arrays[l]);
  }

}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  //const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  //const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  //const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  //double *Nu           = &y_tc[Constants.ind_start_nu_tc]; // neutrinos

  //=============================================================================
  // Set the initial conditions in the tight coupling regime
  //=============================================================================
  // no neutrinos
  double f_v = 0.0;
  // otherwise would be
  // double f_v = Omega_nu0/(Omega_gamma0 + Omega_nu0);
  double Hp = cosmo->Hp_of_x(x);
  double dtau_dx = rec->dtaudx_of_x(x);

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  const double psi = -1.0/(3.0/2.0); // f_v is zero
  Phi = -psi; // f_v is zero
  // no need to store psi since it gets derived from phi in use
  delta_b = -3.0/2.0 * psi;
  delta_cdm = delta_b;

  v_b = -psi*(Constants.c*k)/(2.0*Hp);
  v_cdm = v_b;

  // SET: Photon temperature perturbations (Theta_ell)
  for (int l = 0; l < n_ell_theta_tc; l++){
    if (l == 0) {
      Theta[l] = -(1.0/2.0)*psi;
    }
    if (l == 1) {
      Theta[l] = +((Constants.c*k)/(6.0*Hp))*psi;
    }
    if (l == 2) {
      // no polarization
      Theta[l] = -((20.0*Constants.c*k)/(45.0*Hp*dtau_dx))*Theta[l-1];
    }
    else {
      Theta[l] = -(l/(2.0*l+1))*((Constants.c*k)/(Hp*dtau_dx))*Theta[l-1];
    }
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // not needed
  }

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x,
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  //const int n_ell_thetap        = Constants.n_ell_thetap;
  //const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  //const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  //const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  //double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  //double *Nu              = &y[Constants.ind_start_nu_tc];

  //=============================================================================
  // fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================
  double Hp = cosmo->Hp_of_x(x);
  double dtau_dx = rec->dtaudx_of_x(x);

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  // tc ells < full theta ells, so set the ones that exist and use generic condition for rest
  for (int l = 0; l < n_ell_theta; l++){
    if (l <= n_ell_theta_tc) {
      Theta[l] = Theta_tc[l];
    }
    else {
      Theta[l] = -(l/(2.0*l+1))*((Constants.c*k)/(Hp*dtau_dx))*Theta[l-1];
    }
  }

  // SET: Photon polarization perturbations (Theta_p_ell)
  if(polarization){
    // nope
  }

  // SET: Neutrino perturbations (N_ell)
  if(neutrinos){
    // nope
  }

  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

std::pair<double, int> Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0; // the returned value for when tight coupling ends
  int x_index = 0;
  bool search_condition_end = false;

  double x_search_start = Constants.x_start;
  double x_search_end = -8.3; // redshift z=4000, don't bother searching beyond this
  Vector x_search_array = Utils::linspace(x_search_start, x_search_end, 100);

  //=============================================================================
  // compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin
  //=============================================================================
  for (auto x : x_search_array){
    double curr_dtau_dx = rec->dtaudx_of_x(x);
    double curr_H_prime = cosmo->Hp_of_x(x);

    if ( curr_dtau_dx > 10 ) {
      search_condition_end = true;
    }
    if ( curr_dtau_dx > 10*((Constants.c*k)/curr_H_prime) ) {
      search_condition_end = true;
    }

    // update with current x, if loop ends then this will be the end-point -8.3 and return
    // if condition is true before that, then breaks early
    x_tight_coupling_end = x;
    if (search_condition_end) {
      break;
    }
    x_index++;
  }

  return std::make_pair(x_tight_coupling_end, x_index);
}

//====================================================
// After integrating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){ // TODO
  Utils::StartTiming("source");

  //=============================================================================
  // Make the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  Vector k_array = exp(Utils::linspace(log(k_min), log(k_max), n_k));
  Vector x_array = Utils::linspace(x_start, x_end, n_x);

  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  //Vector SE_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      // const double Hp       = cosmo->Hp_of_x(x);
      // const double tau      = rec->tau_of_x(x);
      // ...
      // ...

      // Temperatur source
      ST_array[index] = 0.0;

      // Polarization source
      if(Constants.polarization){
        //SE_array[index] = 0.0;
      }
    }
  }

  // Spline the source functions
  ST_spline.create (x_array, k_array, ST_array, "Source_Temp_x_k");

  if(Constants.polarization){
    //SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  //const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  //const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  //double *dNudx           = &dydx[Constants.ind_start_nu_tc];

  //=============================================================================
  // expressions for all the derivatives
  //=============================================================================

  const double a = exp(x);
  double H0 = cosmo->get_H0();
  double Hp = cosmo->Hp_of_x(x);
  double dHp_dx = cosmo->dHpdx_of_x(x);
  double dtau_dx = rec->dtaudx_of_x(x);
  double ddtau_dx = rec->ddtauddx_of_x(x);
  const double Omega_gamma0 = cosmo->get_Omega_gamma(0.0);
  const double Omega_b0 = cosmo->get_Omega_B(0.0);
  const double Omega_CDM0 = cosmo->get_Omega_CDM(0.0);

  const double R = (4.0*Omega_gamma0)/(3.0*Omega_b0*a);
  double ck_HP = (Constants.c*k)/(Hp);

  // SET: Scalar quantities (Phi, delta, v, ...)
  double psi = -Phi - ((12.0*pow(H0,2))/(pow(Constants.c, 2)*pow(k,2)*pow(a,2))) * (Omega_gamma0*Theta[2]);
  double ck_Hp_psi = ck_HP*psi;

  double dphi_factor = (Omega_CDM0*pow(a,-1)*delta_cdm + Omega_b0*pow(a,-1)*delta_b + 4.0*Omega_gamma0*pow(a,-2)*Theta[0]);
  dPhidx = psi - (1.0/3.0)*pow(ck_HP, 2)*Phi + ((pow(H0, 2))/(2.0*pow(Hp,2)))*dphi_factor;

  ddelta_cdmdx = ck_HP*v_cdm - 3.0*dPhidx;
  ddelta_bdx = ck_HP*v_b - 3.0*dPhidx;
  dv_cdmdx = -v_cdm - ck_Hp_psi;
  //dv_bdx = ; // modified for tight coupling

  double q_top_1 = -((1.0 - R)*dtau_dx + (1.0 + R)*ddtau_dx) * (3*Theta[1]+v_b);
  double q_top_2 = (1.0 - (dHp_dx)/(Hp)) * ck_HP * (-Theta[0] + 2.0*Theta[2]);
  double q_top_3 = ck_HP*dThetadx[0];
  double q_top = q_top_1 - ck_Hp_psi + q_top_2 - q_top_3;
  double q_bottom = (1.0 + R)*dtau_dx + (dHp_dx)/(Hp) - 1.0;
  double q = q_top/q_bottom;

  dv_bdx = 1.0/(1.0+R)*(-v_b - ck_Hp_psi + R*(q+ck_HP*(-Theta[0]+2.0*Theta[2]) - ck_Hp_psi));

  // SET: Photon multipoles (Theta_ell)
  for (int l = 0; l < n_ell_theta_tc; l++){
    if (l == 0) {
      dThetadx[l] = -(Constants.c*k)/(Hp)*Theta[1] - dPhidx;
    }
    if (l == 1) {
      dThetadx[l] = (1.0/3.0)*(q-dv_bdx);
    }
    else {
      std::cout << "ERROR: max l too high for tight coupling regime" << std::endl;
    }
  }

  // SET: Neutrino multipoles (Nu_ell)
  if(neutrinos){
    // nope
  }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  //const int n_ell_thetap        = Constants.n_ell_thetap;
  //const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  //const double *Theta_p         = &y[Constants.ind_start_thetap];
  //const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  //double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  //double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  const double a = exp(x);
  double H0 = cosmo->get_H0();
  double Hp = cosmo->Hp_of_x(x);
  double dHp_dx = cosmo->dHpdx_of_x(x);
  const double Omega_gamma0 = cosmo->get_Omega_gamma(0.0);
  const double Omega_b0 = cosmo->get_Omega_B(0.0);
  const double Omega_CDM0 = cosmo->get_Omega_CDM(0.0);
  const double eta = cosmo->eta_of_x(x);

  // Recombination variables
  double dtau_dx = rec->dtaudx_of_x(x);
  double ddtau_dx = rec->ddtauddx_of_x(x);

  //=============================================================================
  // fill in the expressions for all the derivatives
  //=============================================================================
  const double R = (4.0*Omega_gamma0)/(3.0*Omega_b0*a);
  double ck_HP = (Constants.c*k)/(Hp);

  // SET: Scalar quantities (Phi, delta, v, ...)
  double psi = -Phi - ((12.0*pow(H0,2))/(pow(Constants.c, 2)*pow(k,2)*pow(a,2))) * (Omega_gamma0*Theta[2]);
  double ck_Hp_psi = ck_HP*psi;

  double dphi_factor = (Omega_CDM0*pow(a,-1)*delta_cdm + Omega_b0*pow(a,-1)*delta_b + 4.0*Omega_gamma0*pow(a,-2)*Theta[0]);
  dPhidx = psi - (1.0/3.0)*pow(ck_HP, 2)*Phi + ((pow(H0, 2))/(2.0*pow(Hp,2)))*dphi_factor;

  ddelta_cdmdx = ck_HP*v_cdm - 3.0*dPhidx;
  ddelta_bdx = ck_HP*v_b - 3.0*dPhidx;
  dv_cdmdx = -v_cdm - ck_Hp_psi;
  dv_bdx = -v_b - ck_Hp_psi + dtau_dx*R*(3*Theta[1] + v_b);

  // SET: Photon multipoles (Theta_ell)
  for (int l = 0; l < n_ell_theta; l++){
    if (l == 0) {
      dThetadx[l] = -ck_HP*Theta[1] - dPhidx;
    }
    if (l == 1) {
      dThetadx[l] = (1.0/3.0)*(ck_HP*Theta[0]) - (2.0/3.0)*(ck_HP*Theta[2]) + (1.0/3.0)*ck_Hp_psi + dtau_dx*(Theta[1]+(1.0/3.0)*v_b);
    }
    if (l >= 2 && l < n_ell_theta-1) {
      double dthetadx_1 = (l/(2.0*l+1.0))*ck_HP*Theta[l-1];
      double dthetadx_2 = ((l+1.0)/(2.0*l+1.0))*ck_HP*Theta[l+1];
      double Pi_delta_l2 = (l==2) ? Theta[2]: 0.0; // delta_l2 -> 1 if 2 and 0 otherwise
      double dthetadx_3 = dtau_dx*(Theta[l] - (1.0/10.0)*Pi_delta_l2);
      dThetadx[l] = dthetadx_1 - dthetadx_2 + dthetadx_3;
    }
    else {
      dThetadx[l] = ck_HP*Theta[l-1] - Constants.c*(l+1.0)/(Hp*eta)*Theta[l] + dtau_dx*Theta[l];
    }
  }

  // SET: Photon polarization multipoles (Theta_p_ell)
  if(polarization){
    // nope
  }

  // SET: Neutrino multipoles (Nu_ell)
  if(neutrinos){
    // nope
  }

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}
double Perturbations::get_Pi(const double x, const double k) const{
  return Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
// double Perturbations::get_Source_E(const double x, const double k) const{
//   return SE_spline(x,k);
// }
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
// double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
//   return Theta_p_spline[ell](x,k);
// }
// double Perturbations::get_Nu(const double x, const double k, const int ell) const{
//   return Nu_spline[ell](x,k);
// }

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:     " << n_x              << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:     " << n_k              << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string& filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);
  fp << "x"              << "; ";
  fp << "Theta_0"             << "; ";
  fp << "Theta_1"             << "; ";
  fp << "Theta_2"            << "; ";
  fp << "delta_cdm"          << "; ";
  fp << "delta_b"            << "; ";
  fp << "v_cdm"              << "; ";
  fp << "v_b"                << "; ";
  fp << "Phi"        << "; ";
  fp << "Psi"      << "; ";
  fp << "Pi"       << "; ";
  fp << "Source_T_0"  << "; ";
  fp << "Source_T_5"  << "; ";
  fp << "Source_T_50" << "; ";
  fp << "Source_T_500" << "; ";
  fp << "\n";
  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << "; ";
    fp << get_Theta(x,k,0)   << "; ";
    fp << get_Theta(x,k,1)   << "; ";
    fp << get_Theta(x,k,2)   << "; ";
    fp << get_delta_cdm(x,k) << "; ";
    fp << get_delta_b(x,k)   << "; ";
    fp << get_v_cdm(x,k)     << "; ";
    fp << get_v_b(x,k)       << "; ";
    fp << get_Phi(x,k)       << "; ";
    fp << get_Psi(x,k)       << "; ";
    fp << get_Pi(x,k)        << "; ";
    fp << get_Source_T(x,k)  << "; ";
    fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << "; ";
    fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << "; ";
    fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << "";
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

