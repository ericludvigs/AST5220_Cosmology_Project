#ifndef _BACKGROUNDCOSMOLOGY_HEADER
#define _BACKGROUNDCOSMOLOGY_HEADER
#include <iostream>
#include <fstream>
#include "Utils.h"

using Vector = std::vector<double>;

class BackgroundCosmology{
  private:
   
    // Cosmological parameters
    double h;                       // Little h = H0/(100km/s/Mpc)
    double Omega_B0;                // Baryon density today
    double Omega_CDM0;              // CDM density today
    double Omega_Lambda0;           // Dark energy density today
    double N_eff;                    // Effective number of relativistic species (3.046 or 0 if ignoring neutrinos)
    double T_CMB0;                  // Temperature of the CMB today in Kelvin
   
    // Derived parameters
    double Omega_gamma0;            // Photon density today (follows from T_CMB)
    double Omega_Nu0;               // Neutrino density today (follows from T_CMB and N_eff)
    double Omega_K0;                // Curvature density = 1 - OmegaM - OmegaR - OmegaNu - OmegaLambda
    double H0;                      // The Hubble parameter today H0 = 100h km/s/Mpc

    // Helpers to make code look better
    double gamma_term1;
    double gamma_term2;
    double gamma_term3;

    // Start and end of x-integration (can be changed)
    double x_start = Constants.x_start;
    double x_end   = Constants.x_end;
    double num_x_points = Constants.num_x_points;

    // Splines to be made
    Spline eta_of_x_spline{"eta"};
    Spline t_of_x_spline{"t(x)"};
 
  public:

    // Constructors 
    BackgroundCosmology() = delete;
    BackgroundCosmology(
        double h, 
        double Omega_B0,
        double Omega_CDM0,
        double Omega_K0,
        double N_eff,
        double T_CMB0
        );

    // Print some useful info about the class
    void info() const;

    // Do all the solving
    void solve();

    // Output some results to file
    void output(const std::string& filename) const;

    // Get functions that we must implement
    double eta_of_x(double x) const;
    double t_of_x(double x) const;
    double H_of_x(double x) const;
    double Hp_of_x(double x) const;
    double dHpdx_of_x(double x) const;
    double ddHpddx_of_x(double x) const;
    double get_Omega_B(double x = 0.0) const;
    double get_Omega_M(double x = 0.0) const;
    double get_Omega_gamma(double x = 0.0) const;
    double get_Omega_R_tot(double x = 0.0) const;
    double get_Omega_Nu(double x = 0.0) const;
    double get_Omega_CDM(double x = 0.0) const;
    double get_Omega_Lambda(double x = 0.0) const;
    double get_Omega_K(double x = 0.0) const;

    double get_sum_DensityParams(double x = 0.0) const;

    double get_Omega_Mnu(double x = 0.0) const;
    double get_H0() const;
    double get_h() const;
    double get_N_eff() const;
    double get_T_CMB(double x = 0.0) const;

    // Distance measures
    double get_r_distance_of_x(double x) const;
    double get_angular_diameter_distance_of_x(double x) const;
    double get_luminosity_distance_of_x(double x) const;
    double get_comoving_distance_of_x(double x) const;

};

#endif
