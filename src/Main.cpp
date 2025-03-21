#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // Background parameters
  double h           = 0.67;
  double Omega_B0      = 0.05;
  double Omega_CDM0    = 0.267;
  double Omega_K0      = 0.0;
  double N_eff        = 3.046;
  double T_CMB0        = 2.7255;

  // Recombination parameters
  double Yp          = 0.245;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, Omega_B0, Omega_CDM0, Omega_K0, N_eff, T_CMB0);
  cosmo.solve();
  cosmo.info("results/cosmology_params_today.json");
  
  // Output background evolution quantities
  cosmo.output("results/cosmology.csv");

  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  mcmc_fit_to_supernova_data("data/supernovadata.txt", "results/results_supernovafitting.txt");

  //=========================================================================
  // Module II
  //=========================================================================
  
  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  std::cout << "Start solving for recombination";
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("results/recombination.csv");
  
  // Remove when module is completed
  return 0;

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");
  
  // Remove when module is completed
  return 0;
  
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0;

  Utils::EndTiming("Everything");
}
