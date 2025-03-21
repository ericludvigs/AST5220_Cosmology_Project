import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append("..")
import python_methods.constants as constants

show_plots = False

cosmology_data_filename = Path("../results/cosmology.csv")
supernova_data_filename = Path("../data/supernovadata.txt").resolve()
supernova_fitting_data_filename = Path("../results/results_supernovafitting.txt").resolve()
plots_folder = Path("Plots/")

cosmology_df = pd.read_csv(cosmology_data_filename, sep="|", header=0, index_col=0)
cosmology_df.columns = cosmology_df.columns.str.strip()
print(cosmology_df)

# make an array with redshifts, reverse it and cut it down to match with the input data
z_array = constants.redshift(cosmology_df.index.to_numpy())[::-1]
# print(f"{z_array=}")

supernovadata = np.loadtxt(supernova_data_filename)
supernovadata_z = supernovadata[:, 0]
supernovadata_dL = supernovadata[:, 1]
supernovadata_error = supernovadata[:, 2]
print(f"{supernovadata=}")
print(f"{supernovadata_z=}")

z_array = z_array[z_array < np.max(supernovadata_z) + 0.02]
print(f"{z_array=}")

# get the matching lum distances we calculated, also in Gpc
# need to reverse, and only select the cosmologically "closest" elements to match with measured distances
calculated_dL = (cosmology_df["d_L"].to_numpy()[::-1] / constants.Gpc)[: len(z_array)]

print(f"{calculated_dL=}")

N = supernovadata.shape[0]
print(f"Number of datapoints: {N}")

# skip rows with burn in time, skip acceptrate column with % (numpy does not know how to handle percents)
supernova_fitting = np.loadtxt(supernova_fitting_data_filename, usecols=(0, 1, 2, 3), skiprows=200)
supernova_fitting_chi2 = supernova_fitting[:, 0]
supernova_fitting_h = supernova_fitting[:, 1]
supernova_fitting_H0 = supernova_fitting_h * 100
supernova_fitting_OmegaM = supernova_fitting[:, 2]
supernova_fitting_OmegaK = supernova_fitting[:, 3]
supernova_fitting_OmegaLambda = 1 - supernova_fitting_OmegaM - supernova_fitting_OmegaK

# print(f"{supernova_fitting=}")
chi2_min = np.min(supernova_fitting_chi2)
chi2_min_index = np.argmin(supernova_fitting_chi2)
print(f"Best chi^2 score: {chi2_min}. Should be of order N={N}, chi^2/N = {chi2_min / N:.3f}.")

best_h = supernova_fitting_h[chi2_min_index]
best_H0 = best_h * 100
best_OmegaM = supernova_fitting_OmegaM[chi2_min_index]
best_OmegaK = supernova_fitting_OmegaK[chi2_min_index]
best_OmegaLambda = 1 - best_OmegaM - best_OmegaK

print()  # empty line
print("Best-fit parameters")
print(f"Best h: {best_h}")
print(f"Best H0: {best_H0} km/s/Mpc")
print(f"Best OmegaM: {best_OmegaM}")
print(f"Best OmegaK: {best_OmegaK}")
print(f"Best OmegaLambda: {best_OmegaLambda}")

# sigma constraints from https://www.reid.ai/2012/09/chi-squared-distribution-table-with.html
supernova_fitting_1sigma_index = np.argwhere(supernova_fitting_chi2 - chi2_min < 3.53)
supernova_fitting_chi2_1sigma = supernova_fitting_chi2[supernova_fitting_1sigma_index]
supernova_fitting_h_1sigma = supernova_fitting_h[supernova_fitting_1sigma_index]
supernova_fitting_OmegaM_1sigma = supernova_fitting_OmegaM[supernova_fitting_1sigma_index]
supernova_fitting_OmegaK_1sigma = supernova_fitting_OmegaK[supernova_fitting_1sigma_index]
supernova_fitting_OmegaLambda_1sigma = 1 - supernova_fitting_OmegaM_1sigma - supernova_fitting_OmegaK_1sigma

supernova_fitting_2sigma_index = np.argwhere(supernova_fitting_chi2 - chi2_min < 8.02)
supernova_fitting_chi2_2sigma = supernova_fitting_chi2[supernova_fitting_2sigma_index]
supernova_fitting_h_2sigma = supernova_fitting_h[supernova_fitting_2sigma_index]
supernova_fitting_OmegaM_2sigma = supernova_fitting_OmegaM[supernova_fitting_2sigma_index]
supernova_fitting_OmegaK_2sigma = supernova_fitting_OmegaK[supernova_fitting_2sigma_index]
supernova_fitting_OmegaLambda_2sigma = 1 - supernova_fitting_OmegaM_2sigma - supernova_fitting_OmegaK_2sigma

# flat universe
n_points_flat = 100
supernova_fitting_OmegaK_flat = np.zeros(n_points_flat)
supernova_fitting_OmegaM_flat = np.linspace(0, 1, n_points_flat)
supernova_fitting_OmegaLambda_flat = 1 - supernova_fitting_OmegaM_flat - supernova_fitting_OmegaK_flat

"""
Plotting
"""
# H_prime
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
H_unit = (100 * constants.km / constants.s) / constants.Mpc
ax.plot(cosmology_df.index, cosmology_df["Hp"] / H_unit, label=r"$\mathcal{H}(x)$")
ax.hlines(
    y=constants.input_cosmology_params_dict["h"],
    xmin=cosmology_df.index[0],
    xmax=cosmology_df.index[-1],
    linestyles="dashed",
    label="$h$",
)
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\mathcal{H}$")
ax.set_title(r"$\mathcal{H}(x)$ $\left( \frac{100 km/s}{Mpc} \right)$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "H_prime_of_x.png")
if show_plots:
    plt.show()

# eta
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(cosmology_df.index, cosmology_df["eta"] / constants.Mpc)
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\eta$")
ax.set_title(r"$\eta(x)$ $(Mpc)$")
ax.grid()
fig.savefig(plots_folder / "eta_of_x.png")
if show_plots:
    plt.show()

# conversion factor seconds to gigayears
gigayears = 10**9 * 365.24 * 24 * 60 * 60 * constants.s

# conformal time
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(cosmology_df.index, cosmology_df["eta"] / (constants.c * gigayears))
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\eta/c$")
ax.set_title(r"Conformal time $\frac{\eta(x)}{c}$ $(10^9 \, yrs)$")
ax.grid()
fig.savefig(plots_folder / "conformal_time.png")
if show_plots:
    plt.show()

# cosmic time
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(cosmology_df.index, cosmology_df["t"] / gigayears)
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$t$")
ax.set_title(r"Cosmic time $t$ $(10^9 \, yrs)$")
ax.grid()
fig.savefig(plots_folder / "cosmic_time.png")
if show_plots:
    plt.show()

# cosmic vs. conformal time
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(cosmology_df.index, cosmology_df["eta"] / (constants.c * gigayears), label=r"Conformal time $\eta(x)/c$")
(line,) = ax.plot(cosmology_df.index, cosmology_df["t"] / gigayears, label=r"Cosmic time $t$")
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$time$")
ax.set_title(r"Cosmic vs Conformal time $(10^9 \, yrs)$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "cosmic_vs_conformal_time.png")
if show_plots:
    plt.show()

print()  # blank line
print(f"Calculated current age of the Universe, cosmic time: {cosmology_df['t'][0] / gigayears:.3f} * 10^9 yrs")
print("Actual age of the Universe should be 14.3 billion years ish")

print(
    f"Calculated current age of the Universe, conformal time: {cosmology_df['eta'][0] / (constants.c * gigayears):.3f} * 10^9 yrs"
)

# eta(x) H_prime / c
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(cosmology_df.index, cosmology_df["eta"] * cosmology_df["Hp"] / constants.c)
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$(\eta(x) \mathcal{H})/c$")
ax.set_title(r"$\frac{\eta(x) \mathcal{H}}{c}$")
ax.grid()
fig.savefig(plots_folder / "etaHp_over_c_of_x.png")
if show_plots:
    plt.show()

# omega_i components
omega_R = cosmology_df["Omega_gamma"] + cosmology_df["Omega_Nu"]
omega_M = cosmology_df["Omega_B"] + cosmology_df["Omega_CDM"]
omega_L = cosmology_df["Omega_Lambda"]
total_density = omega_R + omega_M + omega_L  # should be 1

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(cosmology_df.index, omega_R, label=r"$\Omega_{relativistic} = \Omega_r + \Omega_{\nu}$")
ax.plot(cosmology_df.index, omega_M, label=r"$\Omega_{matter} = \Omega_b + \Omega_{CDM}$")
ax.plot(cosmology_df.index, omega_L, label=r"$\Omega_{\Lambda}$")
ax.plot(cosmology_df.index, total_density, label=r"$\sum{\Omega_i} = 1$", linestyle="dashed")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\Omega_i$")
ax.set_title(r"$\Omega_i (x)$")
# ax.legend()
ax.grid()

# from https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width, box.height * 0.95])

# Put a legend below current axis
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.20), fancybox=True, shadow=True, ncol=5)

fig.savefig(plots_folder / "Omega_i_of_x.png")
if show_plots:
    plt.show()

# luminosity distance vs redshift
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
eps = 1e-11  # avoid division by zero
(line,) = ax.plot(z_array, calculated_dL / (z_array + eps), label="Calculated distances $d_L$")
ax.errorbar(
    supernovadata_z,
    supernovadata_dL / (supernovadata_z + eps),
    yerr=supernovadata_error / (supernovadata_z + eps),
    fmt="r.",
    capsize=2,
    label="Betoule et al. 2014",
)
# ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim(3.4, 8.1)
ax.set_xlabel(r"$z$")
ax.set_ylabel(r"$d_L / z$")
ax.set_title(r"Luminosity distance $d_L / z$ $(Gpc)$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "luminosity_distance.png")
if show_plots:
    plt.show()

# scatter plot of confidence regions
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(
    supernova_fitting_OmegaM_flat,
    supernova_fitting_OmegaLambda_flat,
    "k--",
    label=r"Flat Universe",
)
ax.scatter(
    supernova_fitting_OmegaM_2sigma,
    supernova_fitting_OmegaLambda_2sigma,
    label=r"2$\sigma$ constraint",
)
ax.scatter(
    supernova_fitting_OmegaM_1sigma,
    supernova_fitting_OmegaLambda_1sigma,
    label=r"1$\sigma$ constraint",
)
# ax.set_yscale("log")
# ax.set_xscale("log")
ax.set_xlabel(r"$\Omega_M$")
ax.set_ylabel(r"$\Omega_{\Lambda}$")
ax.set_title(r"Parameter constraints from supernova fitting")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "supernovafitting_confidence_regions.png")
if show_plots:
    plt.show()

# assuming posterior is gaussian
parameter_list = {
    "H0": supernova_fitting_H0,
    "OmegaM": supernova_fitting_OmegaM,
    "OmegaK": supernova_fitting_OmegaK,
    "OmegaLambda": supernova_fitting_OmegaLambda,
}
for param_name, parameter_array in parameter_list.items():
    print(f"{param_name}")
    current_std = np.std(parameter_array)
    print(f"STD: {current_std:.5f}")
    current_mean = np.mean(parameter_array)
    print(f"Mean: {current_mean:.5f}")
    current_variance = np.var(parameter_array)
    current_sigma = np.sqrt(current_variance)
    num_bins = 20

    # +0.1 looks slightly better in plot
    param_space_array = np.linspace(np.min(parameter_array), np.max(parameter_array) + 0.1, 1000)
    gaussian_array = (1 / (current_sigma * np.sqrt(2 * np.pi))) * np.exp(
        -(1 / 2) * (param_space_array - current_mean) ** 2 / (current_variance)
    )

    # histogram of parameter selections
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    result = ax.hist(
        parameter_array,
        bins=num_bins,
        density=True,
    )

    # plot equivalent gaussian distribution
    ax.plot(param_space_array, gaussian_array)

    # make gaussian plot with seaborn
    # sns.histplot(parameter_array, bins=num_bins, kde=True, color="blue", common_norm=True, ax=ax)

    ax.axvline(parameter_array[chi2_min_index], color="black", linestyle="dashed", label="Best-fit value")

    ax.set_xlabel(f"{param_name}")
    ax.set_ylabel(r"Bin count / (Total count * bin width)")
    ax.set_title(r"Parameter histogram")
    # ax.grid()
    ax.legend()
    fig.savefig(plots_folder / f"{param_name}_histogram.png")
    if show_plots:
        plt.show()
