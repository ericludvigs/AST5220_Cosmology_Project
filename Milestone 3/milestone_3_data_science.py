import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.append("..")
import python_methods.constants as constants

show_plots = False

perturbation_small_scales_data_filename = Path("../results/perturbations_k0_001.csv")
plots_folder = Path("Plots/")

perturbation_small_scales_df = pd.read_csv(perturbation_small_scales_data_filename, sep=";", header=0, index_col=0)
perturbation_small_scales_df.columns = perturbation_small_scales_df.columns.str.strip()
print(perturbation_small_scales_df)

perturbation_medium_scales_data_filename = Path("../results/perturbations_k0_01.csv")

perturbation_medium_scales_df = pd.read_csv(perturbation_medium_scales_data_filename, sep=";", header=0, index_col=0)
perturbation_medium_scales_df.columns = perturbation_medium_scales_df.columns.str.strip()
print(perturbation_medium_scales_df)

perturbation_large_scales_data_filename = Path("../results/perturbations_k0_1.csv")

perturbation_large_scales_df = pd.read_csv(perturbation_large_scales_data_filename, sep=";", header=0, index_col=0)
perturbation_large_scales_df.columns = perturbation_large_scales_df.columns.str.strip()
print(perturbation_large_scales_df)


# make an array with redshifts
z_array = constants.redshift(perturbation_small_scales_df.index.to_numpy())
# print(f"{z_array=}")

# delta gamma, delta b and cdm plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(perturbation_small_scales_df.index, perturbation_small_scales_df["delta_cdm"], label=r"$k = 0.001/$Mpc", c="blue")
(line,) = ax.plot(perturbation_small_scales_df.index, 4 * perturbation_small_scales_df["Theta_0"], linestyle="dotted", c="blue")
(line,) = ax.plot(perturbation_small_scales_df.index, perturbation_small_scales_df["delta_b"], linestyle="dashed", c="blue")

(line,) = ax.plot(
    perturbation_medium_scales_df.index, perturbation_medium_scales_df["delta_cdm"], label=r"$k = 0.01/$Mpc", c="orange"
)
(line,) = ax.plot(perturbation_medium_scales_df.index, 4 * perturbation_medium_scales_df["Theta_0"], linestyle="dotted", c="orange")
(line,) = ax.plot(perturbation_medium_scales_df.index, perturbation_medium_scales_df["delta_b"], linestyle="dashed", c="orange")

(line,) = ax.plot(perturbation_large_scales_df.index, perturbation_large_scales_df["delta_cdm"], label=r"$k = 0.1/$Mpc", c="green")
(line,) = ax.plot(perturbation_large_scales_df.index, 4 * perturbation_large_scales_df["Theta_0"], linestyle="dotted", c="green")
(line,) = ax.plot(perturbation_large_scales_df.index, perturbation_large_scales_df["delta_b"], linestyle="dashed", c="green")

ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"Perturbation magnitude")
ax.set_title(r"$\delta_{\gamma}, \delta_{CDM}, \delta_b$")
# ax.set_ylim(10**-8, 10**8)
ax.grid()
ax.legend()
fig.savefig(plots_folder / "delta_gamma_delta_b_delta_cdm_plot.png")
if show_plots:
    plt.show()

# velocity b and cdm plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(perturbation_small_scales_df.index, perturbation_small_scales_df["v_cdm"], label=r"$k = 0.001/$Mpc", c="blue")
(line,) = ax.plot(perturbation_small_scales_df.index, perturbation_small_scales_df["v_b"], linestyle="dashed", c="blue")

(line,) = ax.plot(perturbation_medium_scales_df.index, perturbation_medium_scales_df["v_cdm"], label=r"$k = 0.01/$Mpc", c="orange")
(line,) = ax.plot(perturbation_medium_scales_df.index, perturbation_medium_scales_df["v_b"], linestyle="dashed", c="orange")

(line,) = ax.plot(perturbation_large_scales_df.index, perturbation_large_scales_df["v_cdm"], label=r"$k = 0.1/$Mpc", c="green")
(line,) = ax.plot(perturbation_large_scales_df.index, perturbation_large_scales_df["v_b"], linestyle="dashed", c="green")

ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"Perturbation magnitude")
ax.set_title(r"$v_{CDM}, v_b$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "v_b_v_cdm_plot.png")
if show_plots:
    plt.show()

# theta_0 plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(perturbation_small_scales_df.index, perturbation_small_scales_df["Theta_0"], label=r"$k = 0.001/$Mpc", c="blue")
(line,) = ax.plot(
    perturbation_medium_scales_df.index, perturbation_medium_scales_df["Theta_0"], label=r"$k = 0.01/$Mpc", c="orange"
)
(line,) = ax.plot(perturbation_large_scales_df.index, perturbation_large_scales_df["Theta_0"], label=r"$k = 0.1/$Mpc", c="green")
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"Perturbation magnitude")
ax.set_title(r"$\Theta_0$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "theta_0_plot.png")
if show_plots:
    plt.show()

# theta_1 plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(perturbation_small_scales_df.index, perturbation_small_scales_df["Theta_1"], label=r"$k = 0.001/$Mpc", c="blue")
(line,) = ax.plot(
    perturbation_medium_scales_df.index, perturbation_medium_scales_df["Theta_1"], label=r"$k = 0.01/$Mpc", c="orange"
)
(line,) = ax.plot(perturbation_large_scales_df.index, perturbation_large_scales_df["Theta_1"], label=r"$k = 0.1/$Mpc", c="green")
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"Perturbation magnitude")
ax.set_title(r"$\Theta_1$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "theta_1_plot.png")
if show_plots:
    plt.show()

# phi plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(perturbation_small_scales_df.index, perturbation_small_scales_df["Phi"], label=r"$k = 0.001/$Mpc", c="blue")
(line,) = ax.plot(perturbation_medium_scales_df.index, perturbation_medium_scales_df["Phi"], label=r"$k = 0.01/$Mpc", c="orange")
(line,) = ax.plot(perturbation_large_scales_df.index, perturbation_large_scales_df["Phi"], label=r"$k = 0.1/$Mpc", c="green")
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"Perturbation magnitude")
ax.set_title(r"$\Phi$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "phi_plot.png")
if show_plots:
    plt.show()

# phi+psi plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
(line,) = ax.plot(
    perturbation_small_scales_df.index,
    perturbation_small_scales_df["Phi"] + perturbation_small_scales_df["Psi"],
    label=r"$k = 0.001/$Mpc",
    c="blue",
)
(line,) = ax.plot(
    perturbation_medium_scales_df.index,
    perturbation_medium_scales_df["Phi"] + perturbation_small_scales_df["Psi"],
    label=r"$k = 0.01/$Mpc",
    c="orange",
)
(line,) = ax.plot(
    perturbation_large_scales_df.index,
    perturbation_large_scales_df["Phi"] + perturbation_small_scales_df["Psi"],
    label=r"$k = 0.1/$Mpc",
    c="green",
)
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"Perturbation magnitude")
ax.set_title(r"$\Phi + \Psi$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "phi_plus_psi_plot.png")
if show_plots:
    plt.show()
