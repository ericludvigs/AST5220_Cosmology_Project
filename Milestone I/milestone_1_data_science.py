from pathlib import Path

import constants
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

show_plots = False

def redshift(x):
    a = np.exp(x)
    z = 1.0/a - 1
    return z

cosmology_data_filename = Path("../results/cosmology.csv")
supernova_data_filename = Path("../data/supernovadata.txt").resolve()
plots_folder = Path("Plots/")

cosmology_df = pd.read_csv(cosmology_data_filename, sep="|", header=0, index_col=0)
cosmology_df.columns = cosmology_df.columns.str.strip()
print(cosmology_df)

# make an array with redshifts, reverse it and cut it down to match with the input data
z_array = redshift(cosmology_df.index.to_numpy())[::-1]
#print(f"{z_array=}")

supernovadata = np.loadtxt(supernova_data_filename)
supernovadata_z = supernovadata[:,0]
supernovadata_dL = supernovadata[:,1]
supernovadata_error = supernovadata[:,2]
print(f"{supernovadata=}")
print(f"{supernovadata_z=}")

z_array = z_array[z_array < np.max(supernovadata_z)+0.02]
print(f"{z_array=}")

# get the matching lum distances we calculated, also in Gpc
# need to reverse, and only select the cosmologically "closest" elements to match with measured distances
calculated_dL = (cosmology_df["d_L"].to_numpy()[::-1]/constants.Gpc)[:len(z_array)]

print(f"{calculated_dL=}")

# H_prime
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
H_unit = ((100*constants.km/constants.s)/constants.Mpc)
ax.plot(cosmology_df.index, cosmology_df["Hp"]/H_unit, label=r"$\mathcal{H}(x)$")
ax.hlines(y=constants.input_cosmology_params_dict["h"], xmin=cosmology_df.index[0], xmax=cosmology_df.index[-1],
          linestyles="dashed", label="$h$")
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\mathcal{H}$")
ax.set_title(r"$\mathcal{H}(x)$ $\left( \frac{100 km/s}{Mpc} \right)$")
ax.grid()
ax.legend()
fig.savefig(plots_folder/"H_prime_of_x.png")
if show_plots:
    plt.show()

# eta
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
line, = ax.plot(cosmology_df.index, cosmology_df["eta"]/constants.Mpc)
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
line, = ax.plot(cosmology_df.index, cosmology_df["eta"]/(constants.c * gigayears))
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
line, = ax.plot(cosmology_df.index, cosmology_df["t"]/gigayears)
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
line, = ax.plot(cosmology_df.index, cosmology_df["eta"]/(constants.c * gigayears), label=r"Conformal time $\eta(x)/c$")
line, = ax.plot(cosmology_df.index, cosmology_df["t"]/gigayears, label=r"Cosmic time $t$")
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$time$")
ax.set_title(r"Cosmic vs Conformal time $(10^9 \, yrs)$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "cosmic_vs_conformal_time.png")
if show_plots:
    plt.show()

print() # blank line
print(f"Calculated current age of the universe, cosmic time: {cosmology_df["t"][0]/gigayears:.3f} * 10^9 yrs")
print("Actual age of the universe should be 14.3 billion years ish")

print(f"Calculated current age of the universe, conformal time: {cosmology_df["eta"][0]/(constants.c*gigayears):.3f} * 10^9 yrs")

# eta(x) H_prime / c
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
line, = ax.plot(cosmology_df.index, cosmology_df["eta"]*cosmology_df["Hp"]/constants.c)
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
total_density = omega_R + omega_M + omega_L # should be 1

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(cosmology_df.index, omega_R, label=r"$\Omega_{relativistic} = \Omega_r + \Omega_{\nu}$")
ax.plot(cosmology_df.index, omega_M, label=r"$\Omega_{matter} = \Omega_b + \Omega_{CDM}$")
ax.plot(cosmology_df.index, omega_L, label=r"$\Omega_{\Lambda}$")
ax.plot(cosmology_df.index, total_density, label=r"$\sum{\Omega_i} = 1$", linestyle="dashed")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\Omega_i$")
ax.set_title(r"$\Omega_i (x)$")
#ax.legend()
ax.grid()

# from https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot
# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0,
                 box.width, box.height * 0.95])

# Put a legend below current axis
ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.20),
          fancybox=True, shadow=True, ncol=5)

fig.savefig(plots_folder / "Omega_i_of_x.png")
if show_plots:
    plt.show()

# luminosity distance vs redshift
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
line, = ax.plot(z_array, calculated_dL/z_array, label="Calculated distances $d_L$")
ax.errorbar(
    supernovadata_z, supernovadata_dL/supernovadata_z, yerr=supernovadata_error, fmt="r.", capsize=2, label="Betoule et al. 2014"
    )
#ax.set_yscale("log")
ax.set_xscale("log")
ax.set_xlabel(r"$z$")
ax.set_ylabel(r"$d_L$")
ax.set_title(r"Luminosity distance $d_L$ $(Gpc)$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "luminosity_distance.png")
if show_plots:
    plt.show()
