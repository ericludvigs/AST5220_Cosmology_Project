from pathlib import Path

import constants
import matplotlib.pyplot as plt
import pandas as pd

show_plots = False

cosmology_data_filename = Path("../results/cosmology.csv")
plots_folder = Path("Plots/")

cosmology_df = pd.read_csv(cosmology_data_filename, sep="|", header=0, index_col=0)
cosmology_df.columns = cosmology_df.columns.str.strip()
print(cosmology_df)

# H_prime
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
# ((100*constants.km/constants.s)/constants.Mpc)
line, = ax.plot(cosmology_df.index, cosmology_df["Hp"]*constants.Mpc/(100*constants.km/constants.s))
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\mathcal{H}$")
ax.set_title(r"$\mathcal{H}(x)$ $\left( \frac{100 km/s}{Mpc} \right)$")
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
ax.legend()
fig.savefig(plots_folder / "cosmic_vs_conformal_time.png")
if show_plots:
    plt.show()

print(f"Current age of the universe: {cosmology_df["t"][0]/gigayears:.3f} * 10^9 yrs ")
print("Actual age of the universe should be 14.3 billion years ish")

# eta(x) H_prime / c
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
line, = ax.plot(cosmology_df.index, cosmology_df["eta"]*cosmology_df["Hp"]/constants.c)
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$(\eta(x) \mathcal{H})/c$")
ax.set_title(r"$\frac{\eta(x) \mathcal{H}}{c}$")
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
ax.plot(cosmology_df.index, total_density, label=r"$\sum{\Omega_i} = 1$")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\Omega_i$")
ax.set_title(r"$\Omega_i (x)$")
ax.legend()
fig.savefig(plots_folder / "Omega_i_of_x.png")
if show_plots:
    plt.show()
