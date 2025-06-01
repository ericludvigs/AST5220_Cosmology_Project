import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append("..")
import python_methods.constants as constants

show_plots = False

recombination_data_filename = Path("../results/recombination.csv")
plots_folder = Path("Plots/")

recombination_df = pd.read_csv(recombination_data_filename, sep=";", header=0, index_col=0)
recombination_df.columns = recombination_df.columns.str.strip()
print(recombination_df)

# make an array with redshifts
z_array = constants.redshift(recombination_df.index.to_numpy())
# print(f"{z_array=}")

# first index where we go below 0.99 is where we switched approximations
saha_regime_switch_index = int(np.argwhere(recombination_df["Xe"] < 0.99)[0][0])  # numpy wants double extraction for some reason
saha_regime_switch_x_val = recombination_df.index[saha_regime_switch_index]

# kind of arbitrary when we set recombination
recombination_condition = 0.50
recombination_index = int(np.argwhere(recombination_df["Xe"] < recombination_condition)[0][0])
recombination_x_val = recombination_df.index[recombination_index]

recombination_year = recombination_df["t(x)"].loc[recombination_x_val]
recombination_year = recombination_year / constants.years

print(f"Switched from Saha approx at index={saha_regime_switch_index}, with x = {saha_regime_switch_x_val:.3f}")
print()
print(f"Defining recombination when free electron fraction drops below Xe < {recombination_condition:.1f}")
print(f"Recombination at index={recombination_index}, with x = {recombination_x_val:.3f}")
formatted_year = f"{recombination_year:,.3f}".replace(",", " ")
print(f"Recombination time: {formatted_year} years")
print(f"Redshift: {z_array[saha_regime_switch_index]:.3f}")

# tau + derivs plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.axvline(saha_regime_switch_x_val, color="purple", linestyle="dashed", label="Saha regime switch", linewidth=1)
(line,) = ax.plot(recombination_df.index, recombination_df["tau"], label=r"Optical depth $\tau(x)$")
(line,) = ax.plot(recombination_df.index, -recombination_df["dtau_dx"], label=r"$-\tau'(x)$")
(line,) = ax.plot(recombination_df.index, recombination_df["ddtau_ddx"], label=r"$\tau''(x)$")
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\tau$")
ax.set_title(r"Optical depth function and its derivatives [dimensionless]")
ax.set_ylim(10**-8, 10**8)
ax.grid()
ax.legend()
fig.savefig(plots_folder / "tau_and_derivs.png")
if show_plots:
    plt.show()

# X_e plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.axvline(saha_regime_switch_x_val, color="purple", linestyle="dashed", label="Saha regime switch", linewidth=1)
(line,) = ax.plot(recombination_df.index, recombination_df["Xe"], label=r"$X_e(x)$")
# TODO: get X_e with simple reion
# (line,) = ax.plot(recombination_df.index, -recombination_df["dtau_dx"], label=r"$-\tau'(x)$")
ax.set_yscale("log")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\tau$")
ax.set_title(r"Fractional electron density $X_e = n_e / n_H$ [dimensionless]")
# ax.set_ylim(10**-8, 10**8)
ax.grid()
ax.legend()
fig.savefig(plots_folder / "Xe.png")
if show_plots:
    plt.show()

# g_tilde plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.axvline(saha_regime_switch_x_val, color="purple", linestyle="dashed", label="Saha regime switch", linewidth=1)
(line,) = ax.plot(recombination_df.index, recombination_df["g_tilde"], label=r"$\tilde{g}(x)$")
# TODO: get g_tilde with simple reion?
# (line,) = ax.plot(recombination_df.index, -recombination_df["dtau_dx"], label=r"$-\tau'(x)$")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\tilde{g}$")
ax.set_title(r"Visibility function $\tilde{g}(x)$ [PDF]")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "g_tilde.png")
if show_plots:
    plt.show()

# g_tilde deriv plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.axvline(saha_regime_switch_x_val, color="purple", linestyle="dashed", label="Saha regime switch", linewidth=1)
(line,) = ax.plot(recombination_df.index, recombination_df["dg_tilde_dx"], label=r"$\tilde{g}'(x)$")
# TODO: get g_tilde with simple reion?
# (line,) = ax.plot(recombination_df.index, -recombination_df["dtau_dx"], label=r"$-\tau'(x)$")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\tilde{g}'$")
ax.set_title(r"Visibility function derivative $\tilde{g}'(x)$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "g_tilde_deriv.png")
if show_plots:
    plt.show()

# g_tilde double deriv plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.axvline(saha_regime_switch_x_val, color="purple", linestyle="dashed", label="Saha regime switch", linewidth=1)
(line,) = ax.plot(recombination_df.index, recombination_df["ddg_tilde_ddx"], label=r"$\tilde{g}''(x)$")
# TODO: get g_tilde with simple reion?
# (line,) = ax.plot(recombination_df.index, -recombination_df["dtau_dx"], label=r"$-\tau'(x)$")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$\tilde{g}''$")
ax.set_title(r"Visibility function double derivative $\tilde{g}''(x)$")
ax.grid()
ax.legend()
fig.savefig(plots_folder / "g_tilde_double_deriv.png")
if show_plots:
    plt.show()
