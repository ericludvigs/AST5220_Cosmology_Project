# Basic units (here we use SI)
m           = 1.0                         # Length (in meters)
s           = 1.0                         # Time (in seconds)
kg          = 1.0                         # Kilo (in kilos)
K           = 1.0                         # Temperature (in Kelvins)

# Derived units
km          = 1e3 * m                     # Kilometers
N           = kg*m/(s*s)                  # Newton
J           = N*m                         # Joule
W           = J/s                         # Watt
Mpc         = 3.08567758e22 * m           # Megaparsec
eV          = 1.60217653e-19 * J          # Electronvolt

# Physical constants
k_b         = 1.38064852e-23 * J/K        # Bolzmanns constant
m_e         = 9.10938356e-31 * kg         # Mass of electron
m_H         = 1.6735575e-27 * kg          # Mass of hydrogen atom
c           = 2.99792458e8 * m/s          # Speed of light
G           = 6.67430e-11 * N*m*m/(kg*kg) # Gravitational constant
hbar        = 1.054571817e-34 * J*s       # Reduced Plancks constant
sigma_T     = 6.6524587158e-29 * m*m      # Thomas scattering cross-section
lambda_2s1s = 8.227 / s                   # Transition time between 2s and 1s in Hydrogen
H0_over_h   = 100 * km/s/Mpc              # H0 / h
epsilon_0   = 13.605693122994 * eV        # Ionization energy for the ground state of hydrogen
xhi0        = 24.587387 * eV              # Ionization energy for neutral Helium
xhi1        = 4.0 * epsilon_0             # Ionization energy for singly ionized Helium
