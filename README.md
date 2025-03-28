# Cosmology II
## The Cosmic Microwave Background and the Large Scale Structure of our Universe in theory and practice

This repository is for a specific solution to the numerical project in the course AST5220 "Cosmology II" at ITA University of Oslo. It is based on the C++ templates (and general repository) at [HAWinther/AST5220-Cosmology](https://github.com/HAWinther/AST5220-Cosmology).

This project is for making an Einstein-Boltzmann solver (a CAMB like code). This code will ultimately lead to matter and CMB power spectra, after solving theoretical equations. 

For the current version of the course C++ is the main language, this project focuses on the C++ code.

# Fortran

Not using Fortran. PDFs outlining what to do for each milestone kept in Fortran folder.

# Python?

Simple templates for project with Python in `python_template`. Strongly recommended to not use Python though.

# Website
All relevant information about the project and the different milestones can be found on this [website](https://cmb.wintherscoming.no/).

# Compiling

Compile target `cmb` using [ make cmb ] in terminal, as expected. The [GSL library](ftp://ftp.gnu.org/gnu/gsl/) is required. See below for instructions.

# Running
The compiled binary `cmb` used to generate output for each Milestone is copied by me to the respective Milestone X folder before delivery. The binary expects to be in the project root folder, and must be moved there before running. Use `$ ./cmb` to run.

If you compile and run at the same time, it should generate `cmb` in the project root folder and run with the most up to date source code automatically.

TODO: Change instructions below as changes are made

The code runs from Main.cpp and then proceeds to go through the different milestones one by one untill we have the CMB power spectra in the end. If you get it compiled then run it as [ ./cmb ]. It will crash with "Error: Spline eta has not been created". That is fine, it's one of your task to implement this. 

See Examples.cpp - and run the examples as [ make examples; ./examples ; ] - for examples on how to make splines, solve ODEs and the functionality of the stuff supplied with this template.

For the last milestone you need to compute spherical bessel functions (see the function j\_ell(n, x) in Utils.cpp). There are several options for this: if you have a C++17 compiler (use -std=c++17 instead of c++11 in the Makefile) then you can use provided by the C++ standard std::sph\_bessel(n, x). Otherwise the GSL library provides a function gsl\_sf\_bessel\_jl(n, x) for this that we use in the template. This implementation has problems with underflow for small x and large n, but we correct this using known asymptotical expressions. The last option is to use another library like for example the [Complex Bessel](https://github.com/joeydumont/complex_bessel) library.

# How to install GSL

- Use homebrew, mac, `brew install gsl`

- Makefile changed to point to library

```
INC = -I/opt/homebrew/opt/gsl/include
LIBS = -L/opt/homebrew/opt/gsl/lib -lgsl -lgslcblas
```


