# feos
feos is a Fortran library for the calculation of the residual Helmholtz Energy
for mixtures and the thermodynamic properties that are possible to calculate
with it, as well as fluid equilibria.

# What it provides?

## Cubic equations of state
feos includes the following three cubic equations of state, but it also
provides and easy to use interface to extend it with other equations of state.
Taking advantage of a heavy modular design, where each equation is designed as
an object pointing to it's corresponding parameters equations.

The included Cubic EOS are:

- Peng Robinson
- Soave-Redlich-Kwong
- RKPR

### Mixing rules
- Classic Van der Waals
    - Constant $k_{ij}$ and $l_{ij}$
    - Exponentialy dependant of temperature $k_{ij}$
- Cubic mixing rules
    - Constant $k_{ij}$ and $l_{ij}$
    - Exponentialy dependant of temperature $k_{ij}$

# Motivation
A lot of Fortran code is used for thermodynamics calculation, either being
due to the need to work with legacy code or for taking advantage of the
speed performance that it provides.

feos intends to both give an easy to use library for any user intersted in
Cubic EoS usage.

# Installation
This program is intended to be used with fortran-lang's 
[Fortran Package Manager (`fpm`)](https://fpm.fortran-lang.org/en/index.html)

To build it just `cd` into the main directory and run

```sh
fpm build
```

The main program can also be run using:

```
fpm run
```

Or running the executable on the build dir.
