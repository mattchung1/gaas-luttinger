# GaAs Luttinger Hamiltonian Band Structure and DOS

![Description of Image](figures/demo_screenshot.png)

**Check figures folder for all .pngs of output figures.**

This project computes:

- 6×6 Luttinger Hamiltonian for GaAs
- Band structure along high-symmetry k-directions
- 3D iso-energy surfaces for heavy- and light-hole bands
- Density of states using both Monte Carlo sampling and a deterministic grid

## Requirements

- MATLAB R20xx
- No toolboxes required (only base MATLAB functions)

## Usage

```matlab
cd src
main
```

## The script will:

- Plot band structure vs k
- Plot HH and LH iso-surfaces
- Plot DOS for HH and LH bands with a √E trendlin
