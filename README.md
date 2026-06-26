# MATLAB scripts for figure generation

This repository contains the MATLAB scripts used to generate the figures for

> **Contrasting Optical Properties of Different Arctic Sea Ice Regimes**

## Contents

- **1 orthomosaics/**
  - `plot_orthomosaic.m`
- **2 energy budget/**
  - `plot_heat_budget_clean.m`
- **3 map and forcing/**
  - `map_and_forcing.m`
- **data/**
  - Input datasets required by the scripts

## Output

- `orthomosaic_1B.tif`
- `heat_budget_components.png`
- `map_forcing.png`

## Data

The `data` directory contains:

- `AMSR2_SIC/`
- `era5/`
- `Polarstern/`
- `radiation stations/`
- `SIMBA/`
- `SVP buoys/`

The scripts assume this directory structure. Update the file paths if the repository is stored in a different location.

## Requirements

- MATLAB
- Mapping Toolbox
- TEOS-10 GSW Toolbox
- M_Map toolbox
- Optional: Crameri Scientific Colour Maps (`oslo`, `buda`)

---

Developed for:

**Contrasting Optical Properties of Different Arctic Sea Ice Regimes**

F. Zimmer, R. Tao, T. Blei, E. Salganik, M. Nicolaus
