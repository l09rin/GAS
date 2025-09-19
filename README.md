# GAS
**G**PU **A**ccelerated numerical **S**imulations in CUDA–C++.

Currently supports **constant-energy Molecular Dynamics (MD)** of **Lennard–Jones** particle systems.
Originally based on the [MoDyCa](https://github.com/l09rin/MoDyCa) project (serial code running on CPU).

---

## Getting started

### Clone & Compile
```bash
git clone git@github.com:l09rin/GAS.git
cd GAS
make all
```
Before building the main executable, the build system runs a **CUDA compiler test**.
The GPU architecture is **auto-detected**, but if detection fails, you can specify the compute capability manually, e.g. for Ampere GPUs (CC=80):

```bash
make CUDA_CC=80 all
```

### Install (optional)
Move the executable into your `PATH`:
```bash
mv gas.exe ~/.local/bin/
```

---

## Usage
The program accepts two optional arguments:
1. `-in <file>` : Path to input parameter file.
2. `-seed <int>` : Random seed for initialization.

Example:
```bash
./gas.exe -in input.dat -seed 8000008
```

This will create an **output directory** containing results.  

👉 Example inputs are provided in `examples/lennard_jones/`.

---

## Input

The simulation is controlled by a plain-text input file (example: [`examples/lennard_jones/input_2D.gpu`](examples/lennard_jones/input_2D.gpu)).
Each line consists of a **keyword** followed by one or more parameters.  
Comments start with `#`.  

### General Settings
- **`DIMENSION <2|3>`**  :  Set simulation dimensionality (default = 3).  
- **`DEVICE <CPU|GPU>`**  :  Select compute device.  
- **`PARTICLE_CLASS <particle_3D|particle_2D>`**  :  Choose particle class. Default = `particle_3D`.  
- **`ENSEMBLE <NVE>`**  :  Currently only constant-energy (NVE) supported.  

---

### Physical Units
- **`SIGMA <float>`**  :  Length unit (default = 1, not explicitly used in code).  
- **`EPSILON <float>`**  :  Energy unit (default = 1, not explicitly used).  
- **`PARTICLE_MASS <float>`**  :  Particle mass (default = 1).  
- **`SIDE <float>`**  :  Size of a cubic/square box in units of `SIGMA`.  
- **`SIDES <float> <float>`**  :  For anisotropic boxes: specify each dimension separately.  
- **`KbT <float>`**  :  Thermal energy in units of `EPSILON` (currently not used).  

---

### Initial Configuration
Two options:
1. **From file**  
   ``` 
   INITIAL_CONFIGURATION file <xyz|sph|lmp> <file_name> <energy_per_particle>
   ```
   - `xyz` : _XYZ_-like format
   - `sph` : sph data file as dumped by this [Event Driven Molecular Dynamics code](https://github.com/FSmallenburg/EDMD)
   - `lmp` : LAMMPS data file
2. **Generated configuration**  
   ```
   INITIAL_CONFIGURATION <random|fcc|bcc> <number_of_particles> <energy_per_particle>
   ```
   - `random` : random placement  
   - `fcc` : face-centered cubic lattice  
   - `bcc` : body-centered cubic lattice  

- **`INITIAL_VELOCITIES_BY_FILE <flag> <file>`**  :  If `flag=1`, load initial velocities from `file`.  

---

### Interaction Potential
- **`LENNARD_JONES <cutoff> <shift|auto> <verlet_rad_delta>`**
  Defines Lennard–Jones potential.
  - `cutoff` : cutoff distance in `SIGMA` units
  - `shift` : energy shift in `EPSILON` units (or `auto` for automatic)
  - `verlet_rad_delta` : Verlet list skin in `SIGMA` units

---

### Integration
- **`INTEGRATION_TIME_STEP <float>`**  :  Timestep in units of `SIGMA * sqrt(MASS/EPSILON)`.  
- **`EQUILIBRATION_STEPS <int>`**  :  Number of equilibration MD steps.  
- **`PRODUCTION_STEPS <int>`**  :  Number of production MD steps.  

---

### Output Control
- **`SAVING_INTERVALS <global_quantities> <configurations> <backup>`**  :  Output frequency for thermodynamic observables, snapshots, and position/velocity backup files.
- **`SAVE_FORMAT <xyz>`**  :  Format for particle trajectory (`xyz` supported).  
- **`SAVE_ATTRIBUTE <name> <eq_flag> <prod_flag>`**  :  Control saving of attributes in equilibration / production runs (flags can be set to `0` or `1`). Options:  
  - `particles`
  - `forces`
  - `velocities`
  - `energy`
  - `virial`
  - … (others listed in example file)

- **`SAVE_POTENTIAL_ENERGY <eq_flag> <prod_flag>`**  
- **`SAVE_PARTIAL_POTENTIAL_ENERGIES <eq_flag> <prod_flag>`**  
- **`SAVE_KINETIC_ENERGY <eq_flag> <prod_flag>`**  
- **`SAVE_TEMPERATURE <eq_flag> <prod_flag>`**  
- **`SAVE_MSD <eq_flag> <prod_flag> <format>`**  
  - `format = all` (per-particle)  
  - `format = mols_cm` (molecule centers of mass)  
- **`SAVE_VIRIAL_PRESSURE <eq_flag> <prod_flag>`**  

- **`SAVING_DIRECTORY_NAME <string>`**  :  Directory name for simulation results.  

---

## Output
The simulation produces:
- **Trajectory files** (particle positions over time)  
- **Thermodynamic logs** (energy, temperature, pressure)  
- **Simulation parameter files**

All outputs are written into the directory specified by  
`SAVING_DIRECTORY_NAME` (default: current working directory).
