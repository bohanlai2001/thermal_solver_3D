# thermal_solver_3D
python thermal_solver.py --input_file benchmarks/thermal_grid_3D.sp --additional_power additional_power/thermal_grid_3D_additional_power.sp --output results/thermal_grid_3D.temperature

A Python-based thermal solver that extends 2D thermal grid analysis into 3D cube-shaped node structures. The solver constructs the thermal conductance matrix (G matrix) using Modified Nodal Analysis (MNA) and computes node temperatures under steady-state conditions.

The solver supports:

3D cube node structures (with vertical vias mapped into node-to-node connections).

External power maps as additional input.

SPICE-like netlist format (.sp files) for thermal resistances, fixed temperature sources, and heat flux sources.
