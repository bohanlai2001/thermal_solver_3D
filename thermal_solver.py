import argparse
import numpy as np
import os 
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve


def parse_spice_file(filename):
    """
    Parse a .sp file containing thermal resistances, heat flux sources,
    and fixed temperature sources.
    Returns the elements, node list, mapping, and system size.
    """
    thermal_resistances = []
    fixed_temp_sources = []   # Equivalent to voltage sources (fixed T)
    heat_flux_sources = []    # Equivalent to current sources (heat input/output)
    node_set = set()

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if not line or line.startswith("."):
                continue

            tokens = line.split()
            if len(tokens) < 4:
                continue

            element = tokens[0]
            n1, n2 = tokens[1], tokens[2]
            value = float(tokens[3])
            node_set.update([n1, n2])

            if element[0] == 'R':  # Thermal resistance
                thermal_resistances.append((n1, n2, value))
            elif element[0] == 'V':  # Fixed temperature source
                fixed_temp_sources.append((element, n1, n2, value))
            elif element[0] == 'I':  # Heat flux source
                heat_flux_sources.append((n1, n2, value))

    def sort_key(name):
        digits = ''.join(filter(str.isdigit, name))
        return int(digits) if digits else float('inf')

    node_list = sorted([n for n in node_set if n != '0'], key=sort_key)
    node_to_idx = {name: idx for idx, name in enumerate(node_list)}
    num_nodes = len(node_list)
    num_fixedT = len(fixed_temp_sources)
    size = num_nodes + num_fixedT

    return thermal_resistances, fixed_temp_sources, heat_flux_sources, node_list, node_to_idx, size


def build_matrices(thermal_resistances, fixed_temp_sources, heat_flux_sources, node_list, node_to_idx, size):
    """
    Build the thermal conductance matrix (K) and heat flow vector (Q)
    using nodal analysis for the given thermal network.
    """
    row = []
    col = []
    data = []
    Q = np.zeros(size)   # Heat flow vector

    def idx(node):
        return node_to_idx[node]

    # Thermal resistances (like electrical conductance)
    for n1, n2, rth in thermal_resistances:
        g = 1 / rth  # thermal conductance
        if n1 != '0':
            i = idx(n1)
            row.append(i); col.append(i); data.append(g)
        if n2 != '0':
            j = idx(n2)
            row.append(j); col.append(j); data.append(g)
        if n1 != '0' and n2 != '0':
            i, j = idx(n1), idx(n2)
            row.append(i); col.append(j); data.append(-g)
            row.append(j); col.append(i); data.append(-g)

    # Heat flux sources (like current sources)
    for n1, n2, q in heat_flux_sources:
        if n1 != '0':
            Q[idx(n1)] -= q
        if n2 != '0':
            Q[idx(n2)] += q

    # Fixed temperature sources (like voltage sources)
    for k, (name, n1, n2, tval) in enumerate(fixed_temp_sources):
        row_idx = len(node_list) + k
        if n1 != '0':
            i = idx(n1)
            row += [row_idx, i]
            col += [i, row_idx]
            data += [1.0, 1.0]
        if n2 != '0':
            j = idx(n2)
            row += [row_idx, j]
            col += [j, row_idx]
            data += [-1.0, -1.0]
        Q[row_idx] = tval

    # Thermal conductance matrix
    K = csr_matrix((data, (row, col)), shape=(size, size))
    return K, Q


def write_temperature_output(filename, node_list, T, fixed_temp_sources):
    """
    Write the computed node temperatures and branch values
    of fixed temperature sources to an output text file.
    """
    with open(filename, 'w') as f:
        f.write("Node\tTemperature [K]\n----\t-------------\n")
        for i, name in enumerate(node_list):
            f.write(f"{name}\t{T[i]:.6f}\n")

        f.write("Fixed Temp Branches\n------\t-------------\n")
        offset = len(node_list)
        for k, (tname, _, _, _) in enumerate(fixed_temp_sources):
            branch_val = T[offset + k]
            f.write(f"{tname}#branch\t{branch_val:.6f}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Thermal Solver")
    parser.add_argument("--input_file", required=True, help="Path to main thermal SPICE file")
    parser.add_argument("--additional_power", default=None, help="Optional path to additional power SPICE file")
    parser.add_argument("--output_file", default=None, help="Path to output temperature file")
    args = parser.parse_args()
#    input_file = "benchmarks/thermal_grid_1.sp"
#    output_file = input_file
#    output_file = output_file.replace(".sp", ".temperature")
#    output_file = output_file.replace("benchmarks/", "results/")
#    additional_file = "additional_power/thermal_grid_3D_additional_power.sp"

    rth, fixedT, qsrc, node_list, node_to_idx, size = parse_spice_file(args.input_file)

    if args.additional_power is not None:
        add_rth, add_fixedT, add_qsrc, add_node_list, add_node_to_idx, add_size = parse_spice_file(args.additional_power)
        fixedT.extend(add_fixedT)
        qsrc.extend(add_qsrc)
        size= size + add_size - len(add_node_list)

    K, Q = build_matrices(rth, fixedT, qsrc, node_list, node_to_idx, size)
    T = spsolve(K, Q)  # Solve steady-state thermal equation
    os.makedirs("results", exist_ok=True)
    write_temperature_output(args.output_file, node_list, T, fixedT)
