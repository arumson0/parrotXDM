#!/usr/bin/env python3
from Input_Output.parse_yaml_input import parse_config
from Input_Output.config import Config
from Data.atom_info import atom
from Data.set_params import set_params
from Input_Output.readers import qe_reader
from Input_Output.readers import aims_reader
from Physics.cpair import cpair
from Physics.ctriple import ctriple
from Physics.epair import epair
from Physics.etriple import etriple
import numpy as np
import sys
import os

def run_driver(yaml_path):
    # Load YAML -> dict -> Config object
    yaml_dict = parse_config(yaml_path)
    config = Config(yaml_dict)

    # Loop over structures
    for path in config.input_path:
        if len(config.input_path) > 1:
            print('\n=====================================================================')
        print(f'XDM({config.damping.upper()}) calculation for {os.path.basename(path)}')
        if config.pairwise and config.triples:
            print("Atom pairs and triples requested.")
        elif config.pairwise and not config.triples:
            print("Atom pairs requested.")
        else:
            print("No energy requested.")

        # Set the damping parameters:
        a1, a2, zdamp = set_params(config.functional.lower(), config.file_type.lower(),config.triples)

        a1_ov = config.override_a1
        a2_ov = config.override_a2

        if (a1_ov is None) != (a2_ov is None):
            raise ValueError("Both 'override_a1' and 'override_a2' must be provided together.")

        if a1_ov is not None and a2_ov is not None:
            a1, a2 = a1_ov, a2_ov / 0.529177

        if config.override_zdamp is not None:
            zdamp = config.override_zdamp

        # Set the integer representation of the selected damping function.
        if config.damping=='bj':
            damp_type_int = 0
            print(f" BJ damping selected. Product of three f_3(Rij) for triples.")
            print(f"    a1        = {a1}")
            print(f"    a2 (Bohr) = {a2}")
        elif config.damping=='z':
            damp_type_int = 1
            print(f" Z damping selected. Product of three sqrt(f_6 (Rij)) for triples.")
            print(f"    z_damp    = {zdamp}")
        
        # Set E_CONVERT for units:
        unit = config.output_unit.lower()
        if unit == 'ha':
            E_CONVERT = 1.0
        elif unit == 'ry':
            E_CONVERT = 2.0
        elif unit == 'ev':
            E_CONVERT = 27.211407953
        elif unit == 'kcal/mol':
            E_CONVERT = 627.509
        elif unit == 'kj/mol':
            E_CONVERT = 2625.5000
        else:
            raise ValueError(f"Unknown output unit: {config.output_unit}")

        # Read XDM info from the output file
        if config.file_type.lower() == 'qe':
            m1, m2, m3, volscl, alpha, symbol, alat, tau, mtrx = qe_reader(path)
        elif config.file_type.lower() == 'aims':
            m1, m2, m3, volscl, alpha, symbol, alat, tau, mtrx = aims_reader(path)
        else:
            raise ValueError(f"Unknown file type: {config.file_type}")

        # Scale polarizabilities
        alpha = np.array([atom(s).pol for s in symbol])
        alpha_scl = alpha * volscl
        # Make tau useful.
        tau = alat * tau

        # Compute and print C6, C8, C10, Rc, Zinv
        l = len(tau)
        c6, c8, c10, rc, zinv, rmax2 = cpair(alpha_scl, m1, m2, m3, symbol, l, path)
        print("+ Cn Coeffcients")
        print("#  i   j  C6          C8          C10         Rc          Zinv")
        for i in range(l):
            for j in range(i+1):
                print(f" {i+1:3} {j+1:3}  {c6[i,j]:.5e} {c8[i,j]:.5e} {c10[i,j]:.5e} {rc[i,j]:.5e} {zinv[i,j]:.5e}")


        # Initalize XDM energies
        e_xdm_pairwise = 0
        e_xdm_triples  = 0
        # Do the pairwise XDM calculation
        if config.pairwise:
            e_xdm_pairwise = epair(
                tau, mtrx, c6, c8, c10, rc, zinv,
                damp_type_int, rmax2,
                a1, a2, zdamp,
                l, E_CONVERT, config.verbose_conv
            )
        print(f"Pairwise XDM({config.damping.upper()}) energy ({config.output_unit}): {E_CONVERT*e_xdm_pairwise:12.7f}" )

        # Do the triple-wise XDM calculation
        if config.triples:
            # First, compute all the C9's
            c9, c11 = ctriple(alpha_scl,m1,m2,m3,l)
            if config.boolc11:
                print("+ Triple-wise Coeffcients")
                print(f"#  i   j   k  C9     C11 (three components)     ")
                for i in range(l):
                    for j in range(i+1):
                        for k in range(j+1):
                            print(f" {i+1:3} {j+1:3} {k+1:3}  {c9[i,j,k]:.5e} {c11[i,j,k,0]:.5e} {c11[i,j,k,1]:.5e} {c11[i,j,k,2]:.5e}")
            else:
                print("+ C9 Coeffcients")
                print(f"#  i   j   k  C9")
                for i in range(l):
                    for j in range(i+1):
                        for k in range(j+1):
                            print(f" {i+1:3} {j+1:3} {k+1:3}  {c9[i,j,k]:.5e}")
            e_xdm_triples = etriple(tau, mtrx,c6,c8,rc,c9,c11,zinv,a1,a2,zdamp,damp_type_int,rmax2,l,E_CONVERT,(not config.extrapolate_triples),config.verbose_conv, config.r2tol,config.boolc11)
            print(f"Triple-wise XDM({config.damping.upper()}) energy ({config.output_unit}): {E_CONVERT*e_xdm_triples:12.7f}" )

        e_xdm_total = e_xdm_pairwise + e_xdm_triples
        print(f"Total XDM({config.damping.upper()}) energy ({config.output_unit}): {E_CONVERT*e_xdm_total: 12.7f}")

# CLI
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python driver.py <input.yaml>")
        sys.exit(1)
    yaml_file = sys.argv[1]
    run_driver(yaml_file)
