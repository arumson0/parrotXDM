import numpy as np
import time
import sys
import os
os.environ["OMP_NUM_THREADS"] = "32"
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Kernels')))
import cpair_kernel


def epair(tau, mtrx, c6, c8, c10, rc, zinv, damp_type_int,rmax2,a1,a2,zdamp,l,E_CONVERT,verbose):
    # Convergence tools:
    n_vecs = np.array([0,0,0], dtype=int)
    not_converged = True
    e_vdw_change = 1
    vdw_energy_threshold = 3.675E-8
    vdw_energies_list = []
    if verbose:
        print("Starting the pairwise XDM computation...")
    start_time = time.time()

    e6_total   = 0.0
    e8_total   = 0.0
    e10_total  = 0.0
    exdm_total = 0.0
    while not_converged:
        # Call Fortran routine for the current supercell
        exdm_shell, e6, e8, e10 = cpair_kernel.pairwise_xdm(
            tau, mtrx, rc, c6, c8, c10, zinv,
            a1, a2, zdamp,
            damp_type_int, rmax2,
            n_vecs,
            l
        )
        e6_total   += 0.5 * e6
        e8_total   += 0.5 * e8
        e10_total  += 0.5 * e10
        exdm_total += 0.5 * exdm_shell # 0.5 accounts for double counting.
        vdw_energies_list.append(exdm_total)

        # Check convergence
        if len(vdw_energies_list) > 1:
            e_vdw_change = vdw_energies_list[-1] - vdw_energies_list[-2]

        if verbose:
            print(f"Shell {n_vecs[0]}:     E_total = {vdw_energies_list[-1]:5.9f}     ΔE = {e_vdw_change:5.9f}")

        if abs(e_vdw_change) <= vdw_energy_threshold:
            not_converged = False
            if verbose:
                print("Supercell convergence achieved!")

        # Increase supercell size for next iteration
        n_vecs += 1

    end_time = time.time()
    if verbose:
        print('Final supercell size:', n_vecs)
        print(f"Computed in: {end_time - start_time:.2f} seconds")
    e_xdm_pairwise = exdm_total
    return e_xdm_pairwise, e6_total, e8_total, e10_total
