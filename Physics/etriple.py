import numpy as np
import time
import sys
import os
os.environ["OMP_NUM_THREADS"] = "32"
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../Kernels')))
import c9_kernel

def etriple(tau, mtrx,c6,c8,rc,c9,c11,zinv,a1,a2,zdamp,damp_type_int,rmax2,l,E_CONVERT,act_conv,verbose,r2tol,boolc11):
    e9_total = 0
    energies = []
    delta = []
    max_shell = 99 # Like that would ever work lol
    threshold = 3.675E-8 # case for act_conv
    
    if verbose:
       print("Starting the triple-wise XDM computation...")
       if not act_conv:
           print(f"Extrapolation requested with R^2 = {r2tol}")

    start_time = time.time()
    # F90: c9_loops(tau,mtrx,c6,rc,c9,zinv,a1,a2,zdamp,damp,rmax2,n,l,e9_shell)
    for n in range(0, max_shell+1):
        e_shell = c9_kernel.c9_loops(
            tau, mtrx, c6, c8, rc, c9, c11, zinv, 
            a1, a2, zdamp, 
            damp_type_int, rmax2,
            n=n, l=l, boolc11=boolc11
        )
    
        e9_total += (1/6) * e_shell # 1/3 prevents triple counting
        energies.append(e9_total)
    
        if len(energies) > 1:
            delta = np.append(delta, energies[-1] - energies[-2])
            if len(delta)<=7 or act_conv:
                if verbose:
                    print(f"Shell {n}:     E_total = {energies[-1]:5.9f}     ΔE = {delta[-1]:5.9f}")

            if act_conv:
                if abs(delta[-1]) < threshold:
                    break
            else:
            # Extrapolate the convergence early
                if abs(delta[-1]) < threshold:
                    print("Actual convergence achieved before log-linearity. Using the actual value...")
                    break
                y = np.log(abs(delta))
                if len(y) > 7:
                    yy = y[len(y)-7:]
                    nn = np.arange(len(y)-7,len(y))
                    # Perform linear regression
                    X = np.vstack([np.ones(7),nn]).T
                    B = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(yy)
                    y_hat = X.dot(B)
                    # Quality of the regression:
                    SS_res = np.sum((yy - y_hat)**2)
                    SS_tot = np.sum((yy - np.mean(yy))**2)
                    R2 = 1 - SS_res / SS_tot
                    if verbose:
                        print(f"Shell {n}:     E_total = {energies[-1]:5.9f}     ΔE = {delta[-1]:5.9f}     R^2 = {R2:5.3f}")
                    if R2 > r2tol:
                        if verbose:
                            print("LOG-LINEAR CONVERGENCE ACHIEVED! Extrapolating...")
                        r = np.exp(B[1])
                        E_conv = energies[-1] + delta[-1]*r/(1-r)
                        e9_total = E_conv
                        break
    end_time = time.time()
    if verbose:
        print(f"Computed in: {end_time - start_time:.2f} seconds")

    return e9_total
