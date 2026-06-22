import numpy as np

def ctriple(alpha_scl, m1, m2, m3, l):
    # Create C9 array:
    c9   = np.zeros((l,l,l))
    c11  = np.zeros((l,l,l,3))

    # Populate the array
    for i in range(l):
        for j in range(l):
            for k in range(l):
                Hi = m1[i] / alpha_scl[i]
                Hj = m1[j] / alpha_scl[j]
                Hk = m1[k] / alpha_scl[k]
                excit = (Hi + Hj + Hk) / ((Hi + Hj)*(Hi + Hk)*(Hj + Hk))
                c9[i][j][k]  = m1[i]*m1[j]*m1[k] * excit

                c11[i][j][k][0] = (9/15)*m2[i]*m1[j]*m1[k] * excit
                c11[i][j][k][1] = (9/15)*m1[i]*m2[j]*m1[k] * excit
                c11[i][j][k][2] = (9/15)*m1[i]*m1[j]*m2[k] * excit
                
    return c9, c11
