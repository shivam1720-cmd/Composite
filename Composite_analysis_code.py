import math
import numpy as np
import sympy as sym

# Constants

E1 = float(input("Enter E1 in GPa: ")) #38.6
E2 = float(input("Enter E2 in GPa: "))#8.27
V12 = float(input("Enter V12: "))#0.28
V21 = (E2 * V12) / E1
G12 = float(input("Enter G12 in GPa: "))#4.14
t = float(input("Enter THICKNESS OF EACH LAMINA in mm: "))#0.125  # mm



K = (1 - V12 * V21)
Q11 =  E1/K
Q22 =  E2/K
Q12 = V12*E2/K
Q66 = G12
Q21 = Q12
Q16 = Q61 = 0
Q26 = Q62 = 0
sigma1_T_u = float(input("stress_1_Tensile_u in MPa "))#1062
sigma1_T_c = float(input("stress_1_compressive_u in MPa "))#-610
sigma2_T_u = float(input("stress_2_Tensile_u in MPa "))#31
sigma2_T_c = float(input("stress_2_compressive_u in MPa "))#118
shear12 = float(input("shear_stress in MPa "))#72
a1 = float(input("alpha1 in microm/m "))* pow(10, -6)#8.6 * pow(10, -6)
a2 = float(input("alpha1 in microm/m "))* pow(10, -6)#22.1 * pow(10, -6)
del_T =  float(input("change in Temperature ")) #50
n= int(input("Enter total number of lamina excluding symmetric: "))
total_t = float(n)*t
#N = sym.Matrix([[100], [0], [0]])


N_m=100
Maxm=[]

#j_values = [0, math.pi / 2, math.pi / 2, 0]
print("type all value of Lamina in degree value in sequence")
j_values= [0]*2*n
for i in range (2*n):
    j_values[i]=float(input()) * math.pi / 180
    

store_layer =n+1
for loop in range(n):
    N= sym.Matrix([[N_m], [0], [0]])  
    Q = sym.Matrix([[Q11, Q12, Q16],
                    [Q21, Q22, Q26],
                    [Q61, Q62, Q66]])
    R = sym.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 2]])
    

    
    R = sym.Matrix([[1, 0, 0], [0, 1, 0], [0, 0, 2]])

    def T_a(theta):
        c = math.cos(theta)
        s = math.sin(theta)
        T = sym.Matrix([[c**2, s**2, 2 * s * c], [s**2, c**2, -2 * s * c], [-1 * s * c, s * c, c**2 - s**2]])
        return T

    def Q_bar(theta):
        c = math.cos(theta)
        s = math.sin(theta)
        T = sym.Matrix([[c**2, s**2, 2 * s * c], [s**2, c**2, -2 * s * c], [-1 * s * c, s * c, c**2 - s**2]])
        Q_ba = T.inv() * Q * R * T * R.inv()
        return Q_ba

    Q_bar_0 = Q_bar(0)
    A = sym.zeros(3, 3)
    B = sym.zeros(3, 3)
    D = sym.zeros(3, 3)
    C = sym.zeros(6, 6)

    z = sym.zeros(2 * n + 2, 1)
    for i in range(1, (2 * n + 2)):
        z[i] = -(total_t) / 2 + (i - 1) * t
    k = 0    

    for i in range(1, 2 * n + 1):
        k = (j_values[i - 1])
        if k != (store_layer):
            A += Q_bar(k) * (z[i + 1] - z[i])
            B += Q_bar(k) * (z[i + 1] ** 2 - z[i] ** 2) * .5
            D += Q_bar(k) * (z[i + 1] ** 3 - z[i] ** 3) * (1 / 3)
    #sym.pretty_print(A)

    C = sym.Matrix([[A, B], [B, D]])

    A_inv = A.inv()
    Ex = 1 / (1 * A_inv[0, 0])

    
    ex = A.inv() * N * .001

    e = [0] * (n + 1)
    sigma = [0] * (n + 1)
    for i in range(1, n + 1):
        
        k = j_values[i - 1]
        e[i] = R * T_a(k) * R.inv() * ex
    e_t = sym.zeros(3, n)
    for i in range(0, n):
        e_t[0, i] = e[i + 1]

    for i in range(1, n + 1):
         k = (j_values[i - 1])
         #if k != (store_layer):
         sigma[i] = Q * e[i] * 1000
    sigma_ta = sym.zeros(3, n)

    for j in range(0, n):
        sigma_ta[0, j] = sigma[j + 1]
    sigma_t = [0] * (n + 1)
    sigma_t = sigma_ta.T
    #sym.pretty_print(sigma_t)

    sr = sym.zeros(n, 3)

    for i in range(0, 3):
        for j in range(0, n):
            
            if i == 0:
                store = sigma_t[j, i]
                if store > 0:
                    sr[j, i] = store / sigma1_T_u
                else:
                    sr[j, i] = store / sigma1_T_c
            if i == 1:
                store = sigma_t[j, i]
                if store > 0:
                    sr[j, i] = store / sigma2_T_u
                else:
                    sr[j, i] = store / sigma2_T_c
            if i == 2:
                store = sigma_t[j, i]
                sr[j, i] = store / shear12
     
    
    if loop >0:
        for element in Maxm:
            for i in range(0, 3):
                sr[element, i] = 0
          
            
          
    absolute_max = np.amax(np.abs(sr))
    sr_np = np.array(sr).astype(float)
    max_index = np.argmax(np.abs(sr_np))
    max_row, max_col = np.unravel_index(max_index, sr_np.shape)
    Maxm.append(max_row)
    #print("stress ratio",absolute_max)
    #sym.pretty_print(sr)
    #print(Maxm)
   
    

    Nx = N_m / absolute_max

    alpha1 = sym.Matrix([[a1], [a2], [0]])
    alpha = [0] * (2 * n)
    for i in range(1, n + 1):
        k = j_values[i - 1]
        alpha[i - 1] = R * T_a(k).inv() * alpha1
        alpha[2 * n - i] = alpha[i - 1]

    Nxx = sym.zeros(3, 1)
    Mxx = sym.zeros(3, 1)
    for i in range(1, 2 * n + 1):
        k = j_values[i - 1]
        #if k != (store_layer):
        Nxx += del_T * Q_bar(k) * alpha[i - 1] * (z[i + 1] - z[i]) * pow(10, 6)
        Mxx += (del_T / 2) * Q_bar(k) * alpha[i - 1] * (z[i + 1] ** 2 - z[i] ** 2)

    NM = sym.Matrix([[Nxx], [Mxx]])
    e_x = C.inv() * NM * pow(10, -6)

    exx = np.zeros([2 * n, 3, 1])
    e_xy = np.zeros([3, 1])
    k_xy = np.zeros([3, 1])
    for i in range(3):
        e_xy[i] = e_x[i]
        k_xy[i] = e_x[i + 3]

    exx = np.zeros([2 * n, 3, 1])
    for i in range(1, 2 * n + 1):
        exx[i - 1] = e_xy + z[i] * k_xy

    e_tx = [0] * (2 * n)

    for i in range(1, 2 * n + 1):
        e_tx[i - 1] = del_T * alpha[i - 1]

    e_rx = exx - e_tx

    sigma_rs = [0] * (2 * n)
    for i in range(1, 2 * n + 1):
        k = j_values[i - 1]
        #if k != (store_layer):
        sigma_rs[i - 1] = Q_bar(k) * e_rx[i - 1] * 1000

    sigma_rs_p = [0] * (2 * n)
    for i in range(1, 2 * n + 1):
        k = j_values[i - 1]
        
        sigma_rs_p[i - 1] = T_a(k) * sigma_rs[i - 1]

    sigma_t = np.array(sigma_t).astype(float)
    

    sigma_rt = sigma_rs_p[max_row][max_col]

    sigma_nx = sigma_t[max_row][max_col]
    sigma_total = sigma_t[max_row][max_col] + sigma_rs_p[max_row][max_col]
    

    store = sigma_total
    sigma_residual_stress = sigma_rs_p[max_row][max_col]
    if max_col == 0:
        print("LT failure",'at angle ',j_values[max_row]*180/math.pi,'degree')
        if store > 0:
            sigma_n = sigma1_T_u - sigma_residual_stress
        else:
            sigma_n = sigma1_T_c - sigma_residual_stress
    elif max_col == 1:
        print("TT failure",'at angle',j_values[max_row]*180/math.pi,'degree')
        if store > 0:
            sigma_n = sigma2_T_u - sigma_residual_stress
        else:
            sigma_n = sigma2_T_c - sigma_residual_stress
    elif max_col == 2:
        print("stress failure",'at angle',j_values[max_row]*180/math.pi,'degree')
        sigma_n = shear12 - sigma_residual_stress

    FPF_load = N_m / (sigma_nx / sigma_n)
    #print(max_row)
    
    store_layer = j_values[max_row]
    #print(store_layer)    
    N_m=FPF_load
    #print(N_m)
    
    
    print("PF Load:", FPF_load,"KN/m","\n")
    

print("\nLPF Load:",FPF_load ,"KN/m")


# Test the function
