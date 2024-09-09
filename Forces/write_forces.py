import numpy as np

force_file='./forces_cart.out'   ####### Results are in eV/A
dtype = [('Atom', int), ('dir', 'S1'), ('RPA_diag', float), ('RPA_diag_offiag', float)]
forces=np.loadtxt(force_file, dtype=dtype)

RPA_off_diag=np.zeros(len(forces))

for j in range(len(forces)):
    RPA_off_diag[j]=forces[j][3]
        
force_arange=np.reshape(RPA_off_diag,(50,3))


exp=np.savetxt("force_arange.dat", force_arange)
