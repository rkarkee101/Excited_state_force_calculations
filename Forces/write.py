import numpy as np

file_imp=np.loadtxt("forces_phonons_basis.out",skiprows=1)

convert=np.reshape(file_imp[:,2],(50,3))

exp=np.savetxt("phonon_projection.dat", convert)
