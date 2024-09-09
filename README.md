# Excited_state_force_calculations using BerkeleyGW and Quantum ESPRESSO
Use the attached files as your sample input reference.

###########################################
# Workflow for Running Calculations in the `WFN_fi` Directory
###########################################

1. Perform Fine Grid Calculation:
   pw.x < wfn.in > wfn.out

2. Conduct Phonon Calculation:
   ph.x < ph.in > ph.out

3. Convert Wavefunctions for BerkeleyGW:
   pw2bgw.x < pw2bgw.in > pw2bgw.out

4. (Optional) Perform Dynamic Matrix Calculation:
   dynmat.x < dynmat.in > dynmat.out


###########################################
# Important Notes:
###########################################

- Number of Bands (nbnd):
  -> Ensure it includes the total number of valence and conduction bands used in the BSE calculation.

- Smearing:
  -> Apply smearing, even for semiconductors (suggested value: 1e-3 Ry), 
     to allow `ph.x` to calculate electron-phonon (el-ph) coefficients.

- Critical Variables to Check:
  -> nogg=.true.
  -> electron_phonon='simple'
  -> fildvscf='NAME'
  -> trans=.true.

- Convergence:
  -> Typically, tr2_ph=1d-16 or 1d-18 is sufficient.

- Symmetry:
  -> Set search_sym=.false. to ensure that Quantum ESPRESSO computes the el-ph coefficients for all k-points.

- XML File Validation:
  -> Ensure that elph.xx.yy.xml files end with </Root>. If not, delete the last line to correct it.
     You can run the following script to fix these files:
     bash fix_elph_xml_files.bash


###########################################
# Example Input for `ph.x`:
###########################################

phonon_calc
$inputph
    verbosity = 'high'
    prefix = 'LiF'
    outdir = './'
    fildyn = 'dyn'
    fildvscf = 'dvscf'
    electron_phonon='simple'
    trans=.true.
    nogg=.true.
    tr2_ph=1.0d-18
    search_sym=.false.
/
0.0 0.0 0.0


###########################################
# BerkeleyGW Workflow
###########################################

- Epsilon, Sigma, Kernel:
  -> Follow the usual steps.

- Absorption:
  -> Use the momentum operator. The polarization direction does not matter since the BSE solution is required.
  -> WFNq_fi is not needed when using the momentum operator.


###########################################
# Excited State Forces Calculation
###########################################

To run the calculation:
python3 /path_to_script/excited_forces.py

Input File (forces.inp) Variables:
  - iexc = 1
  - jexc = 1  # Generally, iexc = jexc. The code calculates <iexc | dH^bse | jexc>
  - Calculate_Kernel = False  # Set to False to avoid using the kernel.
  - calc_IBL_way = False
  - eqp_file = eqp.dat  # QP energies on the fine grid (produced by the absorption calculation). Point to eqp.dat in the absorption directory.
  - exciton_file = eigenvectors.h5
  - acoutic_sum_rule = True
  - el_ph_dir = ../../5-wfn_fi_6kpts/_ph0/MAPI.phsave/  # Directory where the el-ph coefficients are stored.
  - just_real = True

Output:
  - The calculated forces will be saved in forces_cart.out, containing the excited state forces for each atom in Cartesian coordinates, in units of eV/Å.

Important Note:
  - The total force equals the sum of the excited state force and the DFT force.


###########################################
# Expressing Forces in Phonon Displacement Basis
###########################################

1. Generate Displacement Patterns:
   dynmat.x < dynmat.in > dynmat.out

Example `dynmat.in`:
$input
fildyn='dyn'
asr='crystal'
fileig='eigvecs'
filxsf='displacements.axsf'
/

- The output file displacements.axsf will contain the displacement patterns for each phonon.

2. Convert Forces to Phonon Basis:
   python cart2ph_eigvec.py forces_cart.out displacements.axsf

- The output file forces_phonons_basis.out will contain the forces expressed in the phonon displacement basis. 
  Here, forces_cart.out is the output of the excited state forces calculation.

