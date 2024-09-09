#### USE python forces_merge.py forces_cart.out displacement.axsf 

import numpy as np
import os
import sys

forces_cart_file = sys.argv[1]
eigvecs_file = sys.argv[2]

source_file = eigvecs_file
force_file=forces_cart_file

destination_file = './test_file.axsf'
destination_file2 = './test_dummy.axsf'

dtype = [('Atom', int), ('dir', 'S1'), ('RPA_diag', float), ('RPA_diag_offiag', float)]
forces=np.loadtxt(force_file, dtype=dtype)


with open(source_file, 'r') as source:
    
    while True:
        line=source.readline()
        if line=='':
            break
        if line.find('PRIMCOORD') !=-1:
            line=source.readline()
            linesplit=line.split()
            num_lines=int(linesplit[0])
            break

            
            
with open(source_file, 'r') as source:        
    with open(destination_file, 'w') as destination:        
        lines = source.readlines()[:num_lines+8]
        destination.writelines(lines[8:num_lines])
        
with open(source_file, 'r') as source:        
    with open(destination_file2, 'w') as destination2:        
        lines = source.readlines()[:num_lines+8]
        destination2.writelines(lines[1:8])
        
def add_line_to_file(file_path, line_to_add):
    try:
        # Create a temporary file
        temp_file_path = file_path + '.tmp'

        # Write the line to add and the original file's contents to the temporary file
        with open(temp_file_path, 'w') as temp_file:
            temp_file.write(line_to_add + '\n')
            with open(file_path, 'r') as original_file:
                temp_file.write(original_file.read())

        # Replace the original file with the temporary file
        os.replace(temp_file_path, file_path)
        
    except Exception as e:
        print(f"Error occurred while adding the line to the file: {e}")

# Usage example
file_path = './test_dummy.axsf'
line_to_add = 'ANIMSTEPS 1'

add_line_to_file(file_path, line_to_add)
        

RPA_off_diag=np.zeros(len(forces))

for j in range(len(forces)):
    RPA_off_diag[j]=forces[j][3]
        
force_arange=np.reshape(RPA_off_diag,(num_lines,3))
        
np.savetxt('dummy_force_arange.txt',force_arange,fmt='%.5f')
        

# Usage example
file1_path = 'test_file.axsf'
file2_path = 'dummy_force_arange.txt'
output_file_path = 'dummy_final.axsf'


def replace_columns(file1, file2, output_file):
    with open(file1, 'r') as f1, open(file2, 'r') as f2, open(output_file, 'w') as out:
        for line1, line2 in zip(f1, f2):
            # Split the lines into columns
            columns1 = line1.split()
            columns2 = line2.split()

            # Replace the last three columns of file1 with columns from file2
            columns1[-3:] = columns2

            # Join the columns with whitespace and write the line to the output file
            replaced_line = ' '.join(columns1) + '\n'
            out.write(replaced_line)




replace_columns(file1_path, file2_path, output_file_path)


def merge_files(file1, file2, merged_file):
    with open(file1, 'r') as f1, open(file2, 'r') as f2, open(merged_file, 'w') as merged:
        # Read the contents of file1 and write them to the merged file
        merged.write(f1.read())

        # Read the contents of file2 and write them to the merged file
        merged.write(f2.read())

    print(f"Successfully merged {file1} and {file2} into {merged_file}.")


# Usage example
file1_path = './test_dummy.axsf'
file2_path = './dummy_final.axsf'
merged_file_path = './final_displacement.axsf'

merge_files(file1_path, file2_path, merged_file_path)

os.remove('test_dummy.axsf')
os.remove('dummy_final.axsf')
os.remove('dummy_force_arange.txt')
os.remove('test_file.axsf')
