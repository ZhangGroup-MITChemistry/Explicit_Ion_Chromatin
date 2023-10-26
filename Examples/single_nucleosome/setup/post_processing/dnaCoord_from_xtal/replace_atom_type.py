import numpy as np
import sys

# replace the atom type index

new_atom_type_file_path = sys.argv[1] # path for the file that contains the new atom type index
template_file_path = sys.argv[2] # path for the template file
output_file_path = sys.argv[3] # path fot the output file

# we will replace the atom type index in the template file with atom type index given by the file that contains the new atom type index 

# read the new atom type index
new_atom_type_file = open(new_atom_type_file_path,'r')
new_atom_type_lines = new_atom_type_file.readlines()
n = len(new_atom_type_lines)
for i in range(n):
    if new_atom_type_lines[i][:5] == 'Atoms':
        i_start = i + 2
    if new_atom_type_lines[i][:5] == 'Bonds':
        i_end = i - 1
        
new_atom_type_list = []
for j in range(i_start,i_end):
    atom_line = new_atom_type_lines[j]
    atom_list = atom_line.split()
    new_atom_type_list.append(int(atom_list[2]))

new_atom_type_file.close()

# read template file
template_file = open(template_file_path,'r')
template_file_lines = template_file.readlines()
template_file.close()

# replace old atom type index with new atom type index
n = len(template_file_lines)
for i in range(n):
    if template_file_lines[i][:5] == 'Atoms':
        i_start = i + 2
    if template_file_lines[i][:5] == 'Bonds':
        i_end = i - 1

k = 0
for j in range(i_start,i_end):
    atom_line = template_file_lines[j]
    atom_list = atom_line.split()
    atom_list[3] = str(new_atom_type_list[k])
    k += 1
    template_file_lines[j] = '\t'.join(atom_list)
    template_file_lines[j] = '\t' + template_file_lines[j] + '\n'

# write the output
output_file = open(output_file_path,'w')
for each_line in template_file_lines:
    output_file.write(each_line)
output_file.close()



