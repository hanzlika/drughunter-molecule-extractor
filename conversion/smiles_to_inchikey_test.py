from tools import smiles_to_inchikey
import os

input_path = "recognition/recognized_smiles.txt"
output_path = "conversion/converted_inchikey.txt"
smiles_list = []

with open(input_path, 'r') as file:
    for line in file:
        smiles_list.append(line)

inchi_key_list = []
for smiles in smiles_list:
    inchi_key_list.append(smiles_to_inchikey(smiles))

with open(output_path, 'w') as file:
    for inchi_key in inchi_key_list:
        file.write(inchi_key + "\n")
        
