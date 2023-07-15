from rdkit import Chem
import requests
import pandas as pd


def smiles_to_inchi(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Smiles: {smiles} is not a valid smiles string.")
        return ""
    inchi = Chem.MolToInchi(mol)
    return inchi

def read_smiles_from_file(file_path):
    smiles_list = []
    with open(file_path, 'r') as file:
        for line in file:
            smiles_list.append(line.strip())
    return smiles_list

def smiles_list_to_inchi_list(smiles_list):
    inchi_list = []
    for smiles in smiles_list:
        inchi_list.append(smiles_to_inchi(smiles))
    return inchi_list

file_path = "recognition/recognized_smiles.txt"
smiles_list = read_smiles_from_file(file_path)
inchi_list = smiles_list_to_inchi_list(smiles_list)

url = "https://www.ebi.ac.uk/unichem/api/v1/compounds"

payload = {
    "type": "inchi",
    "compound": ""
}
headers = {"Content-Type": "application/json"}

validated_inchi_list = []
for index, inchi in enumerate(inchi_list):
    payload["compound"] = inchi
    response = requests.request("POST", url, json=payload, headers=headers)
    resp_dict = response.json()
    if resp_dict["response"] != "Not found":
        validated_inchi_list.append(inchi)
    else:
        validated_inchi_list.append("Not found")

results = {"smiles": smiles_list, "inchi": validated_inchi_list}
df = pd.DataFrame(results)
df.to_csv('results.csv')

    
