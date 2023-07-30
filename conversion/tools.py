from rdkit import Chem

def smiles_to_inchikey(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Smiles: {smiles} is not a valid smiles string.")
        return "Invalid SMILES"

    inchi = Chem.MolToInchiKey(mol)
    return inchi

def smiles_to_inchi(smiles):
    """
    Convert a SMILES string to an InChI string representation.

    Params:
    - smiles (str): The SMILES string to convert.

    Returns:
    - str: The InChI string representation of the molecule.

    If the input SMILES string is not valid, it returns "Invalid SMILES".
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print(f"Smiles: {smiles} is not a valid smiles string.")
        return "Invalid SMILES"

    inchi = Chem.MolToInchi(mol)
    return inchi

def smiles_list_to_inchi_list(smiles_list):
    """
    Convert a list of SMILES strings to a list of InChI string representations.

    Params:
    - smiles_list (List[str]): A list of SMILES strings to convert.

    Returns:
    - List[str]: A list of InChI string representations of the molecules.

    The function iterates through each SMILES string in the input list and uses the `smiles_to_inchi` function
    to convert it to the corresponding InChI string representation. The resulting InChI strings are appended
    to the `inchi_list` and returned at the end.
    """
    print("Converting smiles to inchi")
    inchi_list = []
    for smiles in smiles_list:
        inchi_list.append(smiles_to_inchi(smiles))
    return inchi_list

def smiles_list_to_inchikey_list(smiles_list):
    print("Converting smiles to inchikeys")
    inchikey_list = []
    for smiles in smiles_list:
        inchikey_list.append(smiles_to_inchikey(smiles))
    return inchikey_list