from chembl_webresource_client.unichem import unichem_client as unichem
from time import sleep

def validate_inchikey_list(inchikey_list : list[str]) -> list[bool]:
    print("Validating through unichem.connectivity()...")
    """
    Validate a list of InChIKeys using the UniChem web service.

    Parameters:
        inchikey_list (list[str]): A list of InChIKeys to validate.

    Returns:
        list[bool]: A list of validation results where True means the InChIKey is valid, and False means it's not valid.
    """
    validation_list = []
    for inchikey in inchikey_list:
        try:
            unichem.connectivity(inchikey)
            validation_list.append(True)
        except:
            validation_list.append(False)
        sleep(0.1)
    return validation_list

def main():
    """
    Main function to validate a list of InChIKeys from a file.

    The list of InChIKeys should be stored in a text file with one InChIKey per line.
    The function reads the file, extracts the InChIKeys, and then validates them using the UniChem web service.
    Finally, it prints the validation results.
    """
    input_path = "conversion/converted_inchikey.txt"
    inchikey_list = []
    with open(input_path, 'r') as file:
        for line in file:
            inchikey_list.append(line.strip())

    validation_results = validate_inchikey_list(inchikey_list)

    for idx, inchikey in enumerate(inchikey_list):
        valid_str = "Valid" if validation_results[idx] else "Invalid"
        print(f"InChIKey {idx + 1}: {inchikey} - {valid_str}")


if __name__ == '__main__':
    main()
