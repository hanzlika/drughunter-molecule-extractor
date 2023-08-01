from chembl_webresource_client.unichem import unichem_client as unichem

def validate_inchikey_list(inchikey_list):
    validation_list = []
    for inchikey in inchikey_list:
        try:
            ret = unichem.connectivity(inchikey)
            validation_list.append('True')
        except:
            validation_list.append('False')
    return validation_list

def main ():

    input_path = "conversion/converted_inchikey.txt"
    inchikey_list = []
    with open(input_path, 'r') as file:
        for line in file:
            inchikey_list.append(line.strip())


    print(validate_inchikey_list(inchikey_list))


if __name__ == '__main__':
    main()