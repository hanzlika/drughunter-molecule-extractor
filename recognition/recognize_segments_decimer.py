import os

from rdkit import Chem
from time import time


def recognize_segments(image_list: list) -> dict:
    """
    Iterate through the image list and attempt to recognize each image using DECIMER as a SMILES string.

    Params:
    - image_list (List[Image]): A list of images to recognize.

    Returns:
    - dict: A dictionary containing recognized properties for each image:
            - "smiles": A list of recognized SMILES strings with matching indexing.
            - "inchi": A list of recognized InChI strings with matching indexing.
            - "inchikey": A list of recognized InChIKey strings with matching indexing.
    """
    print("Importing decimer recognition...")
    import_start = time()
    from DECIMER import predict_SMILES
    print(f"Importing took:{time() - import_start} s")

    print("Recognizing with DECIMER V2...")
    recognition_start = time()
    output_dict = {}

    # Initialize lists to store recognized properties
    output_dict["smiles"] = []
    output_dict["inchi"] = []
    output_dict["inchikey"] = []

    for img in image_list:
        # Save the image as a temporary file
        img.save("tmp.png")

        # Use DECIMER to predict the SMILES string for the image
        # predict_SMILES function requires a path to the image as input
        SMILES = predict_SMILES("tmp.png")
        output_dict["smiles"].append(SMILES)

    # Remove the temporary image file
    if image_list:
        os.remove("tmp.png")

    # Convert recognized SMILES to InChI and InChIKey
    for smiles in output_dict["smiles"]:
        mol = Chem.MolFromSmiles(smiles)

        if mol:
            output_dict["inchi"].append(Chem.MolToInchi(mol))
            output_dict["inchikey"].append(Chem.MolToInchiKey(mol))
        else:
            # If the SMILES cannot be converted to a molecule, store empty strings
            output_dict["inchi"].append("")
            output_dict["inchikey"].append("")
    
    print(f"Recognition of {len(image_list)} segments with Decimer took: {time() - recognition_start} s\n({(time() - recognition_start)/len(image_list)} s per segment)")
    return output_dict

def main():
    # Example use:
    input_list = []  # Replace this with the actual list of images you want to recognize
    result_dict = recognize_segments(input_list)

    # Print the recognized properties for each image
    for idx, smiles in enumerate(result_dict["smiles"]):
        inchi = result_dict["inchi"][idx]
        inchikey = result_dict["inchikey"][idx]
        print(f"Image {idx + 1}: SMILES={smiles}, InChI={inchi}, InChIKey={inchikey}")

if __name__ == '__main__':
    main()