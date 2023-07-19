from DECIMER import predict_SMILES
import os

def recognize_segments(image_list):
    """
    Iterate through the image list and attempt to recognize each image using DECIMER as a SMILES string.

    Params:
    - image_list (List[Image]): A list of images to recognize.

    Returns:
    - List[str]: A list of recognized SMILES strings with matching indexing.
    """
    print("Recognizing:")
    smiles_list = []
    for img in image_list:
        # Save the image as a temporary file
        img.save("tmp.png")

        # Use DECIMER to predict the SMILES string for the image
        # predict_SMILES function requires a path to the image as input
        SMILES = predict_SMILES("tmp.png")
        smiles_list.append(SMILES)

    # Remove the temporary image file
    if image_list:
        os.remove("tmp.png")

    return smiles_list

def main():
    # example use
    input_list = []
    smiles_list = recognize_segments(image_list)

    for smiles in smiles_list:
        print(smiles)

if __name__ == '__main__':
    main()