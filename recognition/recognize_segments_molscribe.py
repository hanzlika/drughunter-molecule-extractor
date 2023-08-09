import os
from multiprocessing import freeze_support

import numpy as np
from PIL import Image
from rdkit import Chem
from time import time
import torch
from huggingface_hub import hf_hub_download
from chembl_structure_pipeline import standardizer
from molscribe import MolScribe

def recognize_segments(image_list: list) -> dict:
    """
    Recognize each image from a list using MolScribe.

    Params:
        image_list (List[Image]): A list of images to recognize.

    Returns:
        dict: A dictionary with the structure:
            {
                'smiles': [smiles_1, smiles_2, ...],
                'inchi': [inchi_1, inchi_2, ...],
                'inchikey': [inchikey_1, inchikey_2, ...]
            }
    """

    print("Recognizing with MolScribe...")
    recognition_start = time()

    # Necessary for Windows compatibility with multiprocessing
    freeze_support()

    # Download the model checkpoint from Hugging Face model hub
    ckpt_path = hf_hub_download('yujieq/MolScribe', 'swin_base_char_aux_1m.pth')
    model = MolScribe(ckpt_path, device=torch.device('cpu'))

    np_arr_list = [np.asarray(image) for image in image_list]

    results_list = model.predict_images(np_arr_list)

    output_dict = {}

    output_dict['smiles'] = []
    output_dict['inchi'] = []
    output_dict['inchikey'] = []

    for result in results_list:

        output_dict['smiles'].append(result['smiles'])

        mol = Chem.MolFromSmiles(result['smiles'])

        if mol:
            output_dict['inchi'].append(Chem.MolToInchi(mol))
            output_dict['inchikey'].append(Chem.MolToInchiKey(mol))
        else:
            output_dict['inchi'].append("")
            output_dict['inchikey'].append("")

    print(f"Recognition of {len(image_list)} segments with MolScribe took: {time() - recognition_start} s\n({(time() - recognition_start)/len(image_list)} s per segment)")
    return output_dict


def main():
    # Example use: Replace 'path' with the directory containing the segmented images
    path = "segmentation/segments/with_expand"
    files = os.listdir(path)

    image_list = []
    for file in files:
        image = np.asarray(Image.open(os.path.join(path, file)))
        image_list.append(image)

    result_dict = recognize_segments(image_list)

    # Print the recognized properties for each image
    for idx, smiles in enumerate(result_dict['smiles']):
        inchi = result_dict['inchi'][idx]
        inchikey = result_dict['inchikey'][idx]
        print(f"Image {idx + 1}: SMILES={smiles}, InChI={inchi}, InChIKey={inchikey}")


if __name__ == '__main__':
    main()
