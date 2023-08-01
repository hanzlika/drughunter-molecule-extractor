from multiprocessing import freeze_support
import torch
from molscribe import MolScribe
from huggingface_hub import hf_hub_download
from rdkit import Chem
from PIL import Image
import os 
import numpy as np


def recognize_segments(image_list):
    """
    Recognize each image form a list using Molscribe

    In: list of images to recognize

    Out: dict with the structure:
        {
            'smiles' = [smiles_1, smiles_2, ...]
            'inchi' = [inchi_1, inchi_2, ...]
            'inchikey' = [inchikey_1, inchikey_2, ...]
        }
    """

    # Necessary
    freeze_support()

    # download the model
    ckpt_path = hf_hub_download('yujieq/MolScribe', 'swin_base_char_aux_1m.pth')
    model = MolScribe(ckpt_path, device=torch.device('cpu'))

    results_list = model.predict_images(image_list)

    output_dict = {}
    
    output_dict['smiles'] = []
    output_dict['inchi'] = []
    output_dict['inchikey'] = []

    for result in results_list:
        inchi = Chem.MolBlockToInchi(result['molfile'])

        output_dict['smiles'].append(result['smiles'])
        output_dict['inchi'].append(inchi)
        output_dict['inchikey'].append(Chem.InchiToInchiKey(inchi))

    return output_dict

def main():
    path = 'segmentation/segments'

    files = os.listdir(path)

    image_list = []
    for file in files:
        image = np.asarray(Image.open(os.path.join(path, file)))
        image_list.append(image)

    recognize_segments(image_list)

if __name__ == '__main__':
    main()



