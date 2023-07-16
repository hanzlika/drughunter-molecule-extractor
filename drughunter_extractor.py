print("Importing...")

import argparse
import requests
from bs4 import BeautifulSoup
import calendar
from decimer_segmentation import segment_chemical_structures, segment_chemical_structures_from_file
from DECIMER import predict_SMILES
from pdf2image import convert_from_bytes
from PIL import Image
import numpy as np
import os
from rdkit import Chem
import pandas as pd


def download_pdf(target_year: int, target_month: str):
    """Attempt to download a pdf file from drughunter molecules
    of the month based on given target year and month
    
    As of 13/7/2023 as seen in 13-7-2023-test.log
    succesfully downloaded all available molecules of the month pdfs
    (published from february 2020 to may 2023)
    """

    # general target format
    url = f"https://drughunter.com/molecules-of-the-month/{target_year}/{target_month}-{target_year}"

    # headers simulating a browser request are necessary, without them the request status code returns 403 (forbidden)
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36"
    }

    # headers=headers prevents 403 forbidden response.status_code
    response = requests.get(url, headers=headers)

    if response.status_code == 200:
        # parses the web page for the first pdf file (assuming the molecules
        # of the month file is the first one on the web page)
        soup = BeautifulSoup(response.content, 'html.parser')
        pdf_link = soup.select_one('a[href$=".pdf"]')

        if pdf_link:
            file_name = pdf_link['href'].split('/')[-1]
            file_url = pdf_link['href']
            
            pdf_response = requests.get(file_url, headers=headers)

            if pdf_response:
                print(file_name + " downloaded successfully.")
                return pdf_response.content
        else:
            print("PDF link not found on the page.")
    else:
        print("Failed to download the PDF.")
    return None

def segment_pdf(pdf_dict):
    """
    Uses decimer_segmentation to take a pdf page
    from input_path
    and segment it into individual images
    containing chemical structures.

    Segments are collected into a dict and are returned
    with a tuple key in the format (month_number, segment_number) 
    """

    # on windows, this only works if poppler is installed and in PATH
    # otherwise the path to poppler needs to be specified in the poppler_path parameter
    segment_dict = {}
    for key in pdf_dict.keys():
        segment_count = 0
        pdf = pdf_dict[key]
        pages = convert_from_bytes(pdf, 300)

        for page_number, page in enumerate(pages):

            # visualization can be set to True for a visual confirmation of the segmentation
            # expand=True yields better results than expand=False
            segments = segment_chemical_structures(np.array(page),
                                                   expand=True,
                                                   visualization=False)

            for segment in segments:
                segment_count += 1
                image = Image.fromarray(segment)
                segment_dict[(key, segment_count)] = image
            print(f"Found {segment_count} segments in the {calendar.month_name[key]} set.")
    if segment_dict:
        return segment_dict
    else:
        print("No segments featuring chemical structures have been found.")

def recognize_segments(image_dict):
    """Iterates through the input_directory
    and attempts to recognize every file within it
    as a SMILES string and appends it to the output list
    """
    smiles_dict = {}
    for key in image_dict.keys():
        print(f"Recognizing image {key[1]} from the {calendar.month_name[key[0]]} set.")
        tmp_img = image_dict[key]
        tmp_img.save("tmp.png")
        SMILES = predict_SMILES("tmp.png")
        smiles_dict[key] = SMILES
    os.remove("tmp.png")
    return smiles_dict

def smiles_to_inchi(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Smiles: {smiles} is not a valid smiles string.")
        return ""
    inchi = Chem.MolToInchi(mol)
    return inchi

def smiles_dict_to_inchi_dict(smiles_dict):
    inchi_dict = {}
    for key in smiles_dict.keys():
        inchi_dict[key] = smiles_to_inchi(smiles_dict[key])
    return inchi_dict
    
def main(target_year):
    ### PDF EXTRACTION ##########################################################################
    # attempts to download all drugs of the month sets for each month within the year
    pdf_dict = {}
    for index in range(1, 13):
        target_month = calendar.month_name[index]
        print(f"Attempting to download {target_month}-{target_year}")
        result = download_pdf(target_year, target_month)
        if result:
            pdf_dict[index] = result
    
    # no sets were successfully extracted
    if not pdf_dict:
        print("Failed to download any sets.")
        return
    
    ### SEGMENTATION ############################################################################

    segment_dict = segment_pdf(pdf_dict)
    
    if not segment_dict:
        return

    ### RECOGNITION #############################################################################

    smiles_dict = recognize_segments(segment_dict)

    ### VALIDATION AND EXPORT ###################################################################
    
    print("Validating recognized smiles through UNICHEM.")
        
    inchi_dict = smiles_dict_to_inchi_dict(smiles_dict)

    # update in case the api changes again
    url = "https://www.ebi.ac.uk/unichem/api/v1/compounds"

    payload = {
        "type": "inchi",
        "compound": ""
    }

    headers = {"Content-Type": "application/json"}

    validated_dict = {}
    for key in inchi_dict:
        payload["compound"] = inchi_dict[key]
        response = requests.request("POST", url, json=payload, headers=headers)
        resp_dict = response.json()
        if resp_dict["response"] != "Not found":
            validated_dict[key] = True
        else:
            validated_dict[key] = False
    
    print("Exporting results.")    
    year_list, month_list, n_list = [], [], []
    for key in validated_dict.keys():
        year_list.append(target_year)
        month_list.append(key[0])
        n_list.append(key[1])

    results = {"year":  year_list, "month": month_list, "n": n_list, \
               "smiles": smiles_dict.values(), "inchi": inchi_dict.values(), "validated by Unichem": validated_dict.values()}
    df = pd.DataFrame(results)
    df.to_csv('results.csv')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='DrugHunter extractor')
    parser.add_argument('-y', '--year', type=int, help='targeted year of drughunter sets')
    args = parser.parse_args()

    target_year = args.year
    if target_year:
        if len(str(target_year)) != 4:
            print("Invalid year format. Please provide a 4-digit year (YYYY).")
        elif (target_year < 2020):
            print("DrugHunter molecules of the month start at the year 2020, please provide a year equal or greater.")
        else:
            main(target_year)
    else:
        print("No year provided.")