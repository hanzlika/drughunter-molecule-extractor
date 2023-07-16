print("Importing...")

import argparse
import requests
from bs4 import BeautifulSoup
from calendar import month_name
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

    Params: year and month of the target drug of the month set

    Returns: the first file with the .pdf ending that has been found on the
    web page in the byte format
    
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

def segment_pdf(pdf_list, pdf_month_list):
    """Iterate through the pdfs in pdf_list
    Use decimer_segmentation to get segmented images from the pdfs
    
    Params: list of pdfs and list of months with matching indexes
    
    Returns: list of segments and list of months and segment_numbers with matching indexes
    """

    # on windows, this only works if poppler is installed and in PATH
    # otherwise the path to poppler needs to be specified in the poppler_path parameter
    segment_list = []
    segment_month_list = []
    segment_number_list = []
    for pdf_index, pdf in enumerate(pdf_list):
        segment_count = 0
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
                segment_list.append(image)
                segment_month_list.append(pdf_month_list[pdf_index])
                segment_number_list.append(segment_count)
            print(f"Found {segment_count} segments in the {pdf_month_list[pdf_index]} set page {page_number}.")
    if segment_list:
        return segment_list, segment_month_list, segment_number_list
    else:
        print("No segments featuring chemical structures have been found.")
        return None, None, None

def recognize_segments(image_list):
    """Iterates through the image_list
    and attempts to recognize every image within it using decimer
    as a SMILES string and appends it to the output list

    Params: list of images to recognize

    Returns: list of recognized SMILES strings with matching indexing
    """
    print("Recognizing:")
    smiles_list = []
    for img in image_list:
        # unfortunately predict_smiles only accepts path to image as input
        # so creating a temporary file for this purpose is necessary
        img.save("tmp.png")
        # predict_SMILES has no internal failure mode so try: expect: block is skipped
        # the output list will always be created and contain something, whether
        # that is a valid SMILES or not is checked later
        SMILES = predict_SMILES("tmp.png")
        smiles_list.append(SMILES)
    os.remove("tmp.png")
    return smiles_list

def smiles_to_inchi(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Smiles: {smiles} is not a valid smiles string.")
        return "Invalid SMILES"
    inchi = Chem.MolToInchi(mol)
    return inchi

def smiles_list_to_inchi_list(smiles_list):
    print("Converting smiles to inchi")
    inchi_list = []
    for smiles in smiles_list:
        inchi_list.append(smiles_to_inchi(smiles))
    return inchi_list
    
def validate_inchi(inchi_list):
    print("Validating inchi through UNICHEM.")
    # update in case the api changes again
    url = "https://www.ebi.ac.uk/unichem/api/v1/compounds"

    payload = {"type": "inchi", "compound": ""}

    headers = {"Content-Type": "application/json"}

    validated_list = []
    for inchi in inchi_list:
        payload["compound"] = inchi
        response = requests.request("POST", url, json=payload, headers=headers)
        resp_dict = response.json()
        if resp_dict["response"] != "Not found":
            validated_list.append(True)
        else:
            validated_list.append(False)

    return validated_list

def main(target_year):
    ### PDF EXTRACTION ##########################################################################
    # attempts to download all drugs of the month sets for each month within the year
    pdf_list = []
    month_list = []
    for index in range(1, 13):
        print(f"Attempting to download {month_name[index]}-{target_year}")
        result = download_pdf(target_year, month_name[index])
        if result:
            pdf_list.append(result)
            month_list.append(month_name[index])
    
    # no sets were successfully extracted
    if not pdf_list:
        print("Failed to download any sets.")
        return
    
    ### SEGMENTATION ############################################################################

    segment_list, segment_month_list, segment_number_list = segment_pdf(pdf_list, month_list)
    
    if not segment_list:
        return
  
    ### RECOGNITION #############################################################################

    smiles_list = recognize_segments(segment_list)

    ### VALIDATION ##############################################################################

    inchi_list = smiles_list_to_inchi_list(smiles_list)

    validated_list = validate_inchi(inchi_list)

    ### EXPORT ##################################################################################
    
    print("Exporting results.")    

    results = {"year":  [target_year] * len(validated_list), "segment_month": segment_month_list, "segment_number": segment_number_list, \
               "smiles": smiles_list, "inchi": inchi_list, "validated by Unichem": validated_list}
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