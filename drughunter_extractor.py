from argparse import ArgumentParser
from calendar import month_name

from export.export_results import export_to_csv
from pdf_extraction.pdf_extraction import download_pdf
from recognition.recognize_segments_decimer import \
    recognize_segments as recognize_segments_decimer
from recognition.recognize_segments_molscribe import \
    recognize_segments as recognize_segments_molscribe
from segmentation.segment_pdf import segment_pdf
from validation.validate_with_chembl_webresource import validate_inchikey_list


def extract_molecules_from_pdfs(pdfs: list[tuple[str, bytes]]) -> None:
    """
    Extracts molecules from a list of PDFs by performing the following steps:
    1. Downloads PDF files from the specified URLs using the `download_pdf` function.
    2. Segments the PDFs into images using the `segment_pdf` function.
    3. Recognizes SMILES strings from the segmented images using the `recognize_segments` function.
    4. Converts recognized SMILES strings to InChI strings using the `smiles_list_to_inchi_list` function.
    5. Validates the InChI strings using the `validate_inchi` function.

    Parameters:
    - pdfs (List[Tuple[str, bytes]]): A list of tuples containing the source file name and the content of the PDF.
      Each tuple represents one PDF to be processed.

    The function does not return any value, but it saves the results in CSV format to a file named "with_combined_recognition_YYYY_MM_DD-HH_MM_SS.csv",
    where "YYYY_MM_DD-HH_MM_SS" is the current date and time.

    The function also prints the statistics of successful and unsuccessful recognitions using Molscribe and Decimer.
    """

    extraction_results = {}


    # segmentation_results -> list of tuples (source filename, segment)
    segmentation_results = segment_pdf(pdfs, save_to_directory=True, directory="segmentation/segments/2023") 

    segments = [segment for (_, segment) in segmentation_results]
    extraction_results['source'] = [source for source, _ in segmentation_results]
    
    recognition_results_molscribe = recognize_segments_molscribe(segments)
    recognition_results_molscribe['validation'] = validate_inchikey_list(recognition_results_molscribe.get('inchikey'))

    index_list = []

    for index, validation_result in enumerate(recognition_results_molscribe['validation']):
        if validation_result is False:
            index_list.append(index)

    print(f"Molscribe succesfully recognized {recognition_results_molscribe['validation'].count(True)} segments. \
          ({round(recognition_results_molscribe['validation'].count(True)/len(recognition_results_molscribe['validation'])*100, 2)} % success rate)")
    print(f"Molscribe could not recognize {recognition_results_molscribe['validation'].count(False)} segments. Attemping to recognize them with decimer")

    recognition_results_decimer = recognize_segments_decimer([segment for index, segment in enumerate(segments) if index in index_list ])
    recognition_results_decimer['validation'] = validate_inchikey_list(recognition_results_decimer.get('inchikey'))


    # complement molscribe unsuccessful recognitions with successful decimer recognitions
    for index, validation_result in enumerate(recognition_results_decimer['validation']):
        if validation_result is True:
            for key in ['smiles', 'inchi', 'inchikey', 'validation']:
                recognition_results_molscribe[key][index_list[index]] = recognition_results_decimer[key][index]
    print(f"Decimer managed to recognize {recognition_results_decimer['validation'].count(True)} more segments.")
    print(f"{recognition_results_molscribe['validation'].count(True)} segments recognized in total. \
          ({round(recognition_results_molscribe['validation'].count(True)/len(recognition_results_molscribe['validation'])*100, 2)} % success rate)")
    extraction_results.update(recognition_results_molscribe)

    if extraction_results:
        export_to_csv(extraction_results, file_name="with_combined_recognition_")


def extract_molecules_of_the_month(target_year: int) -> None:
    """
    Extract molecules of the month for a specific target year.

    Params:
    - target_year (int): The target year to extract molecules of the month.

    This function iterates over the months (from index 1 to 12) and constructs the URL for each month
    using the `target_year`. Using download_pdf, it constructs a list of pdfs that the molecules need
    to be extracted from.
    It then calls the `extract_molecules_from_pdfs` function to extract molecules
    from the pdf list. The extracted segment, SMILES, InChI, and validation results are appended to their respective lists.

    Finally, the function calls the `export_results` function to export the results to a CSV file.
    """

    pdfs = []
    for index in range(1, 2):
        url = f"https://drughunter.com/molecules-of-the-month/{target_year}/{month_name[index].lower()}-{target_year}"
        pdfs += download_pdf(url, download_all=True)

    extract_molecules_from_pdfs(pdfs)


def main():
    """
    Main function to handle command-line arguments and execute the appropriate extraction based on the arguments.

    The function uses the `argparse` module to parse the command-line arguments.
    - If the `url` argument is provided, it calls the `extract_molecules` function with the provided URL.
    - If the `set` argument is not provided, it prints an error message.
    - If the `year` argument is not provided or is in an invalid format, it prints an error message.
    - If the `set` is "month", it calls the `extract_molecules_of_the_month` function with the provided year.
    - If the `set` is "drug", it calls the `extract_approved_drugs` function with the provided year.
    - If the `set` is not "month" or "drug", it prints an error message.

    Note: The code for the `extract_approved_drugs` function is not provided in the current context.
    """
    parser = ArgumentParser(description='DrugHunter extractor')
    parser.add_argument('-y', '--year', type=int, help='targeted year of drughunter sets')
    parser.add_argument('-s', '--set', type=str, help='targeted set (currently supports "month" for drugs of the month and \"drug\" for approval reviews from that year)')
    parser.add_argument('-u', '--url', type=str, help='url of webpage with targeted set (in case the format of drughunter url changes, which is likely)')
    args = parser.parse_args()
    
    # user requested a specific url to be checked, year and set is ignored
    if args.url:
        extract_molecules_from_pdfs(download_pdf(args.url, download_all=False))
        return

    # user has not picked the set
    if not args.set:
        print("No set or url provided")
        return
    
    # user has not picked the year
    if not args.year:
        print("No year provided")
        return

    if len(str(args.year)) != 4:
        print("Invalid year format. Please provide a 4-digit year (YYYY).")
        return
    
    if args.set == "month":
        if (args.year < 2020):
            print("DrugHunter molecules of the month start at the year 2020, please provide a year equal or greater.")
            return
        extract_molecules_of_the_month(args.year)
        return

    if args.set == "drug":
        if (args.year < 2019):
            print("DrugHunter Drug approvals start at the year 2019, please provide a year equal or greater.")
        extract_approved_drugs(args.year)
        return

    print("The set parameter currently supports \"month\" and \"drug\" options, please use one of them, or provide a url to the web page you want pdf extracted from")

if __name__ == "__main__":
    main()