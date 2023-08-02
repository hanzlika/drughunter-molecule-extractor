from pdf2image import convert_from_bytes, convert_from_path
from PIL import Image
import numpy as np
import os
import PyPDF2

def segment_pdf(pdfs:list(tuple()), save_to_directory:bool = None, directory:str = None, expand:bool = True, visualization:bool = False):
    from decimer_segmentation import segment_chemical_structures
    """
    Iterate through the PDF pages and use decimer_segmentation to get segmented images.

    Params:
    - pdf (bytes): The PDF file content in bytes.
    - save_to_directory (any value): Whether the segments should also be saved into a directory
    - directory (string): path to the directory where the images should be saved 

    Returns:
    - List[Image]: A list of segmented images extracted from the PDF.
    """

    # on windows, this only works if poppler is installed and in PATH
    # otherwise the path to poppler needs to be specified in the poppler_path parameter
    segment_list = []

    for (filename, content) in pdfs:
        print(f"Attempting to segment {filename}")
        sub_segment_list = []
        pages = convert_from_bytes(content, 300)
        
        for page in pages:

            # visualization can be set to True for a visual confirmation of the segmentation
            # expand=True yields better results than expand=False
            segments = segment_chemical_structures(np.array(page),
                                                expand=expand,
                                                visualization=visualization)

            for segment in segments:
                image = Image.fromarray(segment)
                sub_segment_list.append((filename, image))
        
        print(f"Found {len(sub_segment_list)} segments in {filename}.")

        if save_to_directory:
            if directory:
                os.makedirs(directory, exist_ok=True)
            for index, (filename, segment) in enumerate(sub_segment_list):
                segment.save(os.path.join(directory, f"{filename}_{index}.png"))
        
        segment_list += sub_segment_list

    return segment_list

def main():
    # example use:
    filename = "2020-MOTM-April.pdf"
    filepath = "pdf_extraction/pdfs/2020-MOTM-April.pdf"
    with open(filepath, "rb") as f:
        result = segment_pdf([(filename, f.read())], save_to_directory=True, directory="segmentation/segments/experiments", expand=True)

if __name__ == '__main__':
    main()
