from decimer_segmentation import segment_chemical_structures, segment_chemical_structures_from_file
from pdf2image import convert_from_path
from PIL import Image
import numpy as np
import os

def segment_pdf(input_path, output_directory):
    """
    Uses decimer_segmentation to take a pdf page
    from input_path
    and segment it into individual images
    containing chemical structures.

    Afterwards saves the segments into designated output_directory
    """

    os.makedirs(output_directory, exist_ok=True)

    # on windows, this only works if poppler is installed and in PATH
    # otherwise the path to poppler needs to be specified in the poppler_path parameter
    pages = convert_from_path(input_path, 300)

    for page_number, page in enumerate(pages):

        # visualization can be set to True for a visual confirmation of the segmentation
        # expand=True yields better results than expand=False
        segments = segment_chemical_structures(np.array(page),
                                               expand=True,
                                               visualization=False)

        for segment_number, segment in enumerate(segments):
            image = Image.fromarray(segment)
            image.save(f"{output_directory}\page_{page_number}_segment_{segment_number}.png")

def main():
    input_path = 'pdf_extraction\pdfs\DH-MOTM-Poster-May2023.pdf'
    output_directory = 'segmentation\segments'
    segment_pdf(input_path, output_directory)

if __name__ == '__main__':
    main()
