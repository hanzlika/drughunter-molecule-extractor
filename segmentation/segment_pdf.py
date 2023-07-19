from decimer_segmentation import segment_chemical_structures
from pdf2image import convert_from_bytes
from PIL import Image
import numpy as np
import os

def segment_pdf(pdf):
    """
    Iterate through the PDF pages and use decimer_segmentation to get segmented images.

    Params:
    - pdf (bytes): The PDF file content in bytes.

    Returns:
    - List[Image]: A list of segmented images extracted from the PDF.
    """

    # on windows, this only works if poppler is installed and in PATH
    # otherwise the path to poppler needs to be specified in the poppler_path parameter
    segment_list = []
    pages = convert_from_bytes(pdf, 300)
    
    for page_number, page in enumerate(pages):

        # visualization can be set to True for a visual confirmation of the segmentation
        # expand=True yields better results than expand=False
        segments = segment_chemical_structures(np.array(page),
                                               expand=True,
                                               visualization=False)

        for segment in segments:
            image = Image.fromarray(segment)
            segment_list.append(image)

    return segment_list

def main():
    # example use:
    pdf_list = []
    segment_list = []
    for filename, pdf in pdf_list:
        result = segment_pdf(pdf)
        to_append = [(filename, image) for image in result]
        segment_list += to_append

if __name__ == '__main__':
    main()
