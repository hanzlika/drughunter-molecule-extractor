from decimer_segmentation import segment_chemical_structures
from pdf2image import convert_from_bytes
from PIL import Image
import numpy as np
import os

def segment_pdf(pdf, save_to_directory=None, directory=None, expand=True, visualization=False):
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
    pages = convert_from_bytes(pdf, 300)
    
    for page in pages:

        # visualization can be set to True for a visual confirmation of the segmentation
        # expand=True yields better results than expand=False
        segments = segment_chemical_structures(np.array(page),
                                               expand=expand,
                                               visualization=visualization)

        for segment in segments:
            image = Image.fromarray(segment)
            segment_list.append(image)

    if save_to_directory:
        if directory:
            os.makedirs(directory, exist_ok=True)
        for index, segment in enumerate(segment_list):
            segment.save(os.path.join(directory, f"{index}.png"))

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
