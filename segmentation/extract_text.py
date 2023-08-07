import fitz
import os

def scale_bbox(bbox, max_x, max_y):
    return (bbox[1] / 2480) * max_x, (bbox[0] / 3508) * max_y,  (bbox[3] / 2480) * max_x, (bbox[2] / 3508) * max_y

def extract_text(path, bbox):
    doc = fitz.open(path)  # any supported document type
    page = doc[0]  # we want text from this page


    rect = page.rect
    tolerance = 15

    converted_bbox = scale_bbox(bbox, max_x=rect[2], max_y=rect[3])
    if converted_bbox[2] > rect[2] / 2:
        target_text_box = converted_bbox[2] - tolerance, converted_bbox[1] - tolerance, \
        rect[2], converted_bbox[3] + tolerance
    else:
        target_text_box = converted_bbox[2] - tolerance, converted_bbox[1] - tolerance, \
        page.rect[2] / 2, converted_bbox[3] + tolerance
    return page.get_textbox(target_text_box)
