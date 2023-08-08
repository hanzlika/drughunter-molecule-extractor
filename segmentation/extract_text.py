import fitz
import os

def scale_bbox(bbox, max_x, max_y):
    return (bbox[1] / 2480) * max_x, (bbox[0] / 3508) * max_y,  (bbox[3] / 2480) * max_x, (bbox[2] / 3508) * max_y

def extract_text(path : str,  page_num : int, bboxes : list[tuple]) -> list[str]:
    doc = fitz.open(path)  
    
    page = doc[page_num]  

    rect = page.rect
    tolerance = 15

    textboxes = []

    for bbox in bboxes:
        converted_bbox = scale_bbox(bbox, max_x=rect[2], max_y=rect[3])
        if converted_bbox[2] > rect[2] / 2:
            target_text_box = converted_bbox[2] - tolerance, converted_bbox[1] - tolerance, \
            rect[2], converted_bbox[3] + tolerance
        else:
            target_text_box = converted_bbox[2] - tolerance, converted_bbox[1] - tolerance, \
            page.rect[2] / 2, converted_bbox[3] + tolerance
        textboxes.append(page.get_textbox(target_text_box))

    return textboxes
