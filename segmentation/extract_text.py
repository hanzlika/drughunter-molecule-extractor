import fitz
import os

def extract_text_molecules_of_the_month(
        page,
        bboxes : list[tuple[int, int, int, int]],
        tolerance : int = 10
        ) -> list[str]:
    
    textboxes = []
    page_x0, page_y0, page_x1, page_y1 = page.rect

    for bbox in bboxes:
        bbox_x0, bbox_y0, bbox_x1, bbox_y1 = bbox[1] * page_x1, bbox[0] * page_y1, bbox[3] * page_x1, bbox[2] * page_y1
        if bbox_x1 > page_x1 / 2:
            target_text_box = ( 
                bbox_x1 - tolerance, 
                bbox_y0 - tolerance,
                page_x1, 
                bbox_y1 + tolerance )
        else:
            target_text_box = bbox_x1 - tolerance, bbox_y0 - tolerance, \
            page_x1 / 2, bbox_y1 + tolerance
        textboxes.append(page.get_textbox(target_text_box))

    return textboxes

def extract_text_in_direction(
        page,
        bboxes : list[tuple[int, int, int, int]],
        tolerance : int = 15
        ) -> list[str]:
    
    textboxes = []
    page_x0, page_y0, page_x1, page_y1 = page.rect

    for bbox in bboxes:
        bbox_x0, bbox_y0, bbox_x1, bbox_y1 = bbox[1] * page_x1, bbox[0] * page_y1, bbox[3] * page_x1, bbox[2] * page_y1
        target_text_box =  ( 
            bbox_x0 - tolerance, 
            bbox_y1, 
            bbox_x1 + tolerance, 
            2 * bbox_y1 - bbox_y0 )
        textboxes.append(page.get_textbox(target_text_box))
        print(page.rect)
        print(bbox)
        print(bbox_x0, bbox_y0, bbox_x1, bbox_y1)
        print(target_text_box)

    return textboxes

def extract_text(
        path : str, 
        page_num : int, 
        bboxes : list[tuple],
        direction : str = 'right'
        ) -> list[str]:
    
    doc = fitz.open(path)
    page = doc[page_num]
    textboxes = []

    if direction == 'right':
        textboxes = extract_text_molecules_of_the_month(page, bboxes, 10)
    if direction == 'down':
        textboxes = extract_text_in_direction(page, bboxes, 15)

    return textboxes