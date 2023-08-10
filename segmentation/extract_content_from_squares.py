import cv2
import numpy as np
from PIL import Image

def unique_core(core_set, new_core):
    for core in core_set:
        if abs(core[0] - new_core[0]) < 10 and abs(core[1] - new_core[1]) < 10:
            return False
    return True

def sort_segments_bboxes(
        segments : list[np.array], 
        bboxes : list[tuple[int, int, int, int]], #(y0, x0, y1, x1)
        same_row_pctg_threshold = 0.01
        ) -> tuple[np.array, list[tuple[int, int, int, int]]]:
    """
    Sorts segments and bounding boxes in "reading order"

    Args:
        segments - image segments to be sorted
        bboxes - bounding boxes containing edge coordinates of the image segments
        same_row_pixel_threshold - how many pixels apart can two pixels be to be considered "on the same row"
    
    Returns:
        segments and bboxes in reading order
    """

    # Sort by y-coordinate (top-to-bottom reading order)
    sorted_bboxes = sorted(bboxes, key=lambda bounding_box: bounding_box[0])

    # Group bounding boxes by rows based on y-coordinate
    rows = []
    current_row = [sorted_bboxes[0]]
    for bounding_box in sorted_bboxes[1:]:
        if abs(bounding_box[0] - current_row[-1][0]) < same_row_pctg_threshold:  # You can adjust this threshold as needed
            current_row.append(bounding_box)
        else:
            rows.append(sorted(current_row, key=lambda x: x[1]))  # Sort by x-coordinate within each row
            current_row = [bounding_box]
    rows.append(sorted(current_row, key=lambda x: x[1]))  # Sort the last row

    # Flatten the list of rows and return
    sorted_bboxes = [bounding_box for row in rows for bounding_box in row]

    sorted_segments = [segments[bboxes.index(bounding_box)] for bounding_box in sorted_bboxes]
    return sorted_segments, sorted_bboxes

def extract_content_with_borders(page):
    images, bounding_boxes = [], []

    image = np.array(page)  # Convert the PIL image to a NumPy array
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    _, thresh = cv2.threshold(gray,50,255,0)
    contours, _ = cv2.findContours(thresh, 1, 2)
    
    core_set = set()
    for contour in contours:
        approx = cv2.approxPolyDP(contour, 0.01*cv2.arcLength(contour, True), True)

        
        
        if len(approx) == 4:  # Check for quadrilateral shape
            x, y, w, h = cv2.boundingRect(contour)
            aspect_ratio = w / float(h)

            if 0.9 <= aspect_ratio <= 1.1 and w > 200 and w < 500:  # Check for approximately square shape
                if unique_core(core_set, (x, y)):
                    content = image[y:y+h, x:x+w]
                    core_set.add((x, y))
                    images.append(content)
                    bounding_box = (y / float(image.shape[0]), x / float(image.shape[1]), (y + h) / float(image.shape[0]), (x + h)/image.shape[1])
                    bounding_boxes.append((bounding_box))

    sorted_images, sorted_bboxes = sort_segments_bboxes(images, bounding_boxes, 0.01)

    return  sorted_images, sorted_bboxes
