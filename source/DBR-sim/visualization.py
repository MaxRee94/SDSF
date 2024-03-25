from multiprocessing import Value
import cv2
import numpy as np
import matplotlib.pyplot as plt


def visualize_image(img, image_width):
    img_resized = cv2.resize(img, (image_width, image_width),
               interpolation = cv2.INTER_LINEAR)

    cv2.imshow("DBR Simulation (TEST)", img_resized)
    cv2.waitKey(1)
    
def get_image(grid):
    img = grid.get_distribution()
    img *= 255
    img = img.astype(np.uint8)
    
    return img

def visualize(grid, image_width=1000):
    img = get_image(grid)    
    visualize_image(img, image_width)
    
    return img

def visualize_difference(image1, image2, image_width=1000):
    img = np.array([image2, image2, image1], np.uint8)
    img = img.transpose(1, 2, 0)
    
    visualize_image(img, image_width)



    
    
