from multiprocessing import Value
import cv2
import numpy as np
import matplotlib.pyplot as plt

from helpers import *


def get_color_dict(no_values, begin=0.0, end=1.0):
    x_step = (end * no_values - begin * no_values) / no_values
    x_range = [ no_values * begin + x_step * x for x in range(no_values) ]
    color_dict = {v : 
        np.array((
            get_max((255 - (x * (2 * 255/no_values))), 0),
            get_max(255 - get_max(255 - (x * (2 * 255/no_values)), 0) - get_max((-765/3 + (x * (2 * 255/no_values))), 0), 0),
            get_max((-765/3 + (x * (2 * 255/no_values))), 0)
        ), np.uint8) for v, x in zip(range(1, no_values + 1), x_range)
    }
    color_dict[0] = np.array((0, 0, 0), np.uint8)
    color_dict[-5] = np.array((0, 80, 220))
    color_dict[-6] = np.array((0, 0, 110))
    #for key, val in color_dict.items():
    #   print(key, " -- ", val)
    return color_dict

def visualize_image(img, image_width):
    img_resized = cv2.resize(img, (image_width, image_width),
               interpolation = cv2.INTER_LINEAR)

    cv2.imshow("DBR Simulation (TEST)", img_resized)
    cv2.waitKey(1)

def get_image(grid, collect_states, color_dict):
    img = grid.get_distribution(collect_states)
    if (color_dict == {0:0, 1:255}):
        img *= 255
    elif (len(color_dict[0]) == 3):
        new_img = np.zeros((grid.size, grid.size, 3), np.uint8)
        for i in range(grid.size):
            for j in range(grid.size):
                new_img[i, j] = color_dict[img[i, j]]
        img = new_img.copy()
            
    img = img.astype(np.uint8)
    
    return img

def visualize(grid, image_width=1000, collect_states=True, color_dict={0:0, 1:255}):
    img = get_image(grid, collect_states, color_dict)    
    visualize_image(img, image_width)
    
    return img

def save_image(img, path, image_width = 1000):
    img_resized = cv2.resize(img, (image_width, image_width),
               interpolation = cv2.INTER_LINEAR)
    cv2.imwrite(path, img_resized)

def visualize_difference(image1, image2, image_width=1000):
    img = np.array([image2, image2, image1], np.uint8)
    img = img.transpose(1, 2, 0)
    
    visualize_image(img, image_width)



    
    
