from multiprocessing import Value
import cv2
import numpy as np
import matplotlib.pyplot as plt


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



    
    
