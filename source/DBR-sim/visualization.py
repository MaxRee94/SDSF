from multiprocessing import Value
import cv2
import numpy as np
import matplotlib.pyplot as plt


def visualize(grid, image_width=1000):
    print("visualizing..")
    tilesize = 10
    boardheight = grid.size
    boardwidth = grid.size
    screenheight = boardheight * tilesize
    screenwidth = boardwidth * tilesize
    board_center = (int(boardwidth / 2), int(boardheight / 2))
    screen_center = (int(screenwidth / 2), int(screenheight / 2))

    img = grid.get_distribution()
    img *= 255
    img = img.astype(np.uint8)
    
    img_resized = cv2.resize(img, (image_width, image_width),
               interpolation = cv2.INTER_LINEAR)

    cv2.imshow("DBR Simulation (TEST)", img_resized)
    cv2.waitKey(100000)
    cv2.destroyAllWindows()
    
    
