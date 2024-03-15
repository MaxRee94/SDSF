from multiprocessing import Value
import cv2
import numpy as np
import matplotlib.pyplot as plt


def visualize(grid):
    print("visualizing..")
    tilesize = 10
    boardheight = grid.size
    boardwidth = grid.size
    screenheight = boardheight * tilesize
    screenwidth = boardwidth * tilesize
    board_center = (int(boardwidth / 2), int(boardheight / 2))
    screen_center = (int(screenwidth / 2), int(screenheight / 2))

    img = grid.get_distribution() * 255
    img = img.astype(np.uint8)
    cv2.imshow("Simulation (TEST)", img)
    cv2.waitKey(100000)
    cv2.destroyAllWindows()
    
    
