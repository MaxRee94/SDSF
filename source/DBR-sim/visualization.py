import cv2
import numpy as np


def visualize(grid):
    print("visualizing..")
    tilesize = 10
    boardheight = 100
    boardwidth = 100
    screenheight = boardheight * tilesize
    screenwidth = boardwidth * tilesize
    board_center = (int(boardwidth / 2), int(boardheight / 2))
    screen_center = (int(screenwidth / 2), int(screenheight / 2))

    img = np.zeros((screenheight, screenwidth, 1), np.uint8) + 255
    board = np.zeros((boardwidth, boardheight, 1), np.uint8)
