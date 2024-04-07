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
    color_dict[-5] = np.array((0, 150, 220))
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
        new_img = np.zeros((grid.width, grid.width, 3), np.uint8)
        for i in range(grid.width):
            for j in range(grid.width):
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


class Graphs:
    def __init__(self, dynamics, width=6, height=4):
        self.dynamics = dynamics
        self.types = {"histogram": ["Firefree Interval"], "timeseries": ["Tree cover"]}
        self.graphs = []
        for graph_type, datatypes in self.types.items():
            for datatype in datatypes:
                if graph_type == "timeseries":
                    _graph = Graph(dynamics, datatype, width, height)
                elif graph_type == "histogram":
                    _graph = Histogram(dynamics, datatype, width, height)
                self.graphs.append(_graph)
    def update(self):
        for graph in self.graphs:
            graph.update()


class Graph:
    def __init__(self, dynamics, datatype, width=6, height=4):
        self.times = [0]
        self.dynamics = dynamics
        self.datatype = datatype
        self.data = []
        self.add_datapoint()

        # to run GUI event loop
        plt.ion()
 
        # here we are creating sub plots
        self.figure, self.ax = plt.subplots(figsize=(width, height))
        self.line1, = self.ax.plot(self.times, self.data)
 
        # setting title
        plt.title(datatype, fontsize=12)
 
        # setting x-axis label and y-axis label
        plt.xlabel("Timestep (years)")
        plt.ylabel(datatype)
      
    def add_datapoint(self):
        print("datatype: ", self.datatype)
        if self.datatype == "Tree cover":
            self.data.append(self.dynamics.state.grid.get_tree_cover())
    def update(self):
        self.add_datapoint()
        self.times.append(self.dynamics.time)
 
        # updating data values
        self.line1.set_xdata(self.times)
        self.line1.set_ydata(self.data)
        
        plt.xlim(self.times[0], self.times[-1])
        plt.ylim(min(self.data), max(self.data))
 
        # drawing updated values
        self.figure.canvas.draw()
 
        # This will run the GUI event
        # loop until all UI events
        # currently waiting have been processed
        self.figure.canvas.flush_events()


class Histogram:
    def __init__(self, dynamics, datatype, width=6, height=4, no_bins = 40):
        self.no_bins = no_bins
        self.dynamics = dynamics
        self.datatype = datatype
        self.data = None
        self.replace_data()

        # to run GUI event loop
        plt.ion()
 
        # here we are creating sub plots
        self.figure, self.ax = plt.subplots(figsize=(width, height))
        N, bins, self.bar_container = self.ax.hist(self.data, bins=self.no_bins)
        self.rect_width = self.bar_container.patches[0].get_width()
        initial_xrange = bins[-1] - bins[0]
        initial_hist_width = 0
        for rect in self.bar_container.patches:
            initial_hist_width += rect.get_width()
        self.width_over_range = initial_xrange / initial_hist_width
 
        # setting title
        plt.title(datatype, fontsize=12)
 
        # setting x-axis label and y-axis label
        plt.xlabel(datatype)
        plt.ylabel("Number of cells")

    def replace_data(self):
        if self.datatype == "Firefree Interval":
            self.data = self.dynamics.get_firefree_intervals()

    def update(self):
        self.replace_data()
 
        counts, bins = np.histogram(self.data)

        # Update x limits
        self.ax.set_xlim([bins[0], bins[-1]])
        self.ax.set_ylim([min(counts), max(counts)])
        
        # updating data values
        xrange = bins[-1] - bins[0]
        x = 0
        rect_width = (xrange * self.width_over_range) / self.no_bins
        for count, rect in zip(counts, self.bar_container.patches):
            rect.set_height(count)
            rect.set_width(rect_width)
            rect.set_x(x)
            x += rect_width
 
        # drawing updated values
        self.figure.canvas.draw()
 
        # This will run the GUI event
        # loop until all UI events
        # currently waiting have been processed
        self.figure.canvas.flush_events()
     



    
    
