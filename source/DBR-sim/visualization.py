from multiprocessing import Value
import cv2
import numpy as np
import matplotlib.pyplot as plt

from helpers import *
from config import *


def get_color_dict(no_values, begin=0.0, end=1.0, distr_type="normal"):
    color_dict = {}
    if distr_type == "normal":
        x_step = (end * no_values - begin * no_values) / no_values
        x_range = [ no_values * begin + x_step * x for x in range(no_values) ]
        
        color_dict = {v : 
            np.array((
                get_max((255 - (x * (2 * 255/no_values))), 0),
                get_max(255 - get_max(255 - (x * (2 * 255/no_values)), 0) - get_max((-765/3 + (x * (2 * 255/no_values))), 0), 0),
                get_max((-765/3 + (x * (2 * 255/no_values))), 0)
            ), np.uint8) for v, x in zip(range(1, no_values + 1), x_range)
        }
    elif distr_type == "blackwhite":
        x_step = (end * no_values - begin * no_values) / no_values
        x_range = [ no_values * begin + x_step * x for x in range(no_values) ]
        
        color_dict = {i : np.array((round(i * 2.55), round(i * 2.55), round(i * 2.55)), np.uint8) for i in range(100)}

    black = np.array((0, 0, 0), np.uint8)
    red = np.array((0, 0, 255), np.uint8)
    savanna_color = np.array((170, 255, 255), np.uint8)
    if distr_type == "recruitment":
        color_dict[0] = black
        color_dict[-5] = black
        color_dict[-6] = black
        color_dict[-7] = np.array((150, 255, 255), np.uint8) # Recruitment
    elif distr_type == "fire_freq":
        color_step = 255 / (no_values + 1)
        for i in range(no_values + 1):
            color_dict[i] = np.array((0, 0, i * color_step), np.uint8)
    elif distr_type == "normal":
        color_dict[0] = savanna_color
        color_dict[-5] = savanna_color
        color_dict[-6] = savanna_color
        color_dict[-7] = savanna_color
        
    return color_dict

def visualize_image(img, image_width, interpolation="none"):
    if interpolation == "none":
        img_resized = cv2.resize(img, (image_width, image_width),
               interpolation = cv2.INTER_NEAREST)
    elif interpolation == "linear":
        img_resized = cv2.resize(img, (image_width, image_width),
               interpolation = cv2.INTER_LINEAR)

    cv2.imshow("DBR Simulation (TEST)", img_resized)
    cv2.waitKey(1)

def get_image(img, color_dict, width, height="width"):
    if height == "width":
        height = width
    if (color_dict == {0:0, 1:255}):
        img *= 255
    elif (len(color_dict[0]) == 3):
        new_img = np.zeros((width, height, 3), np.uint8)
        for i in range(width):
            for j in range(height):
                new_img[i, j] = color_dict[img[i, j]]
        img = new_img.copy()
            
    img = img.astype(np.uint8)
    
    return img

def get_image_from_grid(grid, collect_states, color_dict, invert=False):
    img = grid.get_distribution(collect_states)
    img = get_image(img, color_dict, grid.width)
    if invert:
        img = ~img
    return img

def get_fire_freq_image(fire_freq_arrays, color_dict_fire_freq, grid_width, fire_no_timesteps):
    img = np.zeros((grid_width, grid_width), np.uint8)
    for fire_freq_arr in fire_freq_arrays:
        img += fire_freq_arr
    return get_image(img, color_dict_fire_freq, grid_width) * (10.0 / fire_no_timesteps)

def visualize(grid, image_width=1000, collect_states=True, color_dict={0:0, 1:255}):
    img = get_image_from_grid(grid, collect_states, color_dict)    
    visualize_image(img, image_width)
    
    return img

def save_image(img, path, image_width = 1000):
    img_resized = cv2.resize(img, (image_width, image_width), interpolation = cv2.INTER_LINEAR)
    cv2.imwrite(path, img_resized)

def save_resource_grid_colors(dynamics, species, resource, path, image_width=1000):
    arr = dynamics.get_resource_grid_colors(species, resource, 0)
    arr -= arr.min()
    img = (arr / arr.max() * 255).astype(np.uint8)
    save_image(img, path, image_width)

def visualize_difference(image1, image2, image_width=1000):
    img = np.array([image2, image2, image1], np.uint8)
    img = img.transpose(1, 2, 0)
    
    visualize_image(img, image_width)

def visualize_kernel(kernel, title="Kernel"):
    vals = []
    no_samples = 1000000
    for i in range(no_samples):
        if i % (no_samples / 5) == 0:
            print(f"Sampling... ({i} / {no_samples})")
        val = kernel.get_dist()
        vals.append(val)
    plt.hist(vals, bins=50)
    plt.title(title)
    print("samples: ", vals[0], vals[int(no_samples/5 * 1)], vals[int(no_samples/5 * 2)], vals[int(no_samples/5 * 3)], vals[int(no_samples/5 * 4)], vals[-1])
    plt.show()

def visualize_radial_distribution_function(g_r, radii, iteration):
    fig = plt.figure()
    plt.plot(radii, g_r)
    fig.suptitle(f'Radial distribution function (iteration {iteration})', fontsize=14)
    plt.xlabel('radius (m)', fontsize=11)
    plt.ylabel('g(r)', fontsize=11)
    fig.show()
    try:
        fig.savefig(f'{DATA_OUT_DIR}/radial_distribution_function/g_r_{iteration}.png')
    except:
        pass

def visualize_legend(distr_type="normal"):
    width = 500
    height = 50
    end_value = 90 if distr_type == "normal" else 0
    vals = np.zeros((height, width))
    for i in range(width):
        vals[:, i] = end_value - int(float(i) * (end_value / width))
        print(vals[0, i])
        
    color_dict = get_color_dict(100, begin=0.2, end=0.5, distr_type=distr_type)
    img = get_image(vals, color_dict, height, width)
    cv2.imwrite(os.path.join(LEGEND_PATH, f"Legend_{distr_type}.png"), img)
    cv2.imshow(f"Legend_{distr_type}", img)
    cv2.waitKey(1)


class Graphs:
    def __init__(self, dynamics, width=6, height=4):
        self.dynamics = dynamics
        self.types = {"timeseries": ["Tree cover"]}
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
        if self.datatype == "Tree cover":
            self.data.append(self.dynamics.state.grid.get_tree_cover())
        elif self.datatype == "Seeds dispersed":
            self.data.append(self.dynamics.seeds_produced)
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
        self.rect_width = self.bar_container.patches[0].get_width()
        initial_hist_width = 0
        for rect in self.bar_container.patches:
            initial_hist_width += rect.get_width()
        self.width_over_range = initial_hist_width / initial_xrange
 
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
        rect_width = ((xrange * self.width_over_range) / self.no_bins) * 4.0
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
     



    
    
