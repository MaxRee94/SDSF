from multiprocessing import Value
import cv2
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import colorsys
import random
import controllable_pattern_generator as cpg
from PIL import Image

from helpers import *
from config import *


# Perlin noise global variables
perm = list(range(4096))
random.shuffle(perm)
perm += perm
dirs = [(math.cos(a * 2.0 * math.pi / 4096),
         math.sin(a * 2.0 * math.pi / 4096))
            for a in range(4096)]


def able_to_read_and_resize_image(path):
    try:
        img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
        if img is None:
            return False
        img2 = cv2.resize(img, (1000, 1000), interpolation=cv2.INTER_LINEAR)
        return True
    except Exception as e:
        print(f"Error reading image {path}: {e}")
        return False


def generate_controllable_pattern_image(initial_pattern_image, write=True, recursing=False, **kwargs):
    path = f"{cfg.CONTROLLED_PATTERN_DIR}/" + initial_pattern_image + ".png"
    img, positions, radii, stripe_metadata, benchmark_cover = cpg.create_image(**kwargs)
    cv2.imwrite(path, img)
    if recursing:
        time.sleep(1)
    print("Generated controlled pattern image at ", path)
    if not able_to_read_and_resize_image(path):
        print("Error: Image failed to save or read. Re-trying..")
        time.sleep(1) # Wait before trying to load again
        if able_to_read_and_resize_image(path):
            return img, path

        if recursing:
            raise RuntimeError("Initial image generation failed after re-trying.")
        else:
            return generate_controllable_pattern_image(initial_pattern_image, ctrl_pattern_generator_params, write, recursing=True)
    return img, path, benchmark_cover
    

def reshuffle_perlin_noise():
    global perm
    perm = list(range(4096))
    random.shuffle(perm)
    perm += perm

def number_to_rgb(n: int, i: int) -> tuple[int, int, int]:
    # Normalize i to [0,1]
    hue = i / float(n) 
    # Use full saturation and brightness
    r, g, b = colorsys.hsv_to_rgb(hue, 1, 1)
    # Scale to 0�255
    return int(r * 255), int(g * 255), int(b * 255)

def get_most_distinct_index(existing_colors, maximum_possible_no_colors, offset):
    dist = 1e10
    best_dist = 0
    best_index = offset
    if not existing_colors:
        return offset
    for i in range(offset - maximum_possible_no_colors + 1, offset):
        # Get closest distance to other color indices
        for color in existing_colors:
            _dist = abs(color - i)
            if _dist < dist:
                dist = _dist
            # Account for periodic boundary conditions
            wraparound_dist = abs((color - maximum_possible_no_colors) - i)
            if wraparound_dist < dist:
                dist = wraparound_dist
        if dist > best_dist:
            best_dist = dist
            best_index = i
        dist = 1e10

    return best_index

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
        color_dict[-5] = red
        color_dict[-6] = red
        color_dict[-7] = savanna_color
    elif distr_type == "colored_patches":
        color_dict[0] = savanna_color
        color_dict[-5] = savanna_color
        color_dict[-6] = savanna_color
        color_dict[-7] = savanna_color
        max_no_patches = 100

        for i in range(max_no_patches):         
            color_dict[-10 - i] = number_to_rgb(max_no_patches, i)
        
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
    if (color_dict is False):
        color_dict = {0: 0, 1: 255}
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
    return get_image(img, color_dict_fire_freq, grid_width)

def visualize(grid, image_width=1000, collect_states=True, color_dict=False):
    img = get_image_from_grid(grid, collect_states, color_dict)    
    visualize_image(img, image_width)
    
    return img

def save_image(img, path, image_width = 1000, interpolation="linear"):
    if interpolation == "linear":
        img_resized = cv2.resize(img, (image_width, image_width), interpolation = cv2.INTER_LINEAR)
    elif interpolation == "none":
        img_resized = cv2.resize(img, (image_width, image_width), interpolation = cv2.INTER_NEAREST)
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

def get_thresholded_image(orig_img, white_pixel_sum):
    threshold = 245
    ret, img = cv2.threshold(orig_img, threshold, 255, cv2.THRESH_BINARY) 
    while np.sum(img) < white_pixel_sum:
        threshold -= 1
        ret, img = cv2.threshold(orig_img, threshold, 255, cv2.THRESH_BINARY)
    
    return img

def perlin_noise(x, y, per):
    def surflet(gridX, gridY):
        distX, distY = abs(x-gridX), abs(y-gridY)
        polyX = 1 - 6*distX**5 + 15*distX**4 - 10*distX**3
        polyY = 1 - 6*distY**5 + 15*distY**4 - 10*distY**3
        hashed = perm[perm[int(gridX)%per] + int(gridY)%per]
        grad = (x-gridX)*dirs[hashed][0] + (y-gridY)*dirs[hashed][1]
        return polyX * polyY * grad
    
    intX, intY = int(x), int(y)
    return (surflet(intX+0, intY+0) + surflet(intX+1, intY+0) +
            surflet(intX+0, intY+1) + surflet(intX+1, intY+1))

def fBm(x, y, per, octs):
    val = 0
    for o in range(octs):
        val += 0.5**o * perlin_noise(x*2**o, y*2**o, per*2**o)
    return val

def generate_perlin_noise_image(path, width=200, frequency=1/32.0, octaves=5, write=True):
    reshuffle_perlin_noise()
    data = []
    for y in range(width):
        row = []
        for x in range(width):
            val = fBm(x*frequency, y*frequency, int(width*frequency), octaves)
            val = min(255, max(0, (val + 0.5) * 255))
            row.append([val, val, val])
        data.append(row)
    img = np.array(data, dtype=np.uint8)
    if write:
        cv2.imwrite(path, img)

    return img

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
        fig.savefig(f'{cfg.DATA_OUT_DIR}/radial_distribution_function/g_r_{iteration}.png')
    except:
        pass

def visualize_legend(distr_type="normal"):
    width = 500
    height = 50
    end_value = 80 if distr_type == "normal" else 0
    vals = np.zeros((height, width))
    for i in range(width):
        vals[:, i] = end_value - int(float(i) * (end_value / width))
        print(vals[0, i])
        
    color_dict = get_color_dict(100, begin=0.2, end=0.5, distr_type=distr_type)
    img = get_image(vals, color_dict, height, width)
    cv2.imwrite(os.path.join(cfg.LEGEND_PATH, f"Legend_{distr_type}.png"), img)
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
     



    
    
