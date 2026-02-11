from multiprocessing import Value
import cv2
import numpy as np
import matplotlib.pyplot as plt
import math
import time
import colorsys
import platform
import subprocess
import re
import ctypes

import disk_pattern_generator as dpg
from helpers import *
from config import *



class Visualiser:

    def __init__(self, cfg):
        self.cfg = cfg
        self.init_perlin_noise_global_variables()
        self.dpg = dpg.DiskPatternGenerator(cfg)

    def init_perlin_noise_global_variables(self):
        # Perlin noise global variables
        self.perm = list(range(4096))
        self.cfg.rng.shuffle(self.perm)
        self.perm += self.perm
        self.dirs = [(math.cos(a * 2.0 * math.pi / 4096),
                    math.sin(a * 2.0 * math.pi / 4096))
                    for a in range(4096)]

    def able_to_read_and_resize_image(self, path):
        try:
            img = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
            if img is None:
                return False
            img2 = cv2.resize(img, (1000, 1000), interpolation=cv2.INTER_LINEAR)
            return True
        except Exception as e:
            print(f"Error reading image {path}: {e}")
            return False

    def generate_disk_pattern(self, initial_pattern_image, write=True, recursing=False, **kwargs):
        path = f"{cfg.CONTROLLED_PATTERN_DIR}/" + initial_pattern_image + ".png"
        img, positions, radii, stripe_metadata, benchmark_cover = self.dpg.create_image(**kwargs)
        cv2.imwrite(path, img)
        if recursing:
            time.sleep(1)
        print("Generated controlled pattern image at ", path)
        if not self.able_to_read_and_resize_image(path):
            print("Error: Image failed to save or read. Re-trying..")
            time.sleep(1) # Wait before trying to load again
            if self.able_to_read_and_resize_image(path):
                return img, path

            if recursing:
                raise RuntimeError("Initial image generation failed after re-trying.")
            else:
                return self.generate_disk_pattern(initial_pattern_image, self.ctrl_pattern_generator_params, write, recursing=True)
        return img, path, benchmark_cover
    
    def reshuffle_perlin_noise(self):
        global perm
        perm = list(range(4096))
        self.cfg.rng.shuffle(perm)
        perm += perm

    def number_to_rgb(self, n: int, i: int) -> tuple[int, int, int]:
        # Normalize i to [0,1]
        hue = i / float(n) 
        # Use full saturation and brightness
        r, g, b = colorsys.hsv_to_rgb(hue, 1, 1)
        # Scale to 0�255
        return int(r * 255), int(g * 255), int(b * 255)

    def get_random_color_index(self, existing_colors, max_no_colors, offset):
        color_idx = offset - self.cfg.rng.integers(0, max_no_colors-1)
        if len(existing_colors) >= max_no_colors:
            # All possible colors are taken; simply return a random one (we assume it's okay that there will be duplicates)
            return color_idx

        max_no_attempts = 1000
        i = 0
        while color_idx in existing_colors and i < max_no_attempts:
            color_idx = offset - self.cfg.rng.integers(0, max_no_colors-1)
            i+=1

        return color_idx

    def get_most_distinct_index(self, existing_colors, maximum_possible_no_colors, offset):
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

    def get_color_dict(self, no_values, begin=0.0, end=1.0, distr_type="normal"):
        color_dict = {}
        black = np.array((0, 0, 0), np.uint8)
        red = np.array((0, 0, 255), np.uint8)
        purple = np.array((255, 0, 255), np.uint8)
        white = np.array((255, 255, 255), np.uint8)
        savanna_color = np.array((170, 255, 255), np.uint8)

        if distr_type == "normal" or distr_type == "normal_with_fire_effects":
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
            color_dict = {i : np.array((round(i * 2.55), round(i * 2.55), round(i * 2.55)), np.uint8) for i in range(100)}
            color_dict[-1] = red
        if distr_type == "recruitment":
            color_dict[0] = black
            color_dict[-5] = black
            color_dict[-6] = black
            color_dict[-7] = np.array((150, 255, 255), np.uint8) # Recruitment
            color_dict[-8] = black
        elif distr_type == "fire_freq":
            color_step = 255 / (no_values + 1)
            for i in range(no_values + 1):
                color_dict[i] = np.array((0, 0, i * color_step), np.uint8)
        elif distr_type == "normal_with_fire_effects":
            color_dict[0] = savanna_color
            color_dict[-5] = red
            color_dict[-6] = black
            color_dict[-7] = savanna_color
            color_dict[-8] = purple
        elif distr_type == "normal":
            color_dict[0] = savanna_color
            color_dict[-5] = savanna_color
            color_dict[-6] = savanna_color
            color_dict[-7] = savanna_color
            color_dict[-8] = savanna_color
        elif distr_type == "colored_patches":
            # Set all colors between -9 and 0 to savanna color; these are not to be used for patches
            for i in range(10):
                color_dict[-i] = savanna_color

            # Now assign distinct colors for patches from -10 downwards
            max_no_patches = 1000
            for i in range(max_no_patches): 
                color_dict[-10 - i] = self.number_to_rgb(max_no_patches, i)
        
        return color_dict

    def visualize_image(self, img, image_width, interpolation="none"):
        if interpolation == "none":
            img_resized = cv2.resize(img, (image_width, image_width),
                   interpolation = cv2.INTER_NEAREST)
        elif interpolation == "linear":
            img_resized = cv2.resize(img, (image_width, image_width),
                   interpolation = cv2.INTER_LINEAR)

        cv2.imshow("DBR Simulation (TEST)", img_resized)
        cv2.waitKey(1)

    def get_image(self, img, color_dict, width, height="width"):
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

    def normalize_distribution(self, distr, _range):
        # Normalize the given distribution and then make it conform to the given value range.

        distr_min = distr.min()
        distr_max = distr.max()
        normalized_distr = (distr - distr_min) / (distr_max - distr_min)
        distr = normalized_distr * (_range[1] - _range[0]) + _range[0]
        return distr

    def get_image_from_grid(self, grid, color_dict, collect_states=None, img_type=None, invert=False):
        if img_type:
            if img_type == "fuel":
                img = grid.get_fuel_distribution()
                img = (img * 99).astype(int) # Normalize to 0-99
            elif img_type == "aggr_tree_LAI":
                img = grid.get_aggr_tree_LAI_distribution()
                img = (img * 14).astype(int) # Assume max LAI is ~7, normalize to 0-99
                img.clip(0, 99, out=img)
            elif img_type == "fuel_penetration":
                fuel = grid.get_fuel_distribution()
                img = (fuel * 99).astype(int) # Normalize to 0-99
                non_aggregated_tree_LAI = grid.get_distribution(1)
                mask = non_aggregated_tree_LAI < 1 # Mask area considered to be savanna.
                img[mask] = 0
            else:
                raise ValueError(f"Unknown img_type: {img_type}")
        else:
            img = grid.get_distribution(collect_states)
        img = self.get_image(img, color_dict, grid.width)

        if invert:
            img = ~img
        return img

    def get_fire_freq_image(self, fire_freq_arrays, color_dict_fire_freq, grid_width, fire_no_timesteps):
        img = np.zeros((grid_width, grid_width), np.uint8)
        for fire_freq_arr in fire_freq_arrays:
            img += fire_freq_arr
        return self.get_image(img, color_dict_fire_freq, grid_width)

    def visualize(self, grid, image_width=1000, collect_states=True, color_dict=False):
        img = self.get_image_from_grid(grid, color_dict, collect_states=collect_states)   
        _, screen_height = self.get_screen_resolution()
        image_width = int(screen_height * 0.6)
        self.visualize_image(img, image_width)
    
        return img

    def _get_windows_resolution(self):
        # Win32 API via ctypes
        user32 = ctypes.windll.user32
        user32.SetProcessDPIAware()  # avoid DPI scaling issues
        width = user32.GetSystemMetrics(0)
        height = user32.GetSystemMetrics(1)
        return width, height


    def _get_macos_resolution(self):
        # Query system_profiler and parse the resolution
        output = subprocess.check_output(
            ["system_profiler", "SPDisplaysDataType"], text=True
        )
        match = re.search(r"Resolution:\s*(\d+)\s*x\s*(\d+)", output)
        if not match:
            raise RuntimeError("Could not determine macOS screen resolution")
        return int(match.group(1)), int(match.group(2))


    def get_screen_resolution(self):
        system = platform.system()

        if system == "Windows":
            return self._get_windows_resolution()

        if system == "Darwin":
            return self._get_macos_resolution()

        raise NotImplementedError(f"Unsupported OS: {system}")

    def save_image(self, img, path, image_width = 1000, interpolation="none"):
        if interpolation == "linear":
            img_resized = cv2.resize(img, (image_width, image_width), interpolation = cv2.INTER_LINEAR)
        elif interpolation == "none":
            img_resized = cv2.resize(img, (image_width, image_width), interpolation = cv2.INTER_NEAREST)
        cv2.imwrite(path, img_resized)

    def save_resource_grid_colors(self, dynamics, species, resource, path, image_width=1000):
        arr = dynamics.get_resource_grid_colors(species, resource, 0)
        arr -= arr.min()
        img = (arr / arr.max() * 255).astype(np.uint8)
        self.save_image(img, path, image_width)

    def visualize_difference(self, image1, image2, image_width=1000):
        img = np.array([image2, image2, image1], np.uint8)
        img = img.transpose(1, 2, 0)
    
        self.visualize_image(img, image_width)

    def get_thresholded_image(self, orig_img, white_pixel_sum):
        threshold = 245
        ret, img = cv2.threshold(orig_img, threshold, 255, cv2.THRESH_BINARY) 
        while np.sum(img) < white_pixel_sum:
            threshold -= 1
            ret, img = cv2.threshold(orig_img, threshold, 255, cv2.THRESH_BINARY)
    
        return img

    def perlin_noise(self, x, y, per):
        def surflet(gridX, gridY):
            distX, distY = abs(x-gridX), abs(y-gridY)
            polyX = 1 - 6*distX**5 + 15*distX**4 - 10*distX**3
            polyY = 1 - 6*distY**5 + 15*distY**4 - 10*distY**3
            hashed = perm[perm[int(gridX)%per] + int(gridY)%per]
            grad = (x-gridX)*self.dirs[hashed][0] + (y-gridY)*self.dirs[hashed][1]
            return polyX * polyY * grad
    
        intX, intY = int(x), int(y)
        return (surflet(intX+0, intY+0) + surflet(intX+1, intY+0) +
                surflet(intX+0, intY+1) + surflet(intX+1, intY+1))

    def fBm(self, x, y, per, octs):
        val = 0
        for o in range(octs):
            val += 0.5**o * self.perlin_noise(x*2**o, y*2**o, per*2**o)
        return val

    def generate_perlin_noise_image(self, path, width=200, frequency=1/32.0, octaves=5, write=True):
        self.reshuffle_perlin_noise()
        data = []
        for y in range(width):
            row = []
            for x in range(width):
                val = self.fBm(x*frequency, y*frequency, int(width*frequency), octaves)
                val = min(255, max(0, (val + 0.5) * 255))
                row.append(val)
            data.append(row)
        img = np.array(data, dtype=np.uint8)
        if write:
            cv2.imwrite(path, img)

        return img

    def visualize_kernel(self, kernel, title="Kernel"):
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

    def visualize_radial_distribution_function(self, g_r, radii, iteration):
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

    def visualize_legend(self, distr_type="normal"):
        width = 500
        height = 50
        end_value = 80 if distr_type == "normal" else 0
        vals = np.zeros((height, width))
        for i in range(width):
            vals[:, i] = end_value - int(float(i) * (end_value / width))
            print(vals[0, i])
        
        color_dict = self.get_color_dict(100, begin=0.2, end=0.5, distr_type=distr_type)
        img = self.get_image(vals, color_dict, height, width)
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
        # loop until all UI events currently waiting have been processed
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
     

def do_visualizations(dynamics, fire_freq_arrays, fire_no_timesteps, verbose, color_dicts, collect_states, visualization_types, patches, patch_count_change, patch_color_ids, cfg):
    if ("recruitment" in visualization_types):
        print("Saving recruitment img...") if verbose else None
        recruitment_img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["recruitment"], collect_states=0)
        imagepath_recruitment = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/recruitment/" + str(dynamics.time) + ".png")
        cfg.vis.save_image(recruitment_img, imagepath_recruitment, get_max(1000, recruitment_img.shape[0]), interpolation="none")
    
    if ("fire_freq" in visualization_types):
        fire_freq_arrays.append(dynamics.state.grid.get_distribution(0) == -5)
        if dynamics.time > fire_no_timesteps:
            print("Saving fire frequency img...") if verbose else None
            fire_freq_img = cfg.vis.get_fire_freq_image(fire_freq_arrays[-fire_no_timesteps:], color_dicts["fire_freq"], dynamics.state.grid.width, fire_no_timesteps)
            imagepath_fire_freq = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fire_frequencies/" + str(dynamics.time) + ".png")
            cfg.vis.save_image(fire_freq_img, imagepath_fire_freq, get_max(1000, fire_freq_img.shape[0]))

    def get_bbox(positions2d):
        min_x = min(positions2d, key=lambda item: item[0])[0]
        max_x = max(positions2d, key=lambda item: item[0])[0]
        min_y = min(positions2d, key=lambda item: item[1])[1]
        max_y = max(positions2d, key=lambda item: item[1])[1]
        return (min_x, min_y), (max_x, max_y)

    if ("colored_patches" in visualization_types):
        print("Creating colored patches image...") if verbose else None

        # Initialize color index array
        patch_colors_indices = np.zeros((dynamics.state.grid.width, dynamics.state.grid.width), dtype=int) -1
        
        # Sort patches so that larger patches are assigned colors first
        patches = sorted(patches, key=lambda patch: patch["area"], reverse=True)

        minimum_considered_patch_size = 0 # m^2
        forest_area = 0
        central_patch = -1
        max_no_colors = 1000
        offset = -10

        for i, patch in enumerate(patches):
            if (len(patches) > 20) and dynamics.time == 0:
                if i >= 100:
                    print("Breaking patch coloring loop at patch no", i, f"for performance reasons (number of patches = {len(patches)}).")
                    break
                if i % 15 == 0:
                    print(f"Coloring patch no {i}/{len(patches)}... Will be cut off at 100.")
            if patch["area"] < minimum_considered_patch_size: # Only consider patches of a certain size
                continue
            patch_id = patch["id"]

            if not patch_color_ids.get(str(patch_id)):
                if (len(patches) > 30) and dynamics.time > 0:
                    color_idx = cfg.vis.get_random_color_index(patch_color_ids.values(), max_no_colors, offset)
                else:
                    color_idx = cfg.vis.get_most_distinct_index(patch_color_ids.values(), max_no_colors, offset)
                patch_color_ids[str(patch_id)] = color_idx # Assign a new color index to the patch
            patch_color_id = patch_color_ids[str(patch_id)]
            col=color_dicts["colored_patches"][patch_color_id]
            forest_area += patch["area"] if patch["type"] == "forest" else 0
            for cell in patch["cells"]:
                patch_colors_indices[cell[1]][cell[0]] = patch_color_id

        colored_patches_img = cfg.vis.get_image(patch_colors_indices, color_dicts["colored_patches"], dynamics.state.grid.width)
        cfg.show_edges = False # Hardcoded for now
        if cfg.show_edges:
            for patch in patches:
                if patch["area"] < minimum_considered_patch_size: # Only consider patches of a certain size
                    continue
                for edge in patch["perimeter"]:
                    point1 = (edge[0][0], edge[0][1])
                    point2 = (edge[1][0], edge[1][1])
                    if get_2d_dist(point1, point2) > 0.5*dynamics.state.grid.width:
                        continue # Skip drawing wrap-around edges for now.
                    cv2.line(colored_patches_img, point1, point2, (0, 0, 0), 1)  # black line, thickness=2
        imagepath_colored_patches = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/colored_patches/" + str(dynamics.time) + ".png")
        cfg.vis.save_image(colored_patches_img, imagepath_colored_patches, get_max(1000, colored_patches_img.shape[0]), interpolation="none")

    print("-- Visualizing image...") if verbose else None
    if cfg.headless:
        # Get a color image representation of the state
        img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["normal"], collect_states=1)
    else:
        # Get a color image representation of the state and show it.
        img = cfg.vis.visualize(
            dynamics.state.grid, cfg.image_width, collect_states=collect_states,
            color_dict=color_dicts["normal"]
        )

    print("-- Saving image...") if verbose else None
    imagepath = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/" + str(dynamics.time) + ".png")
    cfg.vis.save_image(img, imagepath, get_max(1000, img.shape[0]), interpolation="none")
    
    if ("fuel" in visualization_types):
        print("-- Saving fuel image...") if verbose else None
        fuel_img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["blackwhite"], img_type="fuel")
        imagepath_fuel = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fuel/" + str(dynamics.time) + ".png")
        cfg.vis.save_image(fuel_img, imagepath_fuel, get_max(1000, fuel_img.shape[0]))

    if ("aggr_tree_LAI" in visualization_types):
        print("-- Saving aggregated tree LAI image...") if verbose else None
        aggr_tree_LAI_img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["blackwhite"], img_type="aggr_tree_LAI", invert=False)
        imagepath_aggr_tree_LAI = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/aggr_tree_LAI/" + str(dynamics.time) + ".png")
        cfg.vis.save_image(aggr_tree_LAI_img, imagepath_aggr_tree_LAI, get_max(1000, aggr_tree_LAI_img.shape[0]))

    if ("fuel_penetration" in visualization_types):
        print("-- Saving fuel penetration image...") if verbose else None
        fuel_penetration_img = cfg.vis.get_image_from_grid(dynamics.state.grid, color_dicts["blackwhite"], img_type="fuel_penetration")
        imagepath_fuel_penetration = os.path.join(cfg.DATA_OUT_DIR, "image_timeseries/fuel_penetration/" + str(dynamics.time) + ".png")
        cfg.vis.save_image(fuel_penetration_img, imagepath_fuel_penetration, get_max(1000, fuel_penetration_img.shape[0]))

    
    
