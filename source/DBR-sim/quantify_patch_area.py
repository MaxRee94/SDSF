"""Quantify patch areas in binary images with periodic boundary conditions."""

import visualization as vis
from config import *
import cv2
import numpy as np
import statistics as stats


def get_patch_areas(image, domain_width):
    """
    Count and measure connected patches in a binary image with periodic boundary conditions.

    Parameters
    ----------
    image : numpy.ndarray
        Binary image with 1s for patches and 0s for background.
    domain_width : float
        Width of the domain in real-world units.

    Returns
    -------
    list
        List of patch areas in real-world units.
    """
    # Get image dimensions
    height, width = image.shape

    # Create a label matrix for connected components
    labels = np.zeros_like(image, dtype=int)
    
    # Starting label for connected components
    current_label = 1

    # Get area per cell
    cell_width = domain_width / 200
    cell_area = cell_width * cell_width

    def flood_fill(x, y, label):
        """
        Flood fill algorithm to label connected components.

        Uses a stack to fill the component starting from (x, y) with periodic boundary conditions.

        Parameters
        ----------
        x : int
            Starting x-coordinate.
        y : int
            Starting y-coordinate.
        label : int
            Label to apply to the connected component.

        Returns
        -------
        float
            Area of the patch in real-world units.
        """
        stack = [(x, y)]
        patch_area = 1
        while stack:
            cx, cy = stack.pop()
            if labels[cx, cy] == 0 and image[cx, cy] == 1:
                labels[cx, cy] = label
                # Add all 4 neighbors (including wraparound)
                for nx, ny in [(cx-1, cy), (cx+1, cy), (cx, cy-1), (cx, cy+1)]:
                    # Apply periodic boundary conditions
                    nx, ny = nx % height, ny % width
                    if labels[nx, ny] == 0 and image[nx, ny] == 1:
                        stack.append((nx, ny))
                        patch_area += cell_area
        
        return patch_area

    # Loop through each pixel in the original image (without padding)
    patch_areas = []
    for i in range(height):
        for j in range(width):
            if image[i, j] == 1 and labels[i, j] == 0:
                # If it's a patch and not labeled, do a flood fill
                patch_area = flood_fill(i, j, current_label)
                patch_areas.append(patch_area)
                current_label += 1
    
    return patch_areas

def quantify_patch_area_distribution(patch_image, domain_width):
    """
    Quantify the distribution of patch areas in a binary image.

    Parameters
    ----------
    patch_image : numpy.ndarray
        Binary image with patches.
    domain_width : float
        Width of the domain in real-world units.

    Returns
    -------
    list
        List of patch areas.
    """
    # Count the patches with periodic boundary conditions
    patch_areas = get_patch_areas(patch_image, domain_width)

    return patch_areas


def generate_patch_image(patch_width, treecover=0.5, i=0):
    """
    Generate a thresholded Perlin noise image for patch area quantification.

    Parameters
    ----------
    patch_width : float
        Target patch width.
    treecover : float, optional
        Target tree cover fraction (default 0.5).
    i : int, optional
        Trial index for file naming (default 0).

    Returns
    -------
    numpy.ndarray
        Thresholded binary image.
    """
    # Generate unthresholded perlin noise image
    path = f"{cfg.PERLIN_NOISE_DIR}/perlin_noise_for_area_quantification.png"
    noise_frequency = 5.0 / patch_width # Convert patch width to noise frequency
    noise_frequency = round(noise_frequency, 2) # Conform noise frequency to 2 decimal places, to ensure periodicity of the noise pattern
    img = vis.generate_perlin_noise_image(path, frequency=noise_frequency, octaves=5, write=True)

    # Threshold image
    img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    img = vis.get_thresholded_image(img, treecover * img.shape[0] * img.shape[0] * 255 )
    thresholded_path = path.replace(".png", f"_thresholded_{i}.png")
    cv2.imwrite(thresholded_path, img)

    return img

def generate_and_load_patch_image(patch_width, i):
    """
    Generate and normalize a patch image.

    Parameters
    ----------
    patch_width : float
        Target patch width.
    i : int
        Trial index.

    Returns
    -------
    numpy.ndarray
        Normalized image with values in [0, 1].
    """
    img = generate_patch_image(patch_width, i=i)
    img = img / 255
    return img

def main(patch_width, num_trials=100, domain_width=960):
    """
    Run patch area quantification for multiple trials.

    Parameters
    ----------
    patch_width : float
        Target patch width.
    num_trials : int, optional
        Number of trials to run (default 100).
    domain_width : float, optional
        Width of the domain in real-world units (default 960).
    """
    patch_areas = []
    for i in range(num_trials):
        print(f"Processing pattern {i+1} / {num_trials}...")
        patch_image = generate_and_load_patch_image(patch_width, i)
        patch_areas += quantify_patch_area_distribution(patch_image, domain_width)
    
    domain_area = domain_width * domain_width
    print("-"*40)
    print("Forest patch shape distribution as percentage of domain area:")
    print("Mean:", round(100 * stats.mean(patch_areas) / domain_area, 5))
    print("Stdev:", round(100 * stats.stdev(patch_areas) / domain_area, 5))
    print("-"*40)

if __name__ == "__main__":
    main(167, num_trials=100, domain_width=960)
    

