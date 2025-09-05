import cv2
import numpy as np
import pandas as pd
from scipy.spatial import ConvexHull
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import json
import math
import os
import random
from scipy.optimize import newton
from types import SimpleNamespace
from config import *


cpgn = SimpleNamespace(area_normalization_factor=None)

def compute_area(contour, center):
    """
    Calculate the area of a disk.
    """
    if len(contour) < 3:
        raise ValueError("A disk must have at least 3 points")

    area = 0
    n = len(contour)
    c_x, c_y = center
    for i in range(n):
        x1, y1 = contour[i]
        x2, y2 = contour[(i + 1) % n]  # wrap around to the first point
        tri_area = 0.5 * abs(x1 * y2 - x1 * c_y + x2 * c_y - x2 * y1 + c_x * y1 - c_x * y2)
        area += tri_area

    return cv2.contourArea(contour.astype(np.float32))  # Use OpenCV's contour area function for accuracy

def fraction_white_pixels(img):
    
    # White pixel = 255 in grayscale
    white_pixels = np.sum(img == 255)
    total_pixels = img.size

    fraction = white_pixels / total_pixels
    return fraction

def compute_area_normalization_factor(contour, center, base_radius, global_area_normalization_factor):
    area = compute_area(contour, center)
    circle_area = math.pi * base_radius * base_radius
    if global_area_normalization_factor is not None:
        # Use global normalization factor if provided, to correct for small consistencies in the area of the generated pattern.
        cpgn.area_normalization_factor = 1.0 / (math.sqrt(circle_area / (area * max(1e-10, global_area_normalization_factor))))
    else:
        cpgn.area_normalization_factor = 1.0 / (math.sqrt(circle_area / area))

    return cpgn.area_normalization_factor

def normalize_disk(disk, center, base_radius, norm_factor):
    disk = disk - np.array(center)
    disk = disk / norm_factor
    disk = disk + np.array(center)

    #print("circle area:", math.pi * base_radius * base_radius)
    #print("disk area:", compute_area(disk, center))

    return disk.astype(np.int32)

def apply_jitter(positions, image_size, mean_distance, cv_distance):
    """
    Applies random 2D jitter to each position with std = cv_distance * mean_distance.
    Enforces periodic boundary conditions.
    """
    stdev = mean_distance * cv_distance
    box = np.array(image_size)
    jittered = []

    for pos in positions:
        jitter = np.random.normal(0, stdev, 2)
        new_pos = (np.array(pos) + jitter) % box
        jittered.append(tuple(new_pos))

    return jittered

def generate_disk(
        center, base_radius, amp1=0, wave1=1, amp2=0, wave2=2, index=None, rotate_randomly=True, resolution=360,
        global_area_normalization_factor=None, global_rotation_offset=None
    ):
    angles = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    rotational_offset = 0
    if rotate_randomly:
        if global_rotation_offset is not None:
            index *= global_rotation_offset
        random.seed(index)  # Ensure reproducibility for the same index
        rotational_offset = random.uniform(0, 2*math.pi)
    elif global_rotation_offset > 0:
        rotational_offset = global_rotation_offset
    r = base_radius + amp1 * np.sin(wave1 * angles + rotational_offset) + amp2 * np.sin(wave2 * angles + rotational_offset)
    x = center[0] + r * np.cos(angles)
    y = center[1] + r * np.sin(angles)
    disk = np.stack((x, y), axis=-1).astype(np.int32)
    norm_factor = compute_area_normalization_factor(disk, center, base_radius, global_area_normalization_factor)
    disk = normalize_disk(disk, center, base_radius, norm_factor)
    return disk

def draw_disk(img, center, radius, amp1, wave1, amp2, wave2, rotate_randomly, index, global_area_normalization_factor, global_rotation_offset=None):
    contour = generate_disk(center, radius, amp1=amp1, wave1=wave1, amp2=amp2, wave2=wave2, rotate_randomly=rotate_randomly, index=index, 
                            global_area_normalization_factor=global_area_normalization_factor, global_rotation_offset=global_rotation_offset)
    cv2.fillPoly(img, [contour], 255)

def draw_stripe(img, center1, center2, radius, rotate_randomly):
    disk1 = generate_disk(center1, radius, rotate_randomly=rotate_randomly)
    disk2 = generate_disk(center2, radius, rotate_randomly=rotate_randomly)
    points = np.vstack([disk1, disk2])
    hull = ConvexHull(points)
    hull_points = points[hull.vertices]
    cv2.fillPoly(img, [hull_points], 255)

def adjust_mean_distance_for_uniform_hex_grid(image_size, target_distance, max_rows=1000, max_cols=1000, tolerance=0.01):
    """
    Adjust mean_distance so that a hexagonal grid of disks fits exactly in the image,
    with all neighbor distances within ±1% of the adjusted mean distance.
    """
    w, h = image_size
    vertical_distance = (math.sqrt(3.0) / 2.0) * target_distance # Height of an equilateral triangle
    target_no_hexagons_along_v_axis = int(h / (vertical_distance * 2.0))  # Two rows to get a full hexagon height
    new_hexagon_height = h / target_no_hexagons_along_v_axis
    new_target_distance = (new_hexagon_height / 2.0) / (math.sqrt(3.0) / 2.0) # reversed formula for height of equilateral triangle

    return new_target_distance

def adjust_mean_distance_for_uniform_square_grid(image_size, target_distance, max_rows=1000, max_cols=1000, tolerance=0.01):
    """
    Adjust mean_distance so that a square grid of disks fits exactly in the image,
    with all neighbor distances within ±1% of the adjusted mean distance.
    """
    w, h = image_size
    target_no_squares = int(w / target_distance)
    new_target_distance = w / target_no_squares

    return new_target_distance

def generate_square_grid(image_size, mean_distance):
    """
    Generates a square grid of disks without any filtering.
    """
    w, h = image_size
    dx = dy = mean_distance
    box = np.array([w, h])

    nx = int(np.ceil(w / dx))
    ny = int(np.ceil(h / dy))

    positions = []
    for j in range(ny):
        for i in range(nx):
            x = i * dx
            y = j * dy
            positions.append((x % w, y % h))

    return positions


def generate_hex_grid_with_filter(image_size, mean_distance):
    """
    Generate a hexagonal grid and remove overlaps.
    """
    w, h = image_size
    dx = mean_distance
    dy = mean_distance * np.sqrt(3) / 2
    box = np.array([w, h])

    nx = int(np.ceil(w / dx)) + 3
    ny = int(np.ceil(h / dy)) + 3

    raw_positions = []
    for j in range(ny):
        for i in range(nx):
            x = i * dx + (dx / 2 if j % 2 else 0)
            y = j * dy
            pos = np.array([x, y]) % box
            raw_positions.append(pos)

    filtered_positions = []
    for pos in raw_positions:
        too_close = False
        for prev in reversed(filtered_positions):
            delta = np.abs(pos - np.array(prev))
            delta = np.where(delta > 0.5 * box, box - delta, delta)
            dist = np.linalg.norm(delta)
            if dist < (0.98 * mean_distance):  # Allow a small tolerance
                too_close = True
                break
        if not too_close:
            filtered_positions.append(tuple(pos))

    return filtered_positions

def build_parallel_stripes(positions, image_size, stripe_angle_deg, mean_length, std_length):
    """
    Return stripe pairs with a fixed direction and approximately fixed length.
    """
    angle_rad = np.deg2rad(stripe_angle_deg)
    dir_vec = np.array([np.cos(angle_rad), np.sin(angle_rad)])
    used = set()
    stripes = []
    stripe_metadata = []

    # Project each position onto the stripe axis
    proj = [np.dot(pos, dir_vec) for pos in positions]
    sort_idx = np.argsort(proj)
    positions_sorted = [positions[i] for i in sort_idx]

    for i in range(len(positions_sorted) - 1):
        if i in used:
            continue
        for j in range(i + 1, len(positions_sorted)):
            if j in used:
                continue
            p1 = np.array(positions_sorted[i])
            p2 = np.array(positions_sorted[j])
            vec = p2 - p1
            length = np.linalg.norm(vec)
            angle_diff = np.abs(np.arccos(np.dot(vec, dir_vec) / (np.linalg.norm(vec) + 1e-9)))
            if angle_diff < 0.1 and abs(length - mean_length) < 3 * std_length:
                used.add(i)
                used.add(j)
                stripes.append((tuple(p1), tuple(p2)))
                stripe_metadata.append({
                    "x1": p1[0], "y1": p1[1],
                    "x2": p2[0], "y2": p2[1],
                    "length": length
                })
                break

    return stripes, stripe_metadata

def export_metadata(positions, radii, parameters, stripe_metadata=None, output_dir=cfg.CPG_OUTPUT_DIR):
    os.makedirs(output_dir, exist_ok=True)
    df = pd.DataFrame({
        "x": [p[0] for p in positions],
        "y": [p[1] for p in positions],
        "radius": radii
    })
    df.to_csv(os.path.join(output_dir, "disk_metadata.csv"), index=False)

    with open(os.path.join(output_dir, "disk_parameters.json"), "w") as f:
        json.dump(parameters, f, indent=4)

    if stripe_metadata:
        df_stripes = pd.DataFrame(stripe_metadata)
        df_stripes.to_csv(os.path.join(cfg.CPG_OUTPUT_DIR, "stripe_metadata.csv"), index=False)

def draw_sinusoidal_stripe(img, p1, p2, radius, amplitude, wavelength, n_points=100):
    """
    Draws a sinusoidal stripe (a ribbon-like convex hull) between p1 and p2.
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    vec = p2 - p1
    length = np.linalg.norm(vec)
    if length < 1e-3:
        return

    direction = vec / length
    normal = np.array([-direction[1], direction[0]])

    # Interpolate along the center line
    t = np.linspace(0, 1, n_points)
    line = p1[np.newaxis, :] + t[:, np.newaxis] * vec[np.newaxis, :]

    # Sinusoidal offset perpendicular to direction
    phase = 2 * np.pi * t * length / wavelength
    offset = amplitude * np.sin(phase)[:, np.newaxis] * normal[np.newaxis, :]

    upper = line + radius * direction + offset
    lower = line - radius * direction + offset[::-1]

    ribbon = np.vstack([upper, lower]).astype(np.int32)
    cv2.fillPoly(img, [ribbon], 255)


def create_image(**kwargs):
    args = SimpleNamespace(**kwargs)
    args.global_rotation_offset = random.randint(0, 100000)
    
    if args.grid_type == "hex":
        if args.enforce_distance_uniformity:
            original_distance = args.mean_distance
            args.mean_distance = adjust_mean_distance_for_uniform_hex_grid(args.image_size, args.mean_distance)
            if not math.isclose(original_distance, args.mean_distance, abs_tol=1):
                if args.suppress_distance_warning:
                    print(f"Warning: Mean_distance {original_distance:.2f} adjusted to {args.mean_distance:.2f} for uniform hex grid.")
                else:
                    raise ValueError(f"Mean_distance {original_distance:.2f} needs to be changed to {args.mean_distance-0.5} to have a perfect hex grid.")
        base_positions = generate_hex_grid_with_filter(args.image_size, args.mean_distance)
    elif args.grid_type == "square":
        if args.enforce_distance_uniformity:
            original_distance = args.mean_distance
            args.mean_distance = adjust_mean_distance_for_uniform_square_grid(args.image_size, args.mean_distance)
            if not math.isclose(original_distance, args.mean_distance, abs_tol=1):
                if args.suppress_distance_warning:
                    print(f"Warning: Mean_distance {original_distance:.2f} adjusted to {args.mean_distance:.2f} for uniform square grid.")
                else:
                    raise ValueError(f"Mean_distance {original_distance:.2f} needs to be changed to {args.mean_distance-0.5:.2f} to have a perfect square grid.")
        base_positions = generate_square_grid(args.image_size, args.mean_distance)
    else:
        raise ValueError("grid_type must be 'hex' or 'square'")

    positions = apply_jitter(base_positions, args.image_size, args.mean_distance, args.cv_distance)
    img = np.zeros(args.image_size, dtype=np.uint8)
    radii = [max(2, np.random.normal(args.mean_radius, args.cv_radius * args.mean_radius)) for _ in positions]

    # Draw disks with periodic wrap
    shifts = [
        (0, 0), (args.image_size[0], 0),(0, args.image_size[1]), (args.image_size[0], args.image_size[1]), (-args.image_size[0], -args.image_size[1]),
        (args.image_size[0], -args.image_size[1]), (-args.image_size[0], args.image_size[1])
    ]
    ids = np.random.randint(0, high=10000, size=len(positions))
    for i, (center, radius) in enumerate(zip(positions, radii)):
        x, y = center
        for dx, dy in shifts:
            draw_center = (int(x + dx), int(y + dy))
            draw_disk(
                img, draw_center, radius, args.sine_amp1, args.sine_wave1, args.sine_amp2,
                args.sine_wave2, args.rotate_randomly, int(ids[i]), args.global_area_normalization_factor,
                global_rotation_offset = args.global_rotation_offset
            )

    stripe_metadata = []
    if args.draw_stripes:
        stripe_pairs, stripe_metadata = build_parallel_stripes(
            positions, args.image_size, args.stripe_angle_deg, args.stripe_mean_length, args.stripe_std_length
        )
        for p1, p2 in stripe_pairs:
            if args.sin_stripe:
                draw_sinusoidal_stripe(
                    img, p1, p2, args.mean_radius,
                    amplitude=args.sin_stripe_amp,
                    wavelength=args.sin_stripe_wavelength
                )
            else:
                draw_stripe(img, p1, p2, args.mean_radius, args.rotate_randomly)

    benchmark_cover = -1 # Set default to 0, this will ensure that main program will not use the benchmark cover (instead, the fraction of white pixels will be used)
    if args.enforce_area_constancy:
        if args.cur_image_fraction_pixels is None and args.global_area_normalization_factor is None:
            args.cur_image_fraction_pixels = fraction_white_pixels(img)

            # Reset sine parameters to generate a circular pattern (i.e., the same pattern as 'cur_image', but without sine waves)
            sine_amp1 = args.sine_amp1
            sine_amp2 = args.sine_amp2
            args.sine_amp1 = 0
            args.sine_amp2 = 0

            # Generate the circular pattern
            img, _, _, _, _ = create_image(**vars(args))
            args.circular_image_fraction_pixels = fraction_white_pixels(img)
            benchmark_cover = args.circular_image_fraction_pixels

            # Reset sine parameters to original values
            args.sine_amp1 = sine_amp1
            args.sine_amp2 = sine_amp2

            iterative_correction_factor = 0
            args.global_area_normalization_factor = args.cur_image_fraction_pixels / args.circular_image_fraction_pixels
            error = 100000
            stepsize = 20
            best_version = [100000, iterative_correction_factor, img, positions, radii, stripe_metadata, benchmark_cover]
            idx = 0
            max_iters = 100
            while abs(error) > 0.0001 and idx < max_iters:
                # Generate a revised version of the pattern with the area normalization factor
                img, positions, radii, stripe_metadata, _ = create_image(**vars(args))
                args.cur_image_fraction_pixels = fraction_white_pixels(img)

                # Update error and iterative correction factor
                error = 1 - args.cur_image_fraction_pixels / args.circular_image_fraction_pixels
                _stepsize = stepsize * abs(error) * abs(error) * abs(error) # Adjust step size based on error magnitude
                error_sign = error / max(0.00000000001, abs(error))
                if abs(error) < best_version[0]:
                    best_version = [abs(error), iterative_correction_factor, img, positions, radii, stripe_metadata, benchmark_cover]
                iterative_correction_factor += error_sign * _stepsize
                iterative_correction_factor = max(-0.9, iterative_correction_factor)
                idx += 1
                if (idx+1) % 10 == 0:
                    print(f"Optimizing area correction iteratively... iteration {idx+1} (maximum is {max_iters}).")

                # Derive the next global area normalization factor
                args.global_area_normalization_factor /= 1 + iterative_correction_factor  # Adjust based on sine amplitude (heuristic, needed because of artifacts in image generation)
            
            best_version.pop(1)
            print("Relative deviation from target area ratio:", best_version.pop(0))
            print("Fraction of white pixels in generated image:", fraction_white_pixels(best_version[0]))
            
            return best_version

    return img, positions, radii, stripe_metadata, benchmark_cover

def periodic_distance(p1, p2, box):
    delta = np.abs(np.array(p1) - np.array(p2))
    delta = np.where(delta > 0.5 * np.array(box), np.array(box) - delta, delta)
    return np.sqrt((delta ** 2).sum())

def plot_periodic_neighbor_distances(positions, image_size, mean_distance, k=6):
    """
    Plot disk positions and histogram of periodic distances to k nearest neighbors.
    """
    box = np.array(image_size)
    positions = np.array(positions)

    # Tile image in 3x3 to account for wraparound
    tiles = []
    shifts = [-box, [0, -box[1]], [box[0], -box[1]],
              [-box[0], 0], [0, 0], [box[0], 0],
              [-box[0], box[1]], [0, box[1]], [box]]
    for shift in shifts:
        tiles.append(positions + shift)
    tiled_positions = np.vstack(tiles)

    tree = cKDTree(tiled_positions)
    distances, indices = tree.query(positions, k=k + 1)  # k+1 to skip self (0 distance)

    all_distances = distances[:, 1:].flatten()

    # Scatter plot with connections
    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.scatter(positions[:, 0], positions[:, 1], s=10, label='Disks')
    for i in range(len(positions)):
        for j in range(1, k + 1):
            p1 = positions[i]
            p2 = tiled_positions[indices[i, j]]
            delta = p2 - p1
            # Apply periodic correction
            delta = np.where(np.abs(delta) > box / 2, -np.sign(delta) * (box - np.abs(delta)), delta)
            p2_corrected = p1 + delta
            plt.plot([p1[0], p2_corrected[0]], [p1[1], p2_corrected[1]], color='gray', alpha=0.3)
    plt.title("Disk Positions and Nearest Neighbors")
    plt.axis("equal")
    plt.xlim(0, image_size[0])
    plt.ylim(0, image_size[1])

    # Histogram
    plt.subplot(1, 2, 2)
    plt.hist(all_distances, bins=30, color='skyblue', edgecolor='black')
    plt.axvline(mean_distance, color='red', linestyle='--', label=f'Mean Target = {mean_distance}')
    plt.title("Histogram of Neighbor Distances (Periodic)")
    plt.xlabel("Distance")
    plt.ylabel("Count")
    plt.legend()

    plt.tight_layout()
    plt.show()



# Example usage
if __name__ == "__main__":
    params = {
        "image_size": (1000, 1000),
        "mean_radius": 50,
        "cv_radius": 0,
        "mean_distance": 200,
        "cv_distance": 0,
        "sine_amp1": 35,
        "sine_wave1": 6,
        "sine_amp2": 0,
        "sine_wave2": 0,
        "draw_stripes": False,
        "stripe_angle_deg": 30,
        "stripe_mean_length": 200,
        "stripe_std_length": 15,
        "sin_stripe": False,
        "sin_stripe_amp": 20,
        "sin_stripe_wavelength": 80,
        "enforce_distance_uniformity": True,
        "enforce_area_constancy": True,
        "grid_type": "square",
        "rotate_randomly": True,
        "cur_image_fraction_pixels": None,
        "circular_image_fraction_pixels": None,
        "global_area_normalization_factor": None,
        "global_rotation_offset": None
    }

    img, positions, radii, stripe_metadata, benchmark_cover = create_image(**params)
    export_metadata(positions, radii, params, stripe_metadata)
    cv2.imwrite(os.path.join(cfg.CPG_OUTPUT_DIR, "generated_pattern_sine80.png"), img)
    cv2.imshow("Generated pattern", img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    #plot_periodic_neighbor_distances(positions, params["image_size"], params["mean_distance"])

