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


OUTPUT_DIR = r"F:\Development\DBR-sim\data_out\pattern_metadata" # TODO: Modify to use default defined in config.py.

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

def generate_disk(center, base_radius, amp1=0, wave1=1, amp2=0, wave2=2, resolution=360):
    angles = np.linspace(0, 2 * np.pi, resolution, endpoint=False)
    if wave2 <= wave1 and amp2 != 0:
        raise ValueError("Second sine wave wavelength must be shorter than the first.")
    r = base_radius + amp1 * np.sin(wave1 * angles + random.uniform(0, 2*math.pi)) + amp2 * np.sin(wave2 * angles + random.uniform(0, 2*math.pi))
    x = center[0] + r * np.cos(angles)
    y = center[1] + r * np.sin(angles)
    return np.stack((x, y), axis=-1).astype(np.int32)

def draw_disk(img, center, radius, amp1, wave1, amp2, wave2):
    contour = generate_disk(center, radius, amp1, wave1, amp2, wave2)
    cv2.fillPoly(img, [contour], 255)

def draw_stripe(img, center1, center2, radius):
    disk1 = generate_disk(center1, radius)
    disk2 = generate_disk(center2, radius)
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

def export_metadata(positions, radii, parameters, stripe_metadata=None, output_dir=OUTPUT_DIR):
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
        df_stripes.to_csv(os.path.join("output_dir", "stripe_metadata.csv"), index=False)

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


def create_image(
    image_size=(1000, 1000),
    mean_radius=30,
    cv_radius=5,
    mean_distance=200,
    cv_distance=10,
    sine_amp1=5,
    sine_wave1=6,
    sine_amp2=2,
    sine_wave2=12,
    draw_stripes=True,
    stripe_angle_deg=0,
    stripe_mean_length=200,
    enforce_distance_uniformity=True,
    stripe_std_length=20,
    sin_stripe=False,
    sin_stripe_amp=10,
    sin_stripe_wavelength=100,
    grid_type="hex"
):
    
    if grid_type == "hex":
        if enforce_distance_uniformity:
            original_distance = mean_distance
            mean_distance = adjust_mean_distance_for_uniform_hex_grid(image_size, mean_distance)
            if not math.isclose(original_distance, mean_distance, abs_tol=1):
                raise RuntimeError(f"Mean_distance {original_distance:.2f} needs to be changed to {mean_distance-0.5} to have a perfect hex grid.")
        base_positions = generate_hex_grid_with_filter(image_size, mean_distance)
    elif grid_type == "square":
        if enforce_distance_uniformity:
            original_distance = mean_distance
            mean_distance = adjust_mean_distance_for_uniform_square_grid(image_size, mean_distance)
            if not math.isclose(original_distance, mean_distance, abs_tol=1):
                raise RuntimeError(f"Mean_distance {original_distance:.2f} needs to be changed to {mean_distance:.2f} to have a perfect square grid.")
        base_positions = generate_square_grid(image_size, mean_distance)
    else:
        raise ValueError("grid_type must be 'hex' or 'square'")

    positions = apply_jitter(base_positions, image_size, mean_distance, cv_distance)
    img = np.zeros(image_size, dtype=np.uint8)
    radii = [max(2, np.random.normal(mean_radius, cv_radius * mean_radius)) for _ in positions]

    # Draw disks with periodic wrap
    for i, (center, radius) in enumerate(zip(positions, radii)):
        x, y = center
        shifts = [(0, 0), (image_size[0], 0), (-image_size[0], 0),
                  (0, image_size[1]), (0, -image_size[1]),
                  (image_size[0], image_size[1]), (-image_size[0], -image_size[1]),
                  (image_size[0], -image_size[1]), (-image_size[0], image_size[1])]
        for dx, dy in shifts:
            draw_center = (int(x + dx), int(y + dy))
            draw_disk(img, draw_center, radius, sine_amp1, sine_wave1, sine_amp2, sine_wave2)

    stripe_metadata = []
    if draw_stripes:
        stripe_pairs, stripe_metadata = build_parallel_stripes(
            positions, image_size, stripe_angle_deg, stripe_mean_length, stripe_std_length
        )
        for p1, p2 in stripe_pairs:
            if sin_stripe:
                draw_sinusoidal_stripe(
                    img, p1, p2, mean_radius,
                    amplitude=sin_stripe_amp,
                    wavelength=sin_stripe_wavelength
                )
            else:
                draw_stripe(img, p1, p2, mean_radius)


    return img, positions, radii, stripe_metadata

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
        "mean_radius": 80,
        "cv_radius": 0,
        "mean_distance": 288,
        "cv_distance": 0,
        "sine_amp1": 10,
        "sine_wave1": 6,
        "sine_amp2": 0,
        "sine_wave2": 12,
        "draw_stripes": False,
        "stripe_angle_deg": 30,
        "stripe_mean_length": 200,
        "stripe_std_length": 15,
        "sin_stripe": True,
        "sin_stripe_amp": 15,
        "sin_stripe_wavelength": 80,
        "enforce_distance_uniformity": True,
        "grid_type": "hex"
    }

    img, positions, radii, stripe_metadata = create_image(**params)
    export_metadata(positions, radii, params, stripe_metadata)
    cv2.imshow("Generated pattern", img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    plot_periodic_neighbor_distances(positions, params["image_size"], params["mean_distance"])

