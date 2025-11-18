import cv2
import numpy as np


def generate(amplitude=0.3, scale=1, show=False, offset=0, grid_width=1000, **args):
    """
    Generate a 1000x1000 uniform macropixel noise pattern.

    Each macropixel has a uniform random value in [0.5-amplitude, 0.5+amplitude],
    clipped to [0,1], and the macropixel size in the *final* image is exactly
    <scale> pixels.

    Parameters
    ----------
    amplitude : float
        Noise half-range. Values sampled from [0.5 - amplitude, 0.5 + amplitude].
    scale : int
        Final macropixel width (and height) in pixels.
    show : bool
        Whether to show the final 1000x1000 image.
    """

    FINAL_SIZE = grid_width

    # Determine number of macropixels in each direction
    macro_h = FINAL_SIZE // scale
    macro_w = FINAL_SIZE // scale

    # Generate macropixel noise grid
    noise = np.random.uniform(
        low=0.5 - amplitude,
        high=0.5 + amplitude,
        size=(macro_h, macro_w)
    ) + offset

    # Clip to [0,1]
    noise = np.clip(noise, 0.0, 1.0)

    # Upscale macropixel grid to 1000x1000
    img = cv2.resize(
        noise,
        (FINAL_SIZE, FINAL_SIZE),
        interpolation=cv2.INTER_CUBIC
    )

    # Convert to uint8 grayscale image
    img_uint8 = (img * 255).astype(np.uint8)

    if show:
        cv2.imshow("Uniform Macropixel Noise", img_uint8)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    return img_uint8


if __name__ == "__main__":
    # Example usage
    amplitude = 0.3   # uniform noise range is [0.2, 0.8]
    scale = 1        # each macropixel is 50x50 pixels in the final image
    show = True
    grid_width = 1000
    offset = -0.5

    generate(amplitude=amplitude, scale=scale, show=show, grid_width=grid_width, offset=offset)
