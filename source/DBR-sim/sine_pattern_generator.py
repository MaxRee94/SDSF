import cv2
import numpy as np


def generate(dimensions, sine_amplitude=1, sine_wavelength=100,
             sine_type="horizontal", sine_offset=0, show=False, 
             maximum=None, minimum=None, mean=None, **cfg):
    """
    Generate a grayscale image with sine-based patterns.

    Parameters
    ----------
    dimensions : tuple (height, width)
        Output image size.
    sine_amplitude : float
        Amplitude of the sine wave (0–255 scale).
    sine_wavelength : float
        Wavelength in pixels.
    sine_type : str
        "horizontal", "vertical", "diagonal", or "radial".
    sine_offset : float or tuple
        Offset definition depends on sine_type:
          - horizontal: y-location of the peak from top
          - vertical: x-location of the peak from left
          - diagonal: shift along the diagonal axis
          - radial: (x, y) center of radial pattern
    maximum: float
        Maximum pixel value (minimum will be maximum - 2 * amplitude, mean will be maximum - amplitude)
    minimum: float
        Minimum pixel value (maximum will be minimum + 2 * amplitude, mean will be minimum + amplitude)
    mean: float
        Mean pixel value (maximum will be mean + amplitude, minimum will be mean - amplitude)
    show : bool
        Whether to display the image with cv2.imshow.
    """

    H, W = dimensions
    y = np.arange(H)
    x = np.arange(W)
    X, Y = np.meshgrid(x, y)

    if (maximum is not None) + (minimum is not None) + (mean is not None) > 1:
        raise ValueError("Please specify only one of these: maximum, minimum, or mean.")

    if maximum is not None:
        minimum = maximum - 2 * sine_amplitude
        print("setting minimum to", minimum)
    elif minimum is not None:
        maximum = minimum + 2 * sine_amplitude
    elif mean is not None:
        minimum = mean - sine_amplitude
        maximum = mean + sine_amplitude
    else:
        raise ValueError("Please specify exactly one of these: maximum, minimum, or mean.")

    if sine_offset is None and sine_type in ["horizontal", "vertical"]:
        raise ValueError("sine_offset cannot be None for horizontal or vertical sine types.")

    if sine_type == "horizontal":
        # Peak at sine_offset ⇒ phase shift = -2π * offset / wavelength
        phase = -2 * np.pi * (sine_offset / sine_wavelength - 0.25)
        pattern = np.sin(2 * np.pi * (Y / sine_wavelength) + phase)

    elif sine_type == "vertical":
        phase = -2 * np.pi * (sine_offset / sine_wavelength - 0.25)
        pattern = np.sin(2 * np.pi * (X / sine_wavelength) + phase)

    elif sine_type == "radial":
        if not isinstance(sine_offset, (tuple, list)):
            cx, cy = W / 2, H / 2
        else:
            cx, cy = sine_offset
        R = np.sqrt((X - cx)**2 + (Y - cy)**2)

        # Phase shift so maximum amplitude occurs at center (R=0)
        pattern = minimum + np.sin(2 * np.pi * (R / sine_wavelength) + np.pi / 2)

    else:
        raise ValueError(f"Unknown sine_type: {sine_type}")

    # Rescale value range from [-1, 1] to [minimum, maximum]
    pattern = (pattern + 1) / 2 * (maximum - minimum) + minimum

    if show:
        # Scale to grayscale image
        img_rescaled = (pattern + sine_amplitude) * 255
        img_rescaled = np.clip(img_rescaled, 0, 255).astype(np.uint8)

        cv2.imshow("Sine Pattern", img_rescaled)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    return pattern


if __name__ == "__main__":
    # Example test run
    dimensions = (240, 240)
    sine_amplitude = 127
    sine_wavelength = 300
    sine_type = "radial"
    sine_offset = 0
    show = True

    generate(dimensions, sine_amplitude, sine_wavelength,
                          sine_type=sine_type,
                          sine_offset=sine_offset,
                          show=show)


