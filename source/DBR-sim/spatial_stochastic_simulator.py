import gstools as gs
import cv2
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace

class SpatialStochasticSimulator:
    def __init__(self, args: SimpleNamespace):
        """
        Initialize simulator using parameters from args.
        
        Expects args to have:
            - width       : int       → grid width (nx = ny = width)
            - model_type  : str       → "Exponential", "Gaussian", "Spherical"
            - sill        : float
            - nugget      : float
            - vrange      : float
            - global_seed : int       → random seed
        """
        self.nx = self.ny = args.width
        self.model_type = args.model_type
        self.sill = args.sill
        self.nugget = args.nugget
        self.vrange = args.vrange
        self.seed = args.global_seed

    def generate(
        self,
        output_file="variogram_simulation.png",
        plot_variogram=False,
        plot_image=False
    ):
        """
        Generate a 2D random field and save as grayscale image.
        
        Optional plotting:
            plot_variogram : bool → display semivariogram
            plot_image     : bool → display generated image
        """
        # Define covariance model and random field
        model = getattr(gs, self.model_type)(
            dim=2, var=self.sill, len_scale=self.vrange, nugget=self.nugget
        )
        srf = gs.SRF(model, seed=self.seed)
        field = srf((self.nx, self.ny))

        # Normalize to 0–255
        img = ((field - field.min()) / (field.max() - field.min()) * 255).astype(np.uint8)
        cv2.imwrite(output_file, img)
        print(f"Simulation complete: {output_file} written")

        # Optional plotting
        if plot_variogram:
            lags, gamma = model.semivariance(np.linspace(0, self.vrange*2, 100))
            plt.figure(figsize=(5,3))
            plt.plot(lags, gamma, label=f"{self.model_type} variogram")
            plt.xlabel("Lag distance")
            plt.ylabel("Semivariance")
            plt.title("Variogram")
            plt.grid(True)
            plt.legend()
            plt.show()

        if plot_image:
            print(img.shape)
            cv2.imshow("Simulated Random Field", img)
            cv2.waitKey(0)
            cv2.destroyAllWindows()


# ==========================
# Example usage
# ==========================
if __name__ == "__main__":
    # Create an args object with required parameters
    args = SimpleNamespace(
        width=512,
        model_type="Spherical",
        sill=1.0,
        nugget=0.1,
        vrange=40.0,
        global_seed=42
    )

    # Instantiate simulator
    simulator = SpatialStochasticSimulator(args)

    # Generate simulation with plotting enabled
    simulator.generate(plot_variogram=False, plot_image=True)
