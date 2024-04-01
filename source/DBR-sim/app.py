import sys

import cv2
import visualization as vis
import numpy as np
import time

sys.path.append(r"F:\Development\DBR-sim\build")
from x64.Debug import dbr_cpp as cpp
#from x64.Release import dbr_cpp as cpp



def init(timestep=None, gridsize=None, cellsize=None, max_radius=None,
         treecover=None, self_ignition_factor=None, flammability=None,
         rainfall=None, unsuppressed_flammability=None, suppressed_flammability=None,
         verbosity=None, radius_q1=None, radius_q2=None,
         **_):
    cpp.check_communication()
    dynamics = cpp.Dynamics(
        timestep, unsuppressed_flammability, suppressed_flammability,
        self_ignition_factor, rainfall, verbosity
    )
    dynamics.init_state(gridsize, cellsize, max_radius, radius_q1, radius_q2)
    dynamics.state.set_tree_cover(treecover)
    
    return dynamics


def termination_condition_satisfied(dynamics, start_time, user_args):
    satisfied = False
    condition = ""
    if (time.time() - start_time > user_args["timelimit"]):
        condition = f"Exceeded timelimit ({user_args['timelimit']} seconds)."
    if (dynamics.state.population.size() == 0):
        condition = f"Population collapse."
    #print("Pop size: ", dynamics.state.population.size())
    satisfied = len(condition) > 0
    if satisfied:
        print("Simulation terminated. Cause:", condition)

    return satisfied


def updateloop(dynamics, **user_args):
    start = time.time()
    is_within_timelimit = True;
    while not termination_condition_satisfied(dynamics, start, user_args):
        dynamics.save_2_trees()
        vis.visualize(
            dynamics.state.grid, user_args["image_width"], collect_states=False,
            color_dict={
                0: np.array((0, 0, 0), np.uint8),
                1: np.array((0, 255, 0), np.uint8),
                2: np.array((0, 0, 255), np.uint8),
                3: np.array((255, 0, 255), np.uint8),
                4: np.array((255, 0, 0), np.uint8),
                5: np.array((0, 100, 200), np.uint8),
                6: np.array((255, 255, 255), np.uint8),
                7: np.array((100, 0, 200), np.uint8),
                8: np.array((10, 30, 255), np.uint8),
                9: np.array((100, 200, 80), np.uint8),
            }
        )
        time.sleep(4)
        dynamics.update()
        vis.visualize(
            dynamics.state.grid, user_args["image_width"], collect_states=False,
            color_dict={
                0: np.array((0, 0, 0), np.uint8),
                1: np.array((0, 255, 0), np.uint8),
                2: np.array((0, 0, 255), np.uint8),
                3: np.array((255, 0, 255), np.uint8),
                4: np.array((255, 0, 0), np.uint8),
                5: np.array((0, 100, 200), np.uint8),
                6: np.array((255, 255, 255), np.uint8),
                7: np.array((100, 0, 200), np.uint8),
                8: np.array((10, 30, 255), np.uint8),
                9: np.array((100, 200, 80), np.uint8),
            }
        )
        time.sleep(4)
        # print("Repopulating...")
        # dynamics.state.repopulate_grid(0)
        #vis.visualize(dynamics.state.grid, user_args["image_width"])
        # time.sleep(4)

    cv2.destroyAllWindows()


def do_tests(**user_args):
    tests = cpp.Tests()
    tests.run_all(True)


def main(**user_args):
    if user_args["test"] == "all":
        do_tests(**user_args)
        return
        
    dynamics = init(**user_args)
    updateloop(dynamics, **user_args)
 


    

