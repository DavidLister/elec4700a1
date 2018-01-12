# main.py Monte Carlo Modeling of Electron Transport
# For ELEC 4700
# January 2018
# David Lister
#

# Assumptions:
# Distance units are in nanometers
# Time units are in nanoseconds
# Mass is in kg
# The point (0, 0) is playable

import random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REFLECT = "reflect"
CONTINUOUS = "continuous"

ID, X, Y, VX, VY, NEXT_EVENT_TIME, NEXT_EVENT_OBJECT = [x for x in range(7)]

class Region:
    def __init__(self, playable, c1, c2, bc):
        """
        bool playable - Note, the order of stacking will determine overlapping values, the last value will be used.
        (x, y) c1 - corner 1
        (x, y) c2 - corner 2
        (Vertical, Horizontal) bc - Boundary Conditions
        """
        self.playable = playable
        self.min_x = min([c1[0], c2[0]])
        self.min_y = min([c1[1], c2[1]])
        self.max_x = max([c1[0], c2[0]])
        self.max_y = max([c1[1], c2[1]])
        self.bcVertical = bc[0]
        self.bcHorizontal = bc[1]

def get_bounds(world):
    """Returns the bounds of the playable world to get a range of values for generating initial electrons.
        This assumes the origin is within the playable region."""
    min_x = 0
    min_y = 0
    max_x = 0
    max_y = 0
    for obj in world:
        if obj.playable:
            if obj.min_x < min_x:
                min_x = obj.min_x
            if obj.max_x > max_x:
                max_x = obj.max_x
            if obj.min_y < min_y:
                min_x = obj.min_y
            if obj.max_y > max_y:
                max_x = obj.max_y
    return (min_x, max_x, min_y, max_y)

def is_point_in_object(point, object):
    if point[0] > object.min_x and point[0] < object.max_x and point[1] > object.min_y and point[1] < object.min_y:
        return True

def get_speed():
    r = random.gammavariate(2, 2)
    theta = random.uniform(0, 2 * np.pi)
    vx = r * np.cos(theta)
    vy = r * np.sin(theta)
    return vx, vy

def generate_electron(id, world, bounds):
    # Generate position
    over = False
    while not over:
        x = random.uniform(bounds[0], bounds[1])
        y = random.uniform(bounds[2], bounds[3])
        for obj in world:
            if is_point_in_object((x, y), obj):
                over = obj.playable

    vx, vy = get_speed()
    return [id, x, y, vx, vy, -1, -1]

def get_edges_from_object(obj):
    l1 = ((obj.min_x, obj.max_x), (obj.min_y, obj.min_y))
    l2 = ((obj.max_x, obj.max_x), (obj.max_y, obj.min_y))
    l3 = ((obj.min_x, obj.max_x), (obj.max_y, obj.max_y))
    l4 = ((obj.min_x, obj.min_x), (obj.max_y, obj.min_y))
    return [l1, l2, l3, l4]

def define_border(world):
    pass




def set_next_event(electon):
    pass



if __name__ == "__main__":
    world = []
    world.append(Region(playable=True,
                        c1=(-100, -50),
                        c2=(100, 50),
                        bc=(REFLECT, CONTINUOUS)
                        ))

    bounds = get_bounds(world)
    state = []

    edges = get_edges_from_object(world[0])
    for e in edges:
        plt.plot(e[0], e[1])
    plt.savefig("test.png")
