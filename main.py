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

REFLECT = "REFLECT"
CONTINUOUS = "CONTINUOUS"
PARALLEL = "PARALLEL"
IN_RANGE = "IN_RANGE"
NOT_IN_RANGE = "NOT_IN_RANGE"


class Region:
    def __init__(self, playable, p1, p2, bc):
        """
        bool playable - Note, the order of stacking will determine overlapping values, the last value will be used.
        Point Objcet p1 - point 1
        Point Object p2 - point 2
        (Vertical, Horizontal) bc - Boundary Conditions
        """
        self.playable = playable
        self.min_x = min([p1.x, p2.x])
        self.min_y = min([p1.y, p2.y])
        self.max_x = max([p1.x, p2.x])
        self.max_y = max([p1.y, p2.y])
        self.bcVertical = bc[0]
        self.bcHorizontal = bc[1]
        self.edges = get_edges_from_region(self)

    def plot(self):
        for edge in self.edges:
            edge.plot()

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __str__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + ')'

class Vector:
    def __init__(self, point, xmag, ymax):
        self.point = point
        self.x = point.x
        self.y = point.y
        self.vx = xmag
        self.vy = ymag

class Edge:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        # In the form ax + by = c
        self.a = p2.y - p1.y
        self.b = p1.x - p2.x
        self.c = self.a * p1.x + self.b * p1.y

    def __str__(self):
        return str(self.p1) + "---" + str(self.p2)

    def plot(self):
        plt.plot((self.p1.x, self.p2.x), (self.p1.y, self.p2.y))

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

def get_edges_from_region(obj):
    l1 = Edge(Point(obj.min_x, obj.min_y), Point(obj.max_x, obj.min_y))
    l2 = Edge(Point(obj.max_x, obj.max_y), Point(obj.max_x, obj.min_y))
    l3 = Edge(Point(obj.min_x, obj.max_y), Point(obj.max_x, obj.max_y))
    l4 = Edge(Point(obj.min_x, obj.max_y), Point(obj.min_x, obj.min_y))
    return [l1, l2, l3, l4]

def get_intesection_of_lines(e1, e2):
    det = e1.a * e2.b - e2.a * e1.b
    if det == 0:
        return [PARALLEL, Point(0, 0)]
    x = (e2.b * e1.c - e1.b * e2.c)/det
    y = (e1.a * e2.c - e2.a * e1.c)/det
    p = Point(x, y)
    if x >= min([e1.p1.x, e1.p2.x]) and x <= max([e1.p1.x, e1.p2.x]) and y >= min([e1.p1.y, e1.p2.y]) and y <= max([e1.p1.y, e1.p2.y]) and x >= min([e2.p1.x, e2.p2.x]) and x <= max([e2.p1.x, e2.p2.x]) and y >= min([e2.p1.y, e2.p2.y]) and y <= max([e2.p1.y, e2.p2.y]):
        return [IN_RANGE, p]
    return [NOT_IN_RANGE, p]

def get_next_intersection(line, world):
    all_edges = []
    for region in world:
        for edge in region.edges:
            all_edges.append(edge)

    intersections = []
    for edge in all_edges:
        status, point = get_intersection_of_lines(line, edge)
        if status == IN_RANGE:
            intersections.apend(point


def define_border(world):
    min_x = min([region.min_x for region in world])
    initial_edge = Edge(Point(0, 0), Point(2 * min_x, 0))
    border = []
    over = False
    while not over:
        if is_complete(border):
            over = True



def set_next_event(electon):
    pass

def is_complete(border):
    if border[0].p1 == border[-1].p2

if __name__ == "__main__":
    world = []
    world.append(Region(playable=True,
                        p1=Point(-100, -50),
                        p2=Point(100, 50),
                        bc=(REFLECT, CONTINUOUS)
                        ))

    world.append(Region(playable=True,
                        p1=Point(-50, -100),
                        p2=Point(50, 100),
                        bc=(REFLECT, CONTINUOUS)
                        ))

    bounds = get_bounds(world)
    state = []
    for region in world:
        region.plot()
    plt.savefig("test.png")

    e1 = Edge(Point(-100, 0), Point(100,0))
    e2 = Edge(Point(0, 100), Point(0,-100))
    e3 = Edge(Point(-100, 10), Point(100,10))
    e4 = Edge(Point(0, 100), Point(100,10))

    print(get_intesection_of_lines(e1, e2))
    print(get_intesection_of_lines(e1, e3))
    print(get_intesection_of_lines(e1, e4))
