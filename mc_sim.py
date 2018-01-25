# mc_sim.py
#
# Monte Carlo simulation of electrons for ELEC 4700 assignment 1
# David Lister
# January 2018
#

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import random
import numpy as np

# Todo List
#   - Iniitalize electrons
#   - Finish electron object
#   - Does the mean free path include reflections?
#   - Implement heat in different regions

# Classes

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def __str__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + ')'

class Edge:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2
        # Form ax + by = c
        self.a = p2.y - p1.y
        self.b = p1.x - p2.x
        self.c = self.a * p1.x + self.b * p1.y

    def plot(self):
        return plt.plot((self.p1.x, self.p2.x), (self.p1.y, self.p2.y))

    def __str__(self):
        return str(self.p1) + '---' + str(self.p2)

class Electron:
    def __init__(self, world):
        self.position = get_initial_position(world)


# Functions

def get_intersection(e1, e2):
    det = e1.a * e2.b - e2.a * e1.b
    if det == 0:
        return False
    x = (e2.b * e1.c - e1.b * e2.c)/det
    y = (e1.a * e2.c - e2.a * e1.c)/det
    p = Point(x, y)
    if x >= min([e1.p1.x, e1.p2.x]) and x <= max([e1.p1.x, e1.p2.x]) and y >= min([e1.p1.y, e1.p2.y]) and y <= max([e1.p1.y, e1.p2.y]) and x >= min([e2.p1.x, e2.p2.x]) and x <= max([e2.p1.x, e2.p2.x]) and y >= min([e2.p1.y, e2.p2.y]) and y <= max([e2.p1.y, e2.p2.y]):
        return p
    return False


def make_border(point_list):
    pts = [Point(p[0], p[1]) for p in point_list]
    edges = [Edge(pts[i], pts[i+1]) for i in range(len(pts) - 1)] + [Edge(pts[-1], pts[0])]
    return edges

def get_bounds(point_list):
    xlst = [pt[0] for pt in point_list]
    ylst = [pt[1] for pt in point_list]
    return(min(xlst), max(xlst), min(ylst), max(ylst) )

def plot_world(border, name):
    for edge in border:
        edge.plot()
    plt.savefig(name + ".png")

def get_distance(mfp):
    # Returns a random distance given a mean free path
    val = random.random()
    if val == 0:
        val = 1e-4
    val = val + 0j
    return mfp * np.real(np.log(-1/val))

# Main Simulation

if __name__ == "__main__":
    pt_list = [(-200, -100), (-20, -100), (-20, -60), (-10, -50), (10, -50), (20, -60), (20, -100), (200, -100), (200, 100), (20, 100), (20, 60), (10, 50), (-10, 50), (-20, 60), (-20, 100), (-200, 100)]
    border = make_border(pt_list)
    bounds = get_bounds(pt_list)
    plot_world(border, "world")
    dists = []
    for x in range(1000000):
        dists.append(get_distance(10))
    print("mean", np.mean(dists))
