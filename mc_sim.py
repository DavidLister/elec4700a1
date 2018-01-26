# mc_sim.py
#
# Monte Carlo simulation of electrons for ELEC 4700 assignment 1
# David Lister
# January 2018
#

#import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
import random
import numpy as np

# Todo List
#   - Iniitalize electrons
#   - Finish electron object
#   - Does the mean free path include reflections?
#   - Implement heat in different regions

# Symbolic Constants
REFLECT = "REFLECT"
TELEPORT = "TELEPORT"
DELTA = 1E-3
# Useful Lambda Functions
get_radius = lambda x, y: np.sqrt(x**2 + y**2)

# Classes

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, other):
        return self.x == other.x and self.y == other.y

    def as_tuple(self):
        return (self.x, self.y)

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
        self.velocity = get_thermal_vector(self.position)
        self.path_history = [self.position.as_tuple()]
        self.time = 0

    def plot(self, display = False):
        points = np.array(self.path_history).transpose()
        plt.plot(points[0], points[1])


class World:
    def __init__(self, border_points, teleport_edges = []):
        self.point_list = border_points
        self.border = make_border(self.point_list)
        self.bounds = get_bounds(self.point_list)
        self.teleport_edges = []

    def plot(self, display = False):
        plot_world(self.border, "World")
        if display:
            plt.show()

class Velocity:
    def __init__(self, vx, vy):
        self.x = vx
        self.y = vy
        self.r = get_radius(vx, vy)

    def as_numpy(self):
        return np.array((self.x, self.y))

    def __str__(self):
        return '(' + str(self.x) + ', ' + str(self.y) + ", r: " + str(self.r) + ')'

# Functions

def get_initial_position(world):
    xmin, xmax, ymin, ymax = world.bounds
    remote_point_top = Point(0, ymax * 2)
    remote_point_bot = Point(0, ymin * 2)
    over = False
    while not over:
        x = random.uniform(xmin, xmax)
        y = random.uniform(ymin, ymax)
        test_point = Point(x, y)
        count_top = 0
        count_bot = 0
        remote_line_top = Edge(remote_point_top, test_point)
        remote_line_bot = Edge(remote_point_bot, test_point)
        for edge in world.border:
            if get_intersection(remote_line_top, edge) is not None:
                count_top +=1
            if get_intersection(remote_line_bot, edge) is not None:
                count_bot += 1
        if count_top % 2 == 1 and count_bot % 2 == 1:
            over = True
            #print("out of bounds")
    return test_point

def get_thermal_vector(position):
    # todo - make these parameters physical
    # todo - make velocity depend on location
    vx = random.gauss(0, 10)
    vy = random.gauss(0, 10)
    return Velocity(vx, vy)

def get_intersection(e1, e2):
    #print(e1, e2)
    det = e1.a * e2.b - e2.a * e1.b
    if det == 0:
        #print("None")
        return False
    x = (e2.b * e1.c - e1.b * e2.c)/det
    y = (e1.a * e2.c - e2.a * e1.c)/det
    p = Point(x, y)
    if x >= min([e1.p1.x, e1.p2.x]) - DELTA and x <= max([e1.p1.x, e1.p2.x]) + DELTA and y >= min([e1.p1.y, e1.p2.y]) - DELTA and y <= max([e1.p1.y, e1.p2.y]) + DELTA and x >= min([e2.p1.x, e2.p2.x]) - DELTA and x <= max([e2.p1.x, e2.p2.x]) + DELTA and y >= min([e2.p1.y, e2.p2.y]) - DELTA and y <= max([e2.p1.y, e2.p2.y]) + DELTA:
        #print("in", p)
        return p
    #print("out", p)
    return None


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
    #plt.savefig(name + ".png")


def get_distance(mfp):
    # Returns a random distance given a mean free path
    val = random.random()
    if val == 0:
        val = 1e-4
    val = val + 0j
    return mfp * np.real(np.log(-1/val))

def generate_electron_list(world, n):
    lst = []
    for i in range(n):
        if i%1000 == 0:
            print(i)
        lst.append(Electron(world))

    return lst

def plot_electron_list(electron_list, display = False):
    for electron in electron_list:
        electron.plot()
    if display:
        plt.show()

def get_position_from_velocity(point, velocity, distance):
    uvector = velocity.as_numpy()
    uvector = uvector/np.linalg.norm(uvector)
    path = uvector * distance
    return Point(point.x + path[0], point.y + path[1])

def get_intersections_with_border(border, line):
    points = []
    for edge in border:
        if get_intersection(edge, line) is not None:
            points.append(get_intersection(edge, line))
    return points

def update_electron_scatter(electron, new_position):
    electron.position = new_position
    electron.path_history.append(new_position.as_tuple())
    electron.velocity = get_thermal_vector(new_position)
    electron.time += 1 #todo - fix this

def update_electron_reflect(electron, point, reflection):
    electron.position = point
    electron.position.x += reflection[0] * DELTA
    electron.position.y += reflection[1] * DELTA
    electron.velocity.x = reflection[0]
    electron.velocity.y = reflection[1]
    electron.path_history.append(point.as_tuple())
    electron.time += 1 # todo - fix this

def get_closest_point(ref, lst):
    dx = lst[0].x - ref.x
    dy = lst[0].y - ref.y
    minr = get_radius(dx, dy)
    minpt = lst[0]
    for point in lst:
        dx = ref.x - point.x
        dy = ref.y - point.y
        if get_radius(dx, dy) < minr:
            minr = get_radius(dx, dy)
            minpt = point
    return minpt

def is_point_on_edge(point, edge):
    return abs(edge.c - (edge.a * point.x + edge.b * point.y)) < DELTA

def get_reflect_info(world, point):
    for edge in world.teleport_edges:
        if is_point_on_edge(point, edge):
            return TELEPORT, edge

    for edge in world.border:
        if is_point_on_edge(point, edge):
            return REFLECT, edge
    print("ERROR -Point not on edge")

def simulate_electron_(world, electron, end_time, mean_free_path):
    #print(electron.position)
    over = False
    while not over:
        if electron.time >= end_time:
            over = True
        distance = get_distance(mean_free_path)
        test_position = get_position_from_velocity(electron.position, electron.velocity, distance)
        test_line = Edge(electron.position, test_position)
        intersections = get_intersections_with_border(world.border, test_line)
        if len(intersections) == 0:
            #print("Scatter", str(test_position))
            update_electron_scatter(electron, test_position)
        else:
            #for point in intersections:
                #print(point)
            reflect_point = get_closest_point(electron.position, intersections)
            reflect_type, reflect_edge = get_reflect_info(world, reflect_point)
            #todo - implement teleport border

            normal = np.array(((reflect_edge.p2.y - reflect_edge.p1.y), (reflect_edge.p1.x - reflect_edge.p2.x)))
            normal = normal / np.linalg.norm(normal)
            incident = np.array((electron.velocity.x, electron.velocity.y))
            reflection = incident - 2 * np.dot(incident, normal) * normal
            #print('Reflect', str(reflect_point), electron.velocity, reflection, get_radius(reflection[0], reflection[1]))
            update_electron_reflect(electron, reflect_point, reflection)


# Main Simulation

if __name__ == "__main__":
    pt_list = [(-200, -100), (-20, -100), (-20, -60), (-10, -50), (10, -50), (20, -60), (20, -100), (200, -100), (200, 100), (20, 100), (20, 60), (10, 50), (-10, 50), (-20, 60), (-20, 100), (-200, 100)]
    world = World(pt_list)
    electrons = generate_electron_list(world, 10)
    for electron in electrons:
        simulate_electron_(world, electron, 100, 150)
    world.plot()
    plot_electron_list(electrons, True)
    #test_point = Point(12.4805, -64.5816)
    #xmin, xmax, ymin, ymax = world.bounds
    #remote_point = Point(0, ymax * 2)
    #remote_line = Edge(remote_point, test_point)
    #count = 0
    #for edge in world.border:
    #    if get_intersection(remote_line, edge) is not None:
    #        count += 1
    #print(count)
