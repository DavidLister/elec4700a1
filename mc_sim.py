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
from scipy.ndimage import gaussian_filter
from PIL import Image

# Todo List
#   - Reflextion boundaries
#   - Implement heat in different regions

# Symbolic Constants
REFLECT = "REFLECT"
TELEPORT = "TELEPORT"
DELTA = 1E-9
# Useful Lambda Functions
get_radius = lambda x, y: np.sqrt(x**2 + y**2)
scale_to_red_blue = lambda val, maximum: (int(255 * (val/maximum)), 0, int(255 * (1 - val/maximum)))

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
        self.teleport_position_change = (0, 0)

    def plot(self):
        return plt.plot((self.p1.x, self.p2.x), (self.p1.y, self.p2.y))

    def __str__(self):
        return str(self.p1) + '---' + str(self.p2)

class Electron:
    def __init__(self, world):
        self.position = get_initial_position(world)
        self.velocity = get_thermal_vector(self.position)
        self.path_history = [self.position.as_tuple()]
        self.time_history = []
        self.velocity_history = [self.velocity.as_tuple()]
        self.events = 0

    def plot(self):
        segments = []
        work_in_progress = []
        for point in self.path_history:
            if point == TELEPORT:
                segments.append(work_in_progress)
                work_in_progress = []

            else:
                work_in_progress.append(point)
        segments.append(work_in_progress)

        for segment in segments:
            segment = np.array(segment).transpose()
            plt.plot(segment[0], segment[1])


class World:
    def __init__(self, border_points):
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

    def as_tuple(self):
        return (self.x, self.y)

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
    vx = 1e9 * random.gauss(0, thermal_sigma)
    vy = 1e9 * random.gauss(0, thermal_sigma)
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

def get_distance_mft(mft, velocity):
    mfp = mft * velocity
    return get_distance(mfp)

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
    distance =get_path_length(electron.position.as_tuple(), new_position.as_tuple())
    velocity = electron.velocity.r
    dt = distance/velocity
    electron.time_history.append(dt)
    electron.position = new_position
    electron.path_history.append(new_position.as_tuple())
    electron.velocity = get_thermal_vector(new_position)
    electron.velocity_history.append(electron.velocity.as_tuple())
    electron.events += 1

def update_electron_reflect(electron, point, reflection):
    distance = get_path_length(electron.position.as_tuple(), point.as_tuple())
    velocity = electron.velocity.r
    dt = distance / velocity
    electron.time_history.append(dt)
    electron.position = point
    electron.position.x += (reflection[0]/velocity) * DELTA * 2
    electron.position.y += (reflection[1]/velocity) * DELTA * 2
    electron.velocity.x = reflection[0]
    electron.velocity.y = reflection[1]
    electron.path_history.append(point.as_tuple())
    electron.velocity_history.append(electron.velocity.as_tuple())
    electron.events += 1

def update_electron_teleport(electron, point, delta):
    #print("Teleporting")
    electron.path_history.append(point.as_tuple())
    #print(point)
    distance = get_path_length(electron.position.as_tuple(), point.as_tuple())
    velocity = electron.velocity.r
    dt = distance / velocity
    electron.time_history.append(dt)
    electron.position = point
    electron.position.x += delta[0] + (electron.velocity.x/velocity) * DELTA * 2.0
    electron.position.y += delta[1] + (electron.velocity.y/velocity) * DELTA * 2.0
    electron.path_history.append(TELEPORT)
    electron.path_history.append(electron.position.as_tuple())
    electron.velocity_history.append(electron.velocity.as_tuple())
    #print(electron.position)



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

def simulate_electron_(world, electron, event_count, mean_free_time):
    #print(electron.position)
    over = False
    while not over:
        if electron.events >= event_count:
            over = True
        distance = get_distance_mft(mean_free_time, electron.velocity.r)
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

            if reflect_type == TELEPORT:
                #print("Teleport!")
                update_electron_teleport(electron, reflect_point, reflect_edge.teleport_position_change)

            else:
                normal = np.array(((reflect_edge.p2.y - reflect_edge.p1.y), (reflect_edge.p1.x - reflect_edge.p2.x)))
                normal = normal / np.linalg.norm(normal)
                incident = np.array((electron.velocity.x, electron.velocity.y))
                reflection = incident - 2 * np.dot(incident, normal) * normal
                #print('Reflect', str(reflect_point), electron.velocity, reflection, get_radius(reflection[0], reflection[1]))
                update_electron_reflect(electron, reflect_point, reflection)

def run_simulation(world, electrons, iterations, mean_free_path, print_status_every = 10):
    tally = 0
    for electron in electrons:
        if tally % print_status_every == 0:
            print("Simulating " + str(tally) + " of " + str(len(electrons)) + " electrons")
        simulate_electron_(world, electron, iterations, mean_free_path)
        tally += 1

def get_path_length(t1, t2):
    dx = t1[0] - t2[0]
    dy = t1[1] - t2[1]
    return get_radius(dx, dy)

def plot_path_histogram(electrons, bin_count, show=True):
    path_list = []
    skipnext = False
    for electron in electrons:
        last = electron.path_history[0]
        for step in electron.path_history[1:]:
            if step == TELEPORT:
                skipnext = True
            elif skipnext:
                skipnext = False
            else:
                path_list.append(get_path_length(last, step))
                last = step
    n, bins, patches = plt.hist(path_list, bin_count)
    if show:
        plt.title("Path length histogram")
        plt.xlabel("Path length")
        plt.ylabel("Probability")
        plt.show()
    return np.mean(path_list)

def plot_velocity_histogram(electrons, bin_count, show=True):
    velocity_list = []
    velocity_x = []
    velocity_y = []
    for electron in electrons:
        for step in electron.velocity_history:
            velocity_list.append(get_radius(step[0], step[1]) * 1e-9)
            velocity_x.append(step[0] * 1e-9)
            velocity_y.append(step[1] * 1e-9)
    n, bins, patches = plt.hist(velocity_list, bin_count)
    if show:
        plt.title("Velocity histogram")
        plt.xlabel("Velocity (m/s)")
        plt.ylabel("Probability")
        plt.show()
    print("Mean velocity:", np.mean(velocity_list))
    print("Mean x velocity:", np.mean(velocity_x), "Mean y velocity:", np.mean(velocity_y))
    print("Average temperature", mass * np.mean(velocity_list) ** 2 / (2 * boltzmann))

def int_by_direction(value, dir):
    if np.sign(dir) > 0:
        return np.ceil(value)
    else:
        return np.floor(value)

def make_heatmap(world, electron_list, scale = 1):
    x_range = [x for x in range(int(world.bounds[0] * scale), int(world.bounds[1] * scale + 2))]
    x_index = {x_range[i]:i for i in range(len(x_range))}
    y_range = [y for y in range(int(world.bounds[2] * scale), int(world.bounds[3] * scale + 2))]
    y_index = {y_range[i]:i for i in range(len(y_range))}
    frame = [[0 for x in x_range] for y in y_range]
    frame = np.array(frame, dtype="float64")


    skipnext = False
    tally = 0
    for electron in electron_list:
        tally += 1
        if tally %10 == 0:
            print("Processing", tally, "out of", len(electron_list), "electrons")
        last = electron.path_history[0]
        index = 0
        for point in electron.path_history[1:]:
            if point == TELEPORT:
                skipnext = True
            elif skipnext:
                skipnext = False
                last = point
            else:
                if point[0] < last[0]:
                    left = point
                    right = last
                else:
                    left = last
                    right = point
                #plt.plot((left[0], right[0]), (left[1], right[1]))
                #print(index, len(electron.path_history), len(electron.velocity_history), electron.path_history.count(TELEPORT), electron.path_history.index(point))
                velocity = electron.velocity_history[index]
                velocity = get_radius(velocity[0]*scale, velocity[1]*scale)
                dx = right[0]*scale - left[0]*scale
                dy = right[1]*scale - left[1]*scale
                slope = dy/dx
                over = False
                current_x = left[0] * scale
                current_y = left[1] * scale
                while not over:
                    delta_x = np.ceil(current_x) - current_x
                    delta_y = int_by_direction(current_y, slope) - current_y
                    r_squared_x = delta_x**2 +  (delta_x * slope)**2
                    r_squared_y = delta_y**2 + (delta_y/slope)**2
                    r_squared_to_end = (right[0]*scale - current_x)**2 + (right[1]*scale - current_y)**2
                    min_r = np.min([r_squared_x, r_squared_y, r_squared_to_end])
                    time = np.sqrt(min_r)/velocity
                    #print(slope, left, right, current_x, current_y)
                    frame[y_index[np.ceil(current_y)]][x_index[np.ceil(current_x)]] += time
                    if min_r == r_squared_x:
                        current_x = np.ceil(current_x) + DELTA
                        current_y = current_y + slope * (delta_x + DELTA)
                    elif min_r == r_squared_y:
                        current_x = current_x + ((delta_y + DELTA)/slope)
                        current_y = int_by_direction(current_y, slope) + DELTA * np.sign(slope)
                    else:
                        over = True
                last = point
                index += 1
    #plt.show()
    mask = np.ma.masked_equal(frame, 0)
    frame = gaussian_filter(frame, sigma=scale)
    maximum = np.max(frame)

    img = Image.new("RGB", (len(frame[0]), len(frame)))
    for row in range(len(frame)):
        for col in range(len(frame[0])):
            if mask[row][col] != 0:
                img.putpixel((col, row), scale_to_red_blue(frame[row][col], maximum))

            else:
                img.putpixel((col, row), (50, 50, 50))

    img.save("heatmap.png")
    plt.imshow(frame)
    plt.show()





# Main Simulation

if __name__ == "__main__":
    temperature = 300 # K
    mass = 0.26*9.10938356e-31 #kg
    boltzmann = 1.38064852e-23 # m^2 kg s^-2 K^-1
    mean_free_path_time = 0.2e-12
    thermal_sigma = np.sqrt(boltzmann * temperature / mass)
    #thermal_sigma = 7.5e3
    print(thermal_sigma)

    pt_list = [(-200, -100), (-20, -100), (-20, -60), (-10, -50), (10, -50), (20, -60), (20, -100), (200, -100), (200, 100), (20, 100), (20, 60), (10, 50), (-10, 50), (-20, 60), (-20, 100), (-200, 100)]
    island = [(10, 0), (100, 50), (100, -50)]
    world = World(pt_list)
    world.border = world.border + make_border(island)
    left_tele = Edge(Point(-200, -100), Point(-200, 100))
    left_tele.teleport_position_change = (400, 0)
    right_tele = Edge(Point(200, -100), Point(200, 100))
    right_tele.teleport_position_change = (-400, 0)
    world.teleport_edges.append(left_tele)
    world.teleport_edges.append(right_tele)

    electrons = generate_electron_list(world, 1000)
    run_simulation(world, electrons, 100, mean_free_path_time)
    print("Simulation complete")
    world.plot()
    plot_electron_list(electrons, True)
    plot_velocity_histogram(electrons, 100) # expecting 1.8e5
    make_heatmap(world, electrons, 10)
