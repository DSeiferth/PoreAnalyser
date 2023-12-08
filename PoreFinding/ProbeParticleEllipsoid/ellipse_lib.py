import numpy as np
import math
import matplotlib.pyplot as plt

class atom:
    """
    Class representing a 3D atom with coordinates (x, y, z) and a radius (r).

    Parameters:
    - x (float): x-coordinate of the atom.
    - y (float): y-coordinate of the atom.
    - z (float, optional): z-coordinate of the atom. Default is 0.
    - r (float, optional): Radius of the atom. Default is 1.

    Attributes:
    - x (float): x-coordinate of the atom.
    - y (float): y-coordinate of the atom.
    - z (float): z-coordinate of the atom.
    - r (float): Radius of the atom.

    Example:
    >>> my_atom = Atom(x=1.0, y=2.0, z=0.0, r=1.5)
    >>> print(my_atom.x, my_atom.y, my_atom.z, my_atom.r)
    1.0 2.0 0.0 1.5
    """
    def __init__(self, x, y, z=0, r=1):
        self.x = x
        self.y = y
        self.z = z
        self.r = r

class ellipse:
    """
    Class representing a 2D ellipse with parameters (a, b, cx, cy, cz, r, theta).

    Parameters / Attributes:
    - a (float): Length of the semi-major axis.
    - b (float): Length of the semi-minor axis (radius to grow).
    - cx (float): x-coordinate of the center.
    - cy (float): y-coordinate of the center.
    - cz (float, optional): z-coordinate of the center. Default is 0.
    - r (float, optional): Radius of the ellipse. Default is 1.
    - theta (float, optional): Angle of rotation in radians. Default is 0.

    Methods:
    - on_ellipse(x, y): Check if a point (x, y) is on the ellipse.
    - draw(res=0.01): Generate coordinates of the ellipse for plotting.

    Example:
    >>> my_ellipse = Ellipse(a=3.0, b=2.0, cx=0.0, cy=0.0, theta=0.0)
    >>> print(my_ellipse.on_ellipse(1.0, 1.0))
    True
    >>> x_coords, y_coords = my_ellipse.draw(res=0.01)
    """
    def __init__(self, a, b, cx, cy, cz=0, r=1, theta=0):
        self.a = a
        self.b = b ### radius to grow ###
        self.cx = cx
        self.cy = cy
        self.cz = cz
        self.r = r
        self.theta = theta
    def on_ellipse(self, x, y):
        ### rotate
        x1 = x*np.cos(-self.theta) - y*np.sin(-self.theta)
        y1 = x*np.sin(-self.theta) + y*np.cos(-self.theta)
        #return x1*x1/(a*a) + y1*y1/(b*b) <= 1
        c1r = self.cx*np.cos(-self.theta) - self.cy*np.sin(-self.theta)
        c2r = self.cx*np.sin(-self.theta) + self.cy*np.cos(-self.theta)
        return (x1-c1r)*(x1-c1r)/(self.a*self.a) + (y1-c2r)*(y1-c2r)/(self.b*self.b) <= 1
    def draw(self, res=0.01):
        t = np.arange(0, 2*np.pi, res)
        x = self.a * np.cos(t) #+ c1
        y = self.b * np.sin(t) #+ c2

        x1 = x*np.cos(self.theta) - y*np.sin(self.theta) + self.cx
        y1 = x*np.sin(self.theta) + y*np.cos(self.theta) + self.cy
        return x1, y1
    
def distance_ellipse(semi_major, semi_minor, p):
    """
    Calculate the closest point on an ellipse to a given point in 2D space.

    Parameters:
    - semi_major (float): Length of the semi-major axis of the ellipse.
    - semi_minor (float): Length of the semi-minor axis of the ellipse.
    - p (tuple): A tuple representing the coordinates (x, y) of the point in 2D space.

    Returns:
    tuple: A tuple representing the coordinates (x, y) of the closest point on the ellipse to the given point.

    Reference:
    This function is based on the method described in the following blog post:
    "A simple method for distance to ellipse"
    https://blog.chatfield.io/simple-method-for-distance-to-ellipse/

    Example:
    >>> distance_ellipse(3.0, 2.0, (1.0, 1.0))
    (1.2487110341841325, 1.8185123044348084)
    """
    ### https://blog.chatfield.io/simple-method-for-distance-to-ellipse/ ###
    px = abs(p[0])
    py = abs(p[1])

    t = math.pi / 4

    a = semi_major
    b = semi_minor

    for x in range(0, 3):
        x = a * math.cos(t)
        y = b * math.sin(t)

        ex = (a*a - b*b) * math.cos(t)**3 / a
        ey = (b*b - a*a) * math.sin(t)**3 / b

        rx = x - ex
        ry = y - ey

        qx = px - ex
        qy = py - ey

        r = math.hypot(ry, rx)
        q = math.hypot(qy, qx)

        delta_c = r * math.asin((rx*qy - ry*qx)/(r*q))
        delta_t = delta_c / math.sqrt(a*a + b*b - x*x - y*y)

        t += delta_t
        t = min(math.pi/2, max(0, t))

    return (math.copysign(x, p[0]), math.copysign(y, p[1]))

def dist_ellipse_vdwSphere(ellipse, sphere, plot=0):
    """
    Calculate the distance between an ellipse and a sphere, considering the van der Waals (vdW) radius of the sphere.

    Parameters:
    - ellipse (Ellipse): An Ellipse object representing the ellipse in 2D space.
    - sphere (Sphere): A Sphere object representing the sphere in 2D space, with attributes x, y (center coordinates) and r (radius).
    - plot (int, optional): An integer flag indicating whether to plot the result. Default is 0 (no plot).

    Returns:
    float: The distance between the ellipse and the sphere, taking into account the van der Waals radius.

    Notes:
    The function first checks if the center of the sphere is inside the ellipse. If so, the distance is considered as the negative of the sphere's radius.

    The distance is calculated by transforming the coordinates of the sphere to the ellipse's local coordinate system, finding the closest point on the ellipse,
    and then transforming the closest point back to the global coordinate system. The distance is the difference between the transformed closest point and the sphere's center,
    considering the sphere's radius.

    The van der Waals (vdW) radius of the sphere is taken into account when calculating the distance.

    Example 1:
    >>> ellipse = e_lib.ellipse(a=3.0, b=2.0, cx=0.0, cy=0.0, theta=0.0)
    >>> sphere = e_lib.atom(x=1.0, y=1.0, r=0.5)
    >>>  e_lib.dist_ellipse_vdwSphere(ellipse, sphere, plot=0)
    -0.5 # distance between the ellipse and the sphere

    Example 2:
    >>> ellipse = e_lib.ellipse(a=3.0, b=2.0, cx=5.0, cy=5.0, theta=0.0)
    >>> sphere = e_lib.atom(x=1.0, y=1.0, r=0.5)
    >>> e_lib.dist_ellipse_vdwSphere(ellipse, sphere, plot=0)
    2.6950072040653335
    """
    ### test if center in ellispse ###
    if ellipse.on_ellipse(sphere.x,sphere.y):
        #print('vdw center in ellipse')
        return -sphere.r
    
    center = np.array([ellipse.cx, ellipse.cy])
    rot = np.array([[np.cos(ellipse.theta), -np.sin(ellipse.theta)],
                [np.sin(ellipse.theta),np.cos(ellipse.theta)]])
    rot_back = np.array([[np.cos(-ellipse.theta), -np.sin(-ellipse.theta)],
                     [np.sin(-ellipse.theta),np.cos(-ellipse.theta)]])
    ### transform point to ellipse coordinates ###
    ### 1. translate 2. rotate -angle ###
    point = np.array([sphere.x, sphere.y])
    point_t = point - center 
    point_tr = rot_back.dot(point_t) 
    
    c1, c2 = distance_ellipse(semi_major=ellipse.a, semi_minor=ellipse.b, p=point_tr)
    closest = np.array([c1,c2]) ### in ellipse coordinates ###
    ### transform closest point from ellipse to cartesian coordiantes ###
    ### 1. rotate + angle 2. translate ###
    closest_r = rot.dot(closest)
    closest_rt = closest_r + center
    #print('closest_rt', closest_rt)
    if plot: plt.plot(closest_rt[0], closest_rt[1], '-x', color='black')
    
    ### overlap if radius larger than radius ### 
    return np.linalg.norm(closest_rt-point) - sphere.r

def assign_radius(a_type):
    if a_type=='C':
        return 1.85
    elif a_type=='O':
        return 1.65
    elif a_type=='S':
        return 2.00
    elif a_type=='N':
        return 1.75
    elif a_type=='H':
        return 1.00
    elif a_type=='P':
        return 2.10
