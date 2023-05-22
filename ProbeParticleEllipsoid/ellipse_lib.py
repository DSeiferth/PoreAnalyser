import numpy as np
import math
import matplotlib.pyplot as plt

class atom:
    def __init__(self, x, y, z=0, r=1):
        self.x = x
        self.y = y
        self.z = z
        self.r = r

class ellipse:
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
    if plot: ax.plot(closest_rt[0], closest_rt[1], '-x', color='black')
    
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
