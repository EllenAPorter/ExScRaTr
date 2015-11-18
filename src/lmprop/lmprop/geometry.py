import numpy as np


def createOrientedSquare(centerPoint, wh, u):
    """create a wh x wh square centered on centerPoint w/one edge || to u
    """
    u /= mag(u)
    v = np.array((-u[1], u[0]))
    return (centerPoint - wh/2 * u - wh/2 * v,
            centerPoint + wh/2 * u - wh/2 * v,
            centerPoint + wh/2 * u + wh/2 * v,
            centerPoint - wh/2 * u + wh/2 * v)


def createRegularPolygon(centerPoint, radius, nSides, phaseDeg=0.0):
    dTheta = 2 * np.pi / nSides
    phase = phaseDeg * np.pi / 180
    polygon = []
    for i in range(nSides):
        theta = i * dTheta
        dx = radius * np.cos(theta - phase)
        dy = radius * np.sin(theta - phase)
        p = centerPoint + np.array((dx, -dy))
        polygon.append(p)
    return polygon


def lineThroughPoints(p0, p1):
    """returns the line that passes through points p0 and p1

    it returns (p0, n) where p0 is an arbitrary point on the line and n is
    a normal to the the line, which is defined by n . (p - p0) = 0 for any
    point p on the line.
    """
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    # If dx & dy are positive, the positive half-plane is SE of the line.
    mag = (dx**2 + dy**2)**0.5
    n = (dy/mag, -dx/mag)
    return (p0, n)


def mag(v):
    return np.linalg.norm(v)


def normalize(v):
    return v / mag(v)
