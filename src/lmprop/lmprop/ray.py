import collections
import numpy as np

from .geometry import *

EPS = 1.0e-8 # (feeble) attempt to account for roundoff


Intersection = collections.namedtuple('Intersection', ('p', 't'))


class Ray:

    def __init__(self, o, d):
        """the ray is of the form o+dt
        """
        self.o = o
        self.d = d

    def intersectsLine(self, p0, n):
        """returns the intersection of the ray with a line (or None)

        the line is (p-p0).n = 0
        """
        nDotD = np.dot(n, self.d)
        if np.abs(nDotD) < EPS:
            # If the ray begins on the line, its direction doesn't matter.
            if mag(self.o - p0) < EPS:
                return Intersection(self.o, 0.0)
            else:
                return None
        t = np.dot(n, p0 - self.o) / nDotD
        if t < 0: # t must be positive for the ray to intersect
            return None
        p = self.o + self.d * t
        return Intersection(p, t)

    def intersectsLineSegment(self, p0, p1):
        (pOnLine, n) = lineThroughPoints(p0, p1)
        intersection = self.intersectsLine(pOnLine, n)
        if intersection is None:
            return None

        # At this point, the ray intersects the line defined by the
        # points, but to intersect the line segment, the intersection must
        # lie between them. We'll do this by seeing if the vectors from
        # the endpoints to the intersection are antiparallel.
        v0 = intersection.p - p0
        v1 = intersection.p - p1
        result = intersection if np.dot(v0, v1) <= 0 else None
        return result

    def intersectsPolygon(self, polygon):
        """returns a list of intersections of the ray with a polygon (or [])
        """
        n = len(polygon)
        result = []
        for i in range(n):
            p0 = polygon[i]
            p1 = polygon[(i+1) % n]
            intersection = self.intersectsLineSegment(p0, p1)
            if intersection:
                result.append(intersection)
        def getT(intersection):
            return intersection.t
        return sorted(result, key=getT)
