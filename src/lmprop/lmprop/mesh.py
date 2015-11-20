from .geometry import *

class Mesh:
    """This is a 1D mesh analogous to a more conventional 2D or 3D mesh.
    """

    def __init__(self, p0, p1, nSamples):
        self.p0 = p0
        self.p1 = p1
        d = 1 / (nSamples - 1) # spacing between vertices
        self.samplePoints = []
        for i in range(nSamples):
            t = i * d
            p = self.p0 + t * (self.p1 - self.p0)
            self.samplePoints.append(p)

    def line(self):
        """returns the line the wall intersects

        it returns (p0, n) where p0 is an arbitrary point on the line and n is
        a normal to the the line, which is defined by n . (p - p0) = 0 for any
        point p on the line.
        """
        return lineThroughPoints(self.p0, self.p1)

    def intersectsRay(self, ray):
        """returns true iff ray intersects the wall
        """
        return ray.intersectsLineSegment(self.p0, self.p1)
