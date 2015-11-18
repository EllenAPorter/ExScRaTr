import PIL.Image as pilImage
import tkinter as tk
import numpy as np
from pprint import pprint

from .ray import Ray
from .geometry import *

EYE_FNAME = "/usr/local/share/lmprop/eye_rev.png"


class Application(tk.Frame):

    def __init__(self, parent, nCellRows, nCellColumns, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.canvas = Canvas(self, nCellRows, nCellColumns)
        self.canvas.grid(row=0, column=0)


class Canvas(tk.Canvas):

    def __init__(self, parent, nCellRows, nCellColumns):
        self.nCellRows = nCellRows
        self.nCellColumns= nCellColumns
        self.width = self.nCellColumns * Cell.wh + Cell.offset
        self.height = self.nCellRows * Cell.wh + Cell.offset
        super().__init__(parent, 
                         width=self.width, height=self.height,
                         background='white')

        self.corners = (
            (Cell.inset, Cell.inset),
            (Cell.inset, self.height - Cell.inset),
            (self.width - Cell.inset, self.height - Cell.inset),
            (self.width - Cell.inset, Cell.inset))

        self.cells = {}
        for i in range(self.nCellRows):
            for j in range(self.nCellColumns):
                self.cells[i,j] = Cell(self, i, j)
        for j in range(self.nCellColumns):
            self.drawText(((j + 0.5) * Cell.wh + Cell.offset,
                           0.5 * Cell.offset),
                          "{}".format(j), font=("Times", 16))
        for i in range(self.nCellColumns):
            self.drawText((0.5 * Cell.offset,
                          (i + 0.5) * Cell.wh + Cell.offset),
                          "{}".format(i), font=("Times", 16))
            
                    
    def drawArrow(self, p0, p1, **kwargs):
        d = 15 * normalize(p1 - p0)
        dPerp = np.array((d[1], -d[0]))
        polygon = ( p1,
                    p1 - d - dPerp/4,
                    p1 - 0.75 * d,
                    p1 - d + dPerp/4 )
        if 'fill' not in kwargs:
            kwargs['fill'] = 'black'
        self.drawPolygon(polygon, **kwargs)
        self.create_line(p0[0], p0[1], p1[0], p1[1], **kwargs)

    def drawCircle(self, pCtr, r, **kwargs):
        if 'fill' not in kwargs:
            kwargs['fill'] = 'light green' # default
        return self.create_oval((pCtr[0] - r, pCtr[1] - r),
                                (pCtr[0] + r, pCtr[1] + r), **kwargs)

    def drawImage(self, p, image, **kwargs):
        return self.create_image(p[0], p[1], image=image, **kwargs)

    def drawLine(self, p0, p1, **kwargs):
        return self.create_line(p0[0], p0[1], p1[0], p1[1], **kwargs)

    def drawPolygon(self, polygon, **kwargs):
        # convert arrays, if need be
        return self.create_polygon([ tuple(p) for p in polygon ],
                                   **kwargs)

    def drawRay(self, ray, **kwargs):
        edgePoint = None
        for i in range(4):
            p0 = self.corners[i]
            p1 = self.corners[(i+1) % 4]
            (p0, n) = lineThroughPoints(p0, p1)
            intersection = ray.intersectsLine(p0, n)
            if intersection is not None:
                if edgePoint is None \
                       or mag(intersection.p - ray.o) < mag(edgePoint - ray.o):
                    edgePoint = intersection.p
        if edgePoint is not None:
            self.drawArrow(ray.o, edgePoint)

    def drawRectangle(self, pUL, pLR, **kwargs):
        if 'fill' not in kwargs:
            kwargs['fill'] = 'light green' # default
        return self.create_rectangle(pUL[0], pUL[1], pLR[0], pLR[1], **kwargs)

    def drawRayTrace(self, eyePoint, viewDirection, fovDeg, polygon):
        d = mag(viewDirection)
        pViewCenter = eyePoint + viewDirection
        normalizedViewDirection = normalize(viewDirection)

        # perpendicular to the view direction
        normalizedViewDirectionPerp = normalize(np.array((viewDirection[1], 
                                                          -viewDirection[0])))

        fov = fovDeg * np.pi / 180
        halfW = d * np.tan(fov/2)

        p0 = pViewCenter + halfW * normalizedViewDirectionPerp
        p1 = pViewCenter - halfW * normalizedViewDirectionPerp
        viewPlane = Mesh(self, p0, p1, 7)

        nSamples = len(viewPlane.samplePoints)
        w = 2 * halfW
        pixelWH = w / (nSamples - 1)

        for samplePoint in viewPlane.samplePoints:
            pixelPolygon = createOrientedSquare(samplePoint, pixelWH,
                                                viewDirection)
            ray = Ray(eyePoint, samplePoint - eyePoint)
            intersections = ray.intersectsPolygon(polygon)
            if intersections:
                self.drawArrow(eyePoint, intersections[0].p)
                self.drawPolygon(pixelPolygon, fill='blue',
                                 outline='black', width=1)
            else:
                self.drawRay(ray)
                self.drawPolygon(pixelPolygon, fill='gray60',
                                 outline='black', width=1)

        # need to save image reference so it won't be garbage-collected
        self.image = tk.PhotoImage(file=EYE_FNAME)
        self.drawImage(eyePoint, self.image)

    def drawText(self, p, text, **kwargs):
        return self.create_text(p, text=text, **kwargs)

    def xyPosition(self, x, y):
        return np.array((x * Cell.wh + Cell.offset,
                         y * Cell.wh + Cell.offset), dtype=float)

    def xyDirection(self, dx, dy):
        return np.array((dx * Cell.wh, dy * Cell.wh), dtype=float)


class Cell:

    wh = 300 # width and height of (square) cells
    wallMeshResolution = 9
    whRect = 10 # w & h of rectangle noting light source sample
    inset = 1  # distinguish cell walls
    offset = 30 # allow for labels

    def __init__(self, canvas, i, j):
        self.pUL = np.array((j       * Cell.wh + Cell.inset + Cell.offset,
                             i       * Cell.wh + Cell.inset + Cell.offset))
        self.pLR = np.array(((j + 1) * Cell.wh - Cell.inset + Cell.offset,
                             (i + 1) * Cell.wh - Cell.inset + Cell.offset))
        self.canvas = canvas
        self.walls = []
        for (x0, y0, x1, y1) in (
                (self.pUL[0], self.pUL[1], self.pUL[0], self.pLR[1]),
                (self.pUL[0], self.pLR[1], self.pLR[0], self.pLR[1]),
                (self.pLR[0], self.pLR[1], self.pLR[0], self.pUL[1]),
                (self.pLR[0], self.pUL[1], self.pUL[0], self.pUL[1])):
            p0 = np.array((x0, y0))
            p1 = np.array((x1, y1))
            wall = Mesh(self.canvas, p0, p1, Cell.wallMeshResolution)
            self.walls.append(wall)

        pUL = np.array((self.pUL[0], self.pUL[1]))
        pLR = np.array((self.pLR[0], self.pLR[1]))
        self.backgroundId = self.canvas.drawRectangle(pUL, pLR)
        self.setHighlight(False)

    def drawBlocking(self, pLight, wallIndex, blockingPolygon, wallTo):
        intercedingWall = self.walls[wallIndex]
        for i in range(Cell.wallMeshResolution):
            pTo = wallTo.samplePoints[i]
            o = pLight
            d = normalize(pTo - o)
            ray = Ray(o, d)

            wallIntersection = intercedingWall.intersectsRay(ray)
            if wallIntersection is not None:
                self.canvas.drawLine(pLight, wallIntersection.p)
                self.canvas.drawArrow(wallIntersection.p,
                                      wallIntersection.p + 20*d)
                blockingIntersections = ray.intersectsPolygon(blockingPolygon)
                if blockingIntersections:
                    self.canvas.drawArrow(wallIntersection.p, 
                                     blockingIntersections[0].p)
                    self.canvas.drawLine(blockingIntersections[1].p,
                                    pTo, fill='gray60')
                else:
                    self.canvas.drawLine(wallIntersection.p, pTo)
                    self.canvas.drawArrow(pTo, pTo + 20*d)
            else:
                self.canvas.drawLine(pLight, pTo, fill='gray60')

    def setHighlight(self, status):
        self.canvas.itemconfigure(self.backgroundId,
                                  fill= 'white' if status else 'gray80')


class Mesh:
    """This is a 1D mesh analogous to a more conventional 2D or 3D mesh.
    """

    def __init__(self, canvas, p0, p1, nSamples):
        self.canvas = canvas
        self.p0 = p0
        self.p1 = p1
        d = 1 / (nSamples - 1) # spacing between vertices
        self.samplePoints = []
        for i in range(nSamples):
            t = i * d
            p = self.p0 + t * (self.p1 - self.p0)
            self.samplePoints.append(p)
        self.canvas.drawLine(self.p0, self.p1)

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


def drawFrame0(canvas):
    pLight = canvas.xyPosition(0.3, 0.4)
    illuminatedCell = canvas.cells[canvas.nCellRows-1, canvas.nCellColumns-1]
    blockingCell = canvas.cells[0, 1]
    illuminatedWall = illuminatedCell.walls[0]

    blockingCell.setHighlight(True)

    blockingPolygon = createRegularPolygon(canvas.xyPosition(1.2, 0.8), 15, 6)
    canvas.drawPolygon(blockingPolygon, fill='cyan', outline='black', width=1)

    blockingCell.drawBlocking(pLight, 0, blockingPolygon, illuminatedWall)

    canvas.drawCircle(pLight, 5, fill='yellow')

    targetPolygon = createRegularPolygon(canvas.xyPosition(2.4, 1.5), 80, 10,
                                         phaseDeg=10)
    canvas.drawPolygon(targetPolygon, fill='blue', outline='black', width=1)

    canvas.drawRayTrace(eyePoint = canvas.xyPosition(0.3, 1.25),
                        viewDirection = canvas.xyDirection(0.6, 0.03),
                        fovDeg = 30,
                        polygon = targetPolygon)


def main():
    root = tk.Tk()
    root.bind("q", lambda event: root.quit())

    def snapshot(event):
        print("ok")
        application.canvas.postscript(file="frame0.ps")
    root.bind("p", snapshot)

    application = Application(root,
                              nCellRows = 2,
                              nCellColumns = 3)
    application.grid(row=0, column=0)

    drawFrame0(application.canvas)
    root.mainloop()


if __name__ == '__main__':
    main()
