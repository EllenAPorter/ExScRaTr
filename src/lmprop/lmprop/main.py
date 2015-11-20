import PIL.Image as pilImage
import tkinter as tk
import tkinter.filedialog as tkfd
import numpy as np
from pprint import pprint

from .ray import Ray
from .geometry import *

EYE_FNAME = "/usr/local/share/lmprop/eye_rev.png"


class Application(tk.Frame):

    def __init__(self, parent, nCellRows, nCellColumns, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        menuBar = tk.Menu(parent, tearoff=0)
        parent.config(menu=menuBar)

        fileMenu = tk.Menu(menuBar, tearoff=0)
        fileMenu.add_command(label="Export PS ...", command=self.onFileExport)
        fileMenu.add_command(label="Quit", command=self.quit)
        menuBar.add_cascade(label="File", menu=fileMenu)

        viewMenu = tk.Menu(menuBar, tearoff=0)
        self.showVoxelization = tk.BooleanVar(value=True) # initially
        viewMenu.add_checkbutton(label="Voxelization",
                                 onvalue=True,
                                 offvalue=False,
                                 variable=self.showVoxelization,
                                 command=self.onViewVoxelization)
        menuBar.add_cascade(label="View", menu=viewMenu)

        self.canvas = Canvas(self, nCellRows, nCellColumns)
        self.canvas.grid(row=1, column=0)

    def onFileExport(self):
        fname = tkfd.asksaveasfilename(
            defaultextension=".ps",
            filetypes=(
                ("Postscript", "*.ps" ),
                ("All",        "*" ),
            )
        )
        if fname:
            self.canvas.postscript(file=fname)

    def onViewVoxelization(self):
        self.canvas.setVoxelizationVisibility(self.showVoxelization.get())


class Canvas(tk.Canvas):

    VOXELIZATION = 'VOXELIZATION'

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

        self.cells = self.drawCells()

        pLight = self.xyPosition(0.3, 0.4)
        illuminatedCell = self.cells[self.nCellRows-1, self.nCellColumns-1]
        blockingCell = self.cells[0, 1]
        illuminatedWall = illuminatedCell.walls[0]

        blockingPolygon = createRegularPolygon(self.xyPosition(1.2, 0.8),
                                               20, 6)
        self.drawPolygon(blockingPolygon, fill='cyan', outline='black',
                         width=1)

        blockingCell.drawBlocking(pLight, 0, blockingPolygon, illuminatedWall,
                                  tag=Canvas.VOXELIZATION)

        self.drawCircle(pLight, 5, fill='yellow')

        targetPolygon = createRegularPolygon(self.xyPosition(2.4, 1.5), 80, 10,
                                             phaseDeg=10)
        self.drawPolygon(targetPolygon, fill='blue', outline='black', width=1)

        self.drawRayTrace(eyePoint=self.xyPosition(0.3, 1.25),
                          viewDirection=self.xyDirection(0.6, 0.03),
                          pLight=pLight,
                          fovDeg=30,
                          targetPolygon=targetPolygon,
                          blockingPolygon=blockingPolygon)
        self.setVoxelizationVisibility(True)
            
                    
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

    def drawCells(self):
        cells = {}
        for i in range(self.nCellRows):
            for j in range(self.nCellColumns):
                cells[i,j] = Cell(self, i, j, tag=Canvas.VOXELIZATION)
        for j in range(self.nCellColumns):
            self.drawText(((j + 0.5) * Cell.wh + Cell.offset,
                           0.5 * Cell.offset),
                          "{}".format(j), font=("Times", 16),
                          tag=Canvas.VOXELIZATION)
        for i in range(self.nCellColumns):
            self.drawText((0.5 * Cell.offset,
                          (i + 0.5) * Cell.wh + Cell.offset),
                          "{}".format(i), font=("Times", 16),
                          tag=Canvas.VOXELIZATION)
        return cells

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

    def drawRayTrace(self, eyePoint, viewDirection, pLight, fovDeg,
                     targetPolygon, blockingPolygon):
        d = mag(viewDirection)
        pViewCenter = eyePoint + viewDirection
        normalizedViewDirection = normalize(viewDirection)

        # perpendicular to the view direction
        normalizedViewDirectionPerp = normalize(
            np.array((viewDirection[1], -viewDirection[0])))

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
            lightRay = Ray(eyePoint, samplePoint - eyePoint)
            lightRayIntersections = lightRay.intersectsPolygon(targetPolygon)
            if lightRayIntersections:
                self.drawArrow(eyePoint, lightRayIntersections[0].p)
                self.drawPolygon(pixelPolygon, fill='blue',
                                 outline='black', width=1)
                blockingRay = Ray(lightRayIntersections[0].p,
                                  pLight - lightRayIntersections[0].p)
                blockingRayIntersections = \
                        blockingRay.intersectsPolygon(blockingPolygon)
                if 0:
                    print(blockingRay.o)
                    print(blockingRay.d)
                    print(blockingRayIntersections)
                    print()
                if blockingRayIntersections:
                    self.drawArrow(lightRayIntersections[0].p, 
                                   blockingRayIntersections[0].p,
                                   fill='gray60')
                else:
                    self.drawArrow(lightRayIntersections[0].p, pLight)
            else:
                self.drawRay(lightRay)
                self.drawPolygon(pixelPolygon, fill='white',
                                 outline='black', width=1)

        # need to save image reference so it won't be garbage-collected
        self.image = tk.PhotoImage(file=EYE_FNAME)
        self.drawImage(eyePoint, self.image)

    def drawRectangle(self, pUL, pLR, **kwargs):
        if 'fill' not in kwargs:
            kwargs['fill'] = 'light green' # default
        return self.create_rectangle(pUL[0], pUL[1], pLR[0], pLR[1], **kwargs)

    def drawText(self, p, text, **kwargs):
        return self.create_text(p, text=text, **kwargs)

    def setVoxelizationVisibility(self, isVisible):
        if isVisible:
            self.itemconfigure(Canvas.VOXELIZATION, state=tk.NORMAL)
        else:
            self.itemconfigure(Canvas.VOXELIZATION, state=tk.HIDDEN)

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

    def __init__(self, canvas, i, j, **kwargs):
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
            wall = Mesh(self.canvas, p0, p1, Cell.wallMeshResolution, **kwargs)
            self.walls.append(wall)

        pUL = np.array((self.pUL[0], self.pUL[1]))
        pLR = np.array((self.pLR[0], self.pLR[1]))

        self.backgroundId = self.canvas.drawRectangle(pUL, pLR, **kwargs)

        self.highlight = False
        self.setHighlight()

    def drawBlocking(self, pLight, wallIndex, blockingPolygon, wallTo, 
                     **kwargs):
        intercedingWall = self.walls[wallIndex]
        for i in range(Cell.wallMeshResolution):
            pTo = wallTo.samplePoints[i]
            o = pLight
            d = normalize(pTo - o)
            ray = Ray(o, d)

            wallIntersection = intercedingWall.intersectsRay(ray)
            if wallIntersection is not None:
                self.canvas.drawLine(pLight, wallIntersection.p, **kwargs)
                self.canvas.drawArrow(wallIntersection.p,
                                      wallIntersection.p + 20*d, **kwargs)
                blockingIntersections = ray.intersectsPolygon(blockingPolygon)
                if blockingIntersections:
                    self.canvas.drawArrow(wallIntersection.p, 
                                     blockingIntersections[0].p, **kwargs)
                    self.canvas.drawLine(blockingIntersections[1].p,
                                    pTo, fill='gray60', **kwargs)
                else:
                    self.canvas.drawLine(wallIntersection.p, pTo, **kwargs)
                    self.canvas.drawArrow(pTo, pTo + 20*d, **kwargs)
            else:
                self.canvas.drawLine(pLight, pTo, fill='gray60', **kwargs)

    def setHighlight(self):
        self.canvas.itemconfigure(self.backgroundId,
                                  fill= 'white' if self.highlight
                                                else 'gray80')


class Mesh:
    """This is a 1D mesh analogous to a more conventional 2D or 3D mesh.
    """

    def __init__(self, canvas, p0, p1, nSamples, **kwargs):
        self.canvas = canvas
        self.p0 = p0
        self.p1 = p1
        d = 1 / (nSamples - 1) # spacing between vertices
        self.samplePoints = []
        for i in range(nSamples):
            t = i * d
            p = self.p0 + t * (self.p1 - self.p0)
            self.samplePoints.append(p)
        self.canvas.drawLine(self.p0, self.p1, **kwargs)

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


def main():
    root = tk.Tk()
    root.bind("q", lambda event: root.quit())

    application = Application(root,
                              nCellRows = 2,
                              nCellColumns = 3)
    application.grid(row=0, column=0)

    root.mainloop()


if __name__ == '__main__':
    main()
