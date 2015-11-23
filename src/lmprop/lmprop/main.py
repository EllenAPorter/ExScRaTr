import PIL.Image as pilImage
import tkinter as tk
import tkinter.filedialog as tkfd
import numpy as np
from pprint import pprint

from .ray import Ray
from .geometry import *
from .mesh import Mesh

EYE_FNAME = "/usr/local/share/lmprop/eye_rev.png"


# attached to all objeets that participate in voxel visualization
VOXELIZATION_TAG = 'VOXELIZATION'

# attached to all objeets that illustrate a normal (voxel-free) raytrace
RAYTRACE_TAG = 'RAYTRACE'

# all dashed lines drawn with this dash pattern
DASH_PATTERN = (10, 10)

# color schemed (for filled polygons)
BACKGROUND_HIGHLIGHT_COLOR = 'white'
BACKGROUND_NORMAL_COLOR    = 'gray80'
BLOCKING_POLYGON_COLOR     = 'gray40'
TARGET_POLYGON_COLOR       = 'blue'
LIGHT_COLOR                = 'yellow'
BLOCKED_RAY_COLOR          = 'gray60'
MISSED_PIXEL_COLOR         = 'white'


ITEM_LABEL_FONT = ('Times', 14)

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
        self.showVoxelization = tk.BooleanVar(value=False) # initially
        viewMenu.add_checkbutton(label="Voxelization",
                                 onvalue=True,
                                 offvalue=False,
                                 variable=self.showVoxelization,
                                 command=self.onViewVoxelization)
        self.showRaytrace = tk.BooleanVar(value=False) # initially
        viewMenu.add_checkbutton(label="Raytrace",
                                 onvalue=True,
                                 offvalue=False,
                                 variable=self.showRaytrace,
                                 command=self.onViewRaytrace)
        menuBar.add_cascade(label="View", menu=viewMenu)

        
        self.lmpropCanvas = LmpropCanvas(self, nCellRows, nCellColumns,
                viewPoint=LmpropCanvas.xyPosition(0.3, 1.25))
        self.lmpropCanvas.grid(row=1, column=0)

    def onFileExport(self):
        fname = tkfd.asksaveasfilename(
            defaultextension=".ps",
            filetypes=(
                ("Postscript", "*.ps" ),
                ("All",        "*" ),
            )
        )
        if fname:
            self.lmpropCanvas.postscript(file=fname)

    def onViewRaytrace(self):
        self.lmpropCanvas.setRaytraceVisibility(self.showRaytrace.get())

    def onViewVoxelization(self):
        self.lmpropCanvas.setVoxelizationVisibility(
            self.showVoxelization.get())


class ApplicationCanvas(tk.Canvas):
    """
    a wrapper around tk.Canvas

    This fixes some problems that Tkinter has with things like Numpy arrays
    and adds arrows.
    """
    def __init__(self, parent, **kwargs):
        super().__init__(parent, **kwargs)

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

    def drawRectangle(self, pUL, pLR, **kwargs):
        return self.create_rectangle(pUL[0], pUL[1], pLR[0], pLR[1], **kwargs)

    def drawText(self, p, text, **kwargs):
        return self.create_text(tuple(p), text=text, **kwargs)


class Wall(Mesh):

    def __init__(self, k, p0, p1):
        super().__init__(p0, p1, Cell.WALL_MESH_RESOLUTION)
        self.tag = '{}'.format(k)


class Cell:

    SIZE = 300 # width and height of (square) cells
    WALL_MESH_RESOLUTION = 9
    INSET = 1  # distinguish cell walls

    def __init__(self, lmpropCanvas, i, j, **kwargs):
        self.lmpropCanvas = lmpropCanvas
        insetVector = np.array((Cell.INSET, Cell.INSET))
        self.pUL = LmpropCanvas.xyPosition(j,   i  ) + insetVector
        self.pLR = LmpropCanvas.xyPosition(j+1, i+1) - insetVector
        self.walls = []
        self.tag = '{},{}'.format(i, j)
        for k, (x0, y0, x1, y1) in enumerate((
                (self.pUL[0], self.pUL[1], self.pUL[0], self.pLR[1]),
                (self.pUL[0], self.pLR[1], self.pLR[0], self.pLR[1]),
                (self.pLR[0], self.pLR[1], self.pLR[0], self.pUL[1]),
                (self.pLR[0], self.pUL[1], self.pUL[0], self.pUL[1]))):
            p0 = np.array((x0, y0))
            p1 = np.array((x1, y1))
            wall = Wall(k, p0, p1)
            self.walls.append(wall)

        self.backgroundId = self.lmpropCanvas.drawRectangle(self.pUL,
                                                            self.pLR, **kwargs)

        popupMenu = tk.Menu(self.lmpropCanvas, tearoff=0)
        popupMenu.add_checkbutton(label="Toggle Highlight",
                                  onvalue=True, offvalue=False,
                                  command=self.onToggleHighlight)

        def onPopupBackground(evt):
            try:
                popupMenu.tk_popup(evt.x_root, evt.y_root, 0)
            finally:
                popupMenu.grab_release()

        self.lmpropCanvas.tag_bind(self.backgroundId, '<Button-3>',
                                   onPopupBackground)

        self.highlight = False
        self.setHighlight()

    def onToggleHighlight(self):
        self.highlight = not self.highlight
        self.setHighlight()

    def setHighlight(self):
        self.lmpropCanvas.itemconfigure(self.backgroundId,
             fill=BACKGROUND_HIGHLIGHT_COLOR if self.highlight
                  else BACKGROUND_NORMAL_COLOR)


class LmpropCanvas(ApplicationCanvas):
    """an ApplicationCanvas specficially intended for this application
    """

    CELL_ARRAY_OFFSET = 30 # allow for labels
    LIGHT_RADIUS = 10

    def __init__(self, parent, nCellRows, nCellColumns, viewPoint):
        self.nCellRows = nCellRows
        self.nCellColumns= nCellColumns
        self.width = (self.nCellColumns * Cell.SIZE
                      + LmpropCanvas.CELL_ARRAY_OFFSET)
        self.height = (self.nCellRows * Cell.SIZE
                      + LmpropCanvas.CELL_ARRAY_OFFSET)
        super().__init__(parent, 
                         width=self.width, height=self.height,
                         background='white')

        self.corners = (
            (Cell.INSET, Cell.INSET),
            (Cell.INSET, self.height - Cell.INSET),
            (self.width - Cell.INSET, self.height - Cell.INSET),
            (self.width - Cell.INSET, Cell.INSET))

        self.cells = self.drawCells()

        #self.lightPosition = LmpropCanvas.xyPosition(0.3, 0.4)
        self.lightPosition = LmpropCanvas.xyPosition(0.3, 0.41)

        center = LmpropCanvas.xyPosition(1.2, 0.8)
        radius = 20
        self.blockingPolygon = createRegularPolygon(center, radius, 6)
        self.drawPolygon(self.blockingPolygon, fill=BLOCKING_POLYGON_COLOR,
                         outline='black', width=1)
        self.drawText(center + np.array((radius + 40, radius - 60)),
                      "blocking\nobject",
                      font=ITEM_LABEL_FONT, justify=tk.CENTER)

        illuminatedCell = self.cells[self.nCellRows-1, self.nCellColumns-1]
        self.illuminatedWall = illuminatedCell.walls[0]
        self.drawBlocking(0, 1, 0)

        center = LmpropCanvas.xyPosition(2.4, 1.5)
        radius = 80
        self.targetPolygon = createRegularPolygon(center, radius, 
                                                  10, phaseDeg=10)
        self.drawText(center + np.array((radius + 40, 40 - radius)),
                      "target\nobject",
                      font=ITEM_LABEL_FONT, justify=tk.CENTER)
        self.drawPolygon(self.targetPolygon,
                         fill=TARGET_POLYGON_COLOR, outline='black', width=1)

        self.drawRaytrace(viewPoint,
                          viewDirection=LmpropCanvas.xyDisplacement(0.6, 0.03),
                          fovDeg=30)

        self.drawCircle(self.lightPosition,
                        LmpropCanvas.LIGHT_RADIUS, fill=LIGHT_COLOR)
        self.drawText(self.lightPosition - np.array((0, 40)),
                      "point\nlight source",
                      font=ITEM_LABEL_FONT, justify=tk.CENTER)

        # need to save image reference so it won't be garbage-collected
        self.image = tk.PhotoImage(file=EYE_FNAME)
        self.drawText(viewPoint - np.array((0, 30)),
                      "viewer",
                      font=ITEM_LABEL_FONT, justify=tk.CENTER)
        self.drawImage(viewPoint, self.image)

        self.setVoxelizationVisibility(parent.showVoxelization.get())
        self.setRaytraceVisibility(parent.showRaytrace.get())
            
    def drawBlocking(self, iBlock, jBlock, kBlock):
        """display blocking of wall `kBlock` of cell (`iBlock`, `jBlock`)
        """
        blockingCell = self.cells[iBlock, jBlock]
        blockingWall = blockingCell.walls[kBlock]
        tag = blockingCell.tag + ',{}'.format(kBlock)
        for i in range(Cell.WALL_MESH_RESOLUTION):
            destinationPoint = self.illuminatedWall.samplePoints[i]
            o = self.lightPosition
            d = normalize(destinationPoint - o)
            ray = Ray(o, d)

            wallIntersection = blockingWall.intersectsRay(ray)
            if wallIntersection is not None:
                self.drawLine(
                    self.lightPosition, wallIntersection.p,
                    tag=VOXELIZATION_TAG)
                if 0: # if the blocking cell is the terminus
                  self.drawArrow(wallIntersection.p,
                               wallIntersection.p + 20*d, tag=VOXELIZATION_TAG)
                blockingIntersections = ray.intersectsPolygon(
                        self.blockingPolygon)
                if blockingIntersections:
                    self.drawArrow(wallIntersection.p, 
                                     blockingIntersections[0].p,
                                     tag=tag)

                    self.drawLine(blockingIntersections[1].p,
                                  destinationPoint,
                                  fill=BLOCKED_RAY_COLOR,
                                  tag=tag)
                else:
                    self.drawLine(wallIntersection.p,
                                  destinationPoint, tag=VOXELIZATION_TAG)
                    self.drawArrow(destinationPoint,
                                   destinationPoint + 20*d,
                                   tag=tag)
            else:
                self.drawLine(self.lightPosition, destinationPoint,
                              dash=DASH_PATTERN, tag=tag)

    def drawCells(self):
        cells = {}
        for i in range(self.nCellRows):
            for j in range(self.nCellColumns):
                cells[i,j] = Cell(self, i, j,
                            tag=VOXELIZATION_TAG)
        for j in range(self.nCellColumns):
            self.drawText(((j + 0.5) * Cell.SIZE
                           + LmpropCanvas.CELL_ARRAY_OFFSET,
                           0.5 * LmpropCanvas.CELL_ARRAY_OFFSET),
                          "{}".format(j), font=("Times", 16),
                          tag=VOXELIZATION_TAG)
        for i in range(self.nCellColumns):
            self.drawText((0.5 * LmpropCanvas.CELL_ARRAY_OFFSET,
                          (i + 0.5) * Cell.SIZE
                           + LmpropCanvas.CELL_ARRAY_OFFSET),
                          "{}".format(i), font=("Times", 16),
                          tag=VOXELIZATION_TAG)
        return cells

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
            self.drawArrow(ray.o, edgePoint, **kwargs)

    def drawRaytrace(self, viewPoint, viewDirection, fovDeg):
        d = mag(viewDirection)
        imageCenter = viewPoint + viewDirection
        normalizedViewDirection = normalize(viewDirection)

        # perpendicular to the view direction
        normalizedViewDirectionPerp = normalize(
            np.array((viewDirection[1], -viewDirection[0])))

        fov = fovDeg * np.pi / 180
        halfW = d * np.tan(fov/2)

        p0 = imageCenter + halfW * normalizedViewDirectionPerp
        p1 = imageCenter - halfW * normalizedViewDirectionPerp
        viewPlane = Mesh(p0, p1, 7)

        nSamples = len(viewPlane.samplePoints)
        w = 2 * halfW
        pixelSize = w / (nSamples - 1)

        for samplePoint in viewPlane.samplePoints:
            pixelPolygon = createOrientedSquare(samplePoint, pixelSize,
                                                viewDirection)
            viewRay = Ray(viewPoint, samplePoint - viewPoint)
            viewRayIntersections = viewRay.intersectsPolygon(
                    self.targetPolygon)
            if viewRayIntersections:
                self.drawArrow(viewPoint, viewRayIntersections[0].p,
                               tag=RAYTRACE_TAG)
                self.drawPolygon(pixelPolygon, fill=MISSED_PIXEL_COLOR,
                                 outline='black', width=1)
                # Draw it again in the target object color so it will
                # show up when raytracing is shown.
                self.drawPolygon(pixelPolygon, fill=TARGET_POLYGON_COLOR,
                                 outline='black', width=1,
                                 tag=RAYTRACE_TAG)
                blockingRay = Ray(viewRayIntersections[0].p,
                                  self.lightPosition
                                  - viewRayIntersections[0].p)
                blockingRayIntersections = \
                        blockingRay.intersectsPolygon(self.blockingPolygon)
                if blockingRayIntersections:
                    self.drawArrow(viewRayIntersections[0].p, 
                                   blockingRayIntersections[0].p,
                                   fill=BLOCKED_RAY_COLOR,
                                   tag=RAYTRACE_TAG)
                else:
                    endPoint = (self.lightPosition 
                            - LmpropCanvas.LIGHT_RADIUS 
                              * normalize(blockingRay.d))
                    self.drawArrow(viewRayIntersections[0].p,
                                   endPoint,
                                   tag=RAYTRACE_TAG)
            else:
                self.drawRay(viewRay, tag=RAYTRACE_TAG)
                self.drawPolygon(pixelPolygon, fill=MISSED_PIXEL_COLOR,
                                 outline='black', width=1)

        self.drawText(imageCenter + np.array((0, 80)),
                  "image", font=ITEM_LABEL_FONT, justify=tk.CENTER)

    def onMousePress(self, evt):
        print(self.find_closest(evt.x, evt.y))

    def setRaytraceVisibility(self, isVisible):
        if isVisible:
            self.itemconfigure(RAYTRACE_TAG, state=tk.NORMAL)
        else:
            self.itemconfigure(RAYTRACE_TAG, state=tk.HIDDEN)

    def setVoxelizationVisibility(self, isVisible):
        if isVisible:
            self.itemconfigure(VOXELIZATION_TAG, state=tk.NORMAL)
        else:
            self.itemconfigure(VOXELIZATION_TAG, state=tk.HIDDEN)

    @staticmethod
    def xyPosition(x, y):
        return np.array((x * Cell.SIZE + LmpropCanvas.CELL_ARRAY_OFFSET,
                         y * Cell.SIZE + LmpropCanvas.CELL_ARRAY_OFFSET),
                         dtype=float)

    @staticmethod
    def xyDisplacement(dx, dy):
        return np.array((dx * Cell.SIZE, dy * Cell.SIZE), dtype=float)


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
