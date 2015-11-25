import PIL.Image as pilImage
import tkinter as tk
import tkinter.filedialog as tkfd
import numpy as np
from pprint import pprint

from .ray import Ray
from .geometry import *
from .mesh import Mesh

EYE_FNAME = "/usr/local/share/lmprop/eye_rev.png"


# default line width (allow for both publication and presentations)
LINE_WIDTH = 2

# attached to all objeets that participate in voxel visualization
# (light or view)
VOXELIZATION_TAG = 'VOXELIZATION'

# attached to all objeets that participate in voxel visualization of
# light propagation
LIGHT_VOXELIZATION_TAG = 'LIGHT_VOXELIZATION'

# attached to all objeets that participate in voxel visualization of
# the view
VIEW_VOXELIZATION_TAG = 'VIEW_VOXELIZATION'

# attached to all objeets that illustrate a normal (voxel-free) raytrace
RAYTRACE_TAG = 'RAYTRACE'

# all dashed lines drawn with this dash pattern
DASH_PATTERN = (10, 10)

# color schemes (for filled polygons)
BACKGROUND_HIGHLIGHT_COLOR = 'white'
BACKGROUND_NORMAL_COLOR    = 'gray80'
BLOCKING_POLYGON_COLOR     = 'gray40'
TARGET_POLYGON_COLOR       = 'blue'
LIGHT_COLOR                = 'yellow'
BLOCKED_RAY_COLOR          = 'gray50'
MISSED_PIXEL_COLOR         = 'white'

# These are the choices for view mode.
SCENE_ONLY_VIEW_MODE = "Scene Only"
RAYTRACE_VIEW_MODE = "Raytrace"
VOXELIZATION_VIEW_MODE = "Voxelization"
LIGHT_VOXELIZATION_VIEW_MODE = "Voxelized Light Propagation"
VIEW_VOXELIZATION_VIEW_MODE = "Voxelized View Propagation"
        

ITEM_LABEL_FONT = ('Times', 14)

class Application(tk.Frame):

    def __init__(self, parent, root, nCellRows, nCellColumns, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)

        # Note that we can't bind keyboard events to `self`, since
        # Frames can't get keyboard events.
        root.bind("q",
                    lambda event: root.quit())
        root.bind("s",
                    lambda event: self.setViewMode(SCENE_ONLY_VIEW_MODE))
        root.bind("r",
                    lambda event: self.setViewMode(RAYTRACE_VIEW_MODE))
        root.bind("l",
                    lambda event:
                        self.setViewMode(LIGHT_VOXELIZATION_VIEW_MODE))
        root.bind("v",
                    lambda event: 
                        self.setViewMode(VIEW_VOXELIZATION_VIEW_MODE))
        root.bind("V",
                    lambda event: 
                        self.setViewMode(VOXELIZATION_VIEW_MODE))

        menuBar = tk.Menu(root, tearoff=0)
        root.config(menu=menuBar)

        fileMenu = tk.Menu(menuBar, tearoff=0)
        fileMenu.add_command(label="Export PS ...", command=self.onFileExport)
        fileMenu.add_command(label="Quit", command=self.quit)
        menuBar.add_cascade(label="File", menu=fileMenu)

        viewMenu = tk.Menu(menuBar, tearoff=0)

        self.viewModeVar = tk.StringVar(value=SCENE_ONLY_VIEW_MODE)
        viewMenu.add_radiobutton(label=SCENE_ONLY_VIEW_MODE,
                                 variable=self.viewModeVar,
                                 command=self.onViewMode)
        viewMenu.add_radiobutton(label=RAYTRACE_VIEW_MODE,
                                 variable=self.viewModeVar,
                                 command=self.onViewMode)
        viewMenu.add_radiobutton(label=VOXELIZATION_VIEW_MODE,
                                 variable=self.viewModeVar,
                                 command=self.onViewMode)
        viewMenu.add_radiobutton(label=LIGHT_VOXELIZATION_VIEW_MODE,
                                 variable=self.viewModeVar,
                                 command=self.onViewMode)
        viewMenu.add_radiobutton(label=VIEW_VOXELIZATION_VIEW_MODE,
                                 variable=self.viewModeVar,
                                 command=self.onViewMode)
        menuBar.add_cascade(label="View", menu=viewMenu)

        
        self.lmpropCanvas = LmpropCanvas(self, nCellRows, nCellColumns,
                viewPosition=LmpropCanvas.xyPosition(0.3, 1.25))
        self.lmpropCanvas.grid(row=1, column=0)

    def onFileExport(self):
        fname = tkfd.asksaveasfilename(
            defaultextension=".ps",
            initialfile="lmprop.ps",
            filetypes=(
                ("Postscript", "*.ps" ),
                ("All",        "*" ),
            )
        )
        if fname:
            self.lmpropCanvas.postscript(file=fname)

    def setViewMode(self, viewMode):
        self.viewModeVar.set(viewMode)
        if viewMode == SCENE_ONLY_VIEW_MODE:
            self.lmpropCanvas.setRaytraceVisibility(False)
            self.lmpropCanvas.setVoxelizationVisibility(False)
            self.lmpropCanvas.setLightVoxelizationVisibility(False)
            self.lmpropCanvas.setViewVoxelizationVisibility(False)
        elif viewMode == RAYTRACE_VIEW_MODE:
            self.lmpropCanvas.setRaytraceVisibility(True)
            self.lmpropCanvas.setVoxelizationVisibility(False)
            self.lmpropCanvas.setLightVoxelizationVisibility(False)
            self.lmpropCanvas.setViewVoxelizationVisibility(False)
        elif viewMode == VOXELIZATION_VIEW_MODE:
            self.lmpropCanvas.setRaytraceVisibility(False)
            self.lmpropCanvas.setVoxelizationVisibility(True)
            self.lmpropCanvas.setLightVoxelizationVisibility(False)
            self.lmpropCanvas.setViewVoxelizationVisibility(False)
        elif viewMode == LIGHT_VOXELIZATION_VIEW_MODE:
            self.lmpropCanvas.setRaytraceVisibility(False)
            self.lmpropCanvas.setVoxelizationVisibility(True)
            self.lmpropCanvas.setLightVoxelizationVisibility(True)
            self.lmpropCanvas.setViewVoxelizationVisibility(False)
        else:
            assert viewMode == VIEW_VOXELIZATION_VIEW_MODE
            self.lmpropCanvas.setRaytraceVisibility(False)
            self.lmpropCanvas.setVoxelizationVisibility(True)
            self.lmpropCanvas.setLightVoxelizationVisibility(False)
            self.lmpropCanvas.setViewVoxelizationVisibility(True)

    def onViewMode(self):
        self.setViewMode(self.viewModeVar.get())


class ApplicationCanvas(tk.Canvas):
    """
    a wrapper around tk.Canvas

    This fixes some problems that Tkinter has with things like Numpy arrays
    and adds arrows.
    """
    def __init__(self, parent, **kwargs):
        super().__init__(parent, **kwargs)

    def drawArrow(self, p0, p1, **kwargs):
        if 'width' not in kwargs:
            kwargs['width'] = LINE_WIDTH
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
        if 'width' not in kwargs:
            kwargs['width'] = LINE_WIDTH
        return self.create_oval((pCtr[0] - r, pCtr[1] - r),
                                (pCtr[0] + r, pCtr[1] + r), **kwargs)

    def drawImage(self, p, image, **kwargs):
        return self.create_image(p[0], p[1], image=image, **kwargs)

    def drawLine(self, p0, p1, **kwargs):
        if 'width' not in kwargs:
            kwargs['width'] = LINE_WIDTH
        return self.create_line(p0[0], p0[1], p1[0], p1[1], **kwargs)

    def drawPolygon(self, polygon, **kwargs):
        if 'width' not in kwargs:
            kwargs['width'] = LINE_WIDTH
        # convert arrays, if need be
        return self.create_polygon([ tuple(p) for p in polygon ],
                                   **kwargs)

    def drawRectangle(self, pUL, pLR, **kwargs):
        if 'width' not in kwargs:
            kwargs['width'] = LINE_WIDTH
        return self.create_rectangle(pUL[0], pUL[1], pLR[0], pLR[1], **kwargs)

    def drawText(self, p, text, **kwargs):
        return self.create_text(tuple(p), text=text, **kwargs)


class Wall(Mesh):

    def __init__(self, cell, k, p0, p1):
        super().__init__(p0, p1, Cell.WALL_MESH_RESOLUTION)
        self.cell = cell
        d = p1 - p0
        self.normal = normalize(np.array((d[1], -d[0])))


class Cell:

    SIZE = 300 # width and height of (square) cells
    WALL_MESH_RESOLUTION = 9
    INSET = LINE_WIDTH  # distinguish cell walls

    def __init__(self, lmpropCanvas, i, j, **kwargs):
        self.lmpropCanvas = lmpropCanvas
        insetVector = np.array((Cell.INSET, Cell.INSET))
        pUL = lmpropCanvas.xyPosition(j,   i  ) + insetVector
        pLR = lmpropCanvas.xyPosition(j+1, i+1) - insetVector
        self.walls = []
        self.polygon = np.array(
            (pUL,
             (pUL[0], pLR[1]),
             pLR,
             (pLR[0], pUL[1])))
        for k in range(4):
            wall = Wall(self, k, self.polygon[k], self.polygon[(k+1) % 4])
            self.walls.append(wall)

        self.backgroundId = self.lmpropCanvas.drawRectangle(pUL, pLR, **kwargs)

        popupMenu = tk.Menu(self.lmpropCanvas, tearoff=0)
        popupMenu.add_checkbutton(label="Toggle Highlight",
                                  onvalue=True, offvalue=False,
                                  command=self.onToggleHighlight)
        popupMenu.add_checkbutton(label="Show Propagation",
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

    CELL_ARRAY_OFFSET = 30 # allow for row and column labels
    LIGHT_RADIUS = 10
    cellAtRowColumn = {}

    def __init__(self, parent, nCellRows, nCellColumns, viewPosition):
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

        self.lightPosition = self.xyPosition(0.2, 0.3)

        center = self.xyPosition(1.2, 0.8)
        radius = 20
        self.blockingPolygon = createRegularPolygon(center, radius, 6)

        illuminatedCell = LmpropCanvas.cellAtRowColumn[
            self.nCellRows-1, self.nCellColumns-1]
        self.drawLightRayPropagation(self.lightPosition, illuminatedCell)

        self.drawPolygon(self.blockingPolygon, fill=BLOCKING_POLYGON_COLOR,
                         outline='black', width=1)
        self.drawText(center + np.array((radius + 40, radius + 5)),
                      "blocking\nobject",
                      font=ITEM_LABEL_FONT, justify=tk.CENTER)

        fovDeg = 30
        viewDirection = self.xyDisplacement(0.6, 0.03)
        viewMesh = self.getViewMesh(viewPosition, viewDirection, fovDeg)

        self.drawViewRayPropagation(viewPosition, viewMesh, illuminatedCell)

        center = self.xyPosition(2.4, 1.5)
        radius = 80
        self.targetPolygon = createRegularPolygon(center, radius, 
                                                  10, phaseDeg=10)
        self.drawText(center + np.array((radius + 40, 40 - radius)),
                      "target\nobject",
                      font=ITEM_LABEL_FONT, justify=tk.CENTER)
        self.drawPolygon(self.targetPolygon,
                         fill=TARGET_POLYGON_COLOR, outline='black', width=1)

        self.drawRaytrace(viewPosition, viewDirection, fovDeg=fovDeg)

        self.drawCircle(self.lightPosition,
                        LmpropCanvas.LIGHT_RADIUS, fill=LIGHT_COLOR)
        self.drawText(self.lightPosition - np.array((0, 40)),
                      "point\nlight source",
                      font=ITEM_LABEL_FONT, justify=tk.CENTER)

        # need to save image reference so it won't be garbage-collected
        self.image = tk.PhotoImage(file=EYE_FNAME)
        self.drawText(viewPosition - np.array((0, 30)),
                      "viewer",
                      font=ITEM_LABEL_FONT, justify=tk.CENTER)
        self.drawImage(viewPosition, self.image)

        # Initially, only show the scene.
        self.setVoxelizationVisibility(False)
        self.setLightVoxelizationVisibility(False)
        self.setRaytraceVisibility(False)
        self.setViewVoxelizationVisibility(False)
            
    def drawCells(self):
        cells = []
        for i in range(self.nCellRows):
            for j in range(self.nCellColumns):
                cell = Cell(self, i, j, tag=VOXELIZATION_TAG)
                self.cellAtRowColumn[i, j] = cell
                cells.append(cell)
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

    def drawLightRayPropagation(self, position, illuminatedCell):
        for cell in self.cells:
            if cell == illuminatedCell:
                continue # don't propagate rays through the final cell
            for wall in illuminatedCell.walls:
                (p, n) = lineThroughPoints(wall.p0, wall.p1)
                if np.dot(p - position, n) > 0:
                    self.drawPropagationToCellWall(position, cell, wall,
                                                   LIGHT_VOXELIZATION_TAG)

    def drawViewRayPropagation(self, position, viewMesh, illuminatedCell):
        for cell in self.cells:
            if cell == illuminatedCell:
                continue # don't propagate rays through the illuminated cell
            self.drawPropagationToCellWall(position, cell, viewMesh,
                                           VIEW_VOXELIZATION_TAG)

    def drawRayPropagationInCell(self, ray, cell, tag):
        # Compute the intersection of `ray` with two (presumably)
        # of `cell`'s walls.
        intersections = ray.intersectsPolygon(cell.polygon)
        if intersections:
            blockingIntersections = ray.intersectsPolygon(
                    self.blockingPolygon)
            if len(intersections) == 1:
                if blockingIntersections:
                    if blockingIntersections[0].t <= intersections[0].t:
                        self.drawArrow(intersections[0].p,
                                       blockingIntersections[0].p,
                                       tag=tag)
                    if blockingIntersections[1].t <= intersections[0].t:
                        self.drawLine(blockingIntersections[1].p,
                                      intersections[0].p,
                                      fill=BLOCKED_RAY_COLOR,
                                      tag=tag)
                    if intersections[0].t <= blockingIntersections[0].t:
                        self.drawArrow(ray.o,
                                       intersections[0].p,
                                       tag=tag)
                else:
                    self.drawArrow(ray.o,
                                   intersections[0].p,
                                   tag=tag)
            else:
                assert len(intersections) == 2
                if blockingIntersections:
                    if (intersections[0].t <= blockingIntersections[0].t
                             <= intersections[1].t):
                        self.drawArrow(intersections[0].p,
                                       blockingIntersections[0].p,
                                       tag=tag)
                    if (intersections[0].t <= blockingIntersections[1].t
                             <= intersections[1].t):
                        self.drawLine(blockingIntersections[1].p,
                                      intersections[1].p,
                                      fill=BLOCKED_RAY_COLOR,
                                      tag=tag)
                    if blockingIntersections[1].t <= intersections[1].t:
                        self.drawLine(blockingIntersections[1].p,
                                      intersections[1].p,
                                      fill=BLOCKED_RAY_COLOR,
                                      tag=tag)
                else:
                    self.drawArrow(intersections[0].p,
                                   intersections[1].p,
                                   tag=tag)

    def drawPropagationToCellWall(self, position, cell, destinationWall, tag):
        for destinationPoint in destinationWall.samplePoints:
            o = position
            d = normalize(destinationPoint - o)
            ray = Ray(o, d)
            self.drawRayPropagationInCell(ray, cell, tag)

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

    def getViewMesh(self, viewPosition, viewDirection, fovDeg):
        d = mag(viewDirection)
        imageCenter = viewPosition + viewDirection
        normalizedViewDirection = normalize(viewDirection)

        # perpendicular to the view direction
        normalizedViewDirectionPerp = normalize(
            np.array((viewDirection[1], -viewDirection[0])))

        fov = fovDeg * np.pi / 180
        halfW = d * np.tan(fov/2)

        p0 = imageCenter + halfW * normalizedViewDirectionPerp
        p1 = imageCenter - halfW * normalizedViewDirectionPerp
        return Mesh(p0, p1, 7)

    def drawRaytrace(self, viewPosition, viewDirection, fovDeg):
        viewMesh = self.getViewMesh(viewPosition, viewDirection, fovDeg)
        imageCenter = viewPosition + viewDirection

        nSamples = len(viewMesh.samplePoints)
        w = mag(viewMesh.p1 - viewMesh.p0)
        pixelSize = w / (nSamples - 1)

        for samplePoint in viewMesh.samplePoints:
            pixelPolygon = createOrientedSquare(samplePoint, pixelSize,
                                                viewDirection)
            viewRay = Ray(viewPosition, samplePoint - viewPosition)
            viewRayIntersections = viewRay.intersectsPolygon(
                    self.targetPolygon)
            if viewRayIntersections:
                self.drawArrow(viewPosition, viewRayIntersections[0].p,
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

    def setRaytraceVisibility(self, isVisible):
        if isVisible:
            self.itemconfigure(RAYTRACE_TAG, state=tk.NORMAL)
        else:
            self.itemconfigure(RAYTRACE_TAG, state=tk.HIDDEN)

    def setLightVoxelizationVisibility(self, isVisible):
        if isVisible:
            self.itemconfigure(LIGHT_VOXELIZATION_TAG, state=tk.NORMAL)
        else:
            self.itemconfigure(LIGHT_VOXELIZATION_TAG, state=tk.HIDDEN)

    def setViewVoxelizationVisibility(self, isVisible):
        if isVisible:
            self.itemconfigure(VIEW_VOXELIZATION_TAG, state=tk.NORMAL)
        else:
            self.itemconfigure(VIEW_VOXELIZATION_TAG, state=tk.HIDDEN)

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

    application = Application(root, root,
                              nCellRows = 2,
                              nCellColumns = 3)
    application.grid(row=0, column=0)

    root.mainloop()


if __name__ == '__main__':
    main()
