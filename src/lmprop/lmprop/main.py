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
        self.showRaytrace = tk.BooleanVar(value=True) # initially
        viewMenu.add_checkbutton(label="Raytrace",
                                 onvalue=True,
                                 offvalue=False,
                                 variable=self.showRaytrace,
                                 command=self.onViewRaytrace)
        menuBar.add_cascade(label="View", menu=viewMenu)

        self.canvas = VoxelizationCanvas(self, nCellRows, nCellColumns)
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

    def onViewRaytrace(self):
        self.canvas.setRaytraceVisibility(self.showRaytrace.get())

    def onViewVoxelization(self):
        self.canvas.setVoxelizationVisibility(self.showVoxelization.get())


class ApplicationCanvas(tk.Canvas):
    """
    a wrapper around tk.Canvas

    This fixes some problems that Tkinter has with things like Numpy arrays
    and adds arrows.
    """
    CELL_ARRAY_OFFSET = 30 # allow for labels

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
        return self.create_text(p, text=text, **kwargs)


class Cell:

    SIZE = 300 # width and height of (square) cells
    WALL_MESH_RESOLUTION = 9
    INSET = 1  # distinguish cell walls

    def __init__(self, canvas, i, j, **kwargs):
        self.pUL = np.array((
            j * Cell.SIZE + Cell.INSET + ApplicationCanvas.CELL_ARRAY_OFFSET,
            i * Cell.SIZE + Cell.INSET + ApplicationCanvas.CELL_ARRAY_OFFSET))
        self.pLR = np.array((
            (j+1)*Cell.SIZE - Cell.INSET + ApplicationCanvas.CELL_ARRAY_OFFSET,
            (i+1)*Cell.SIZE - Cell.INSET + ApplicationCanvas.CELL_ARRAY_OFFSET)
            )
        self.canvas = canvas
        self.walls = []
        for (x0, y0, x1, y1) in (
                (self.pUL[0], self.pUL[1], self.pUL[0], self.pLR[1]),
                (self.pUL[0], self.pLR[1], self.pLR[0], self.pLR[1]),
                (self.pLR[0], self.pLR[1], self.pLR[0], self.pUL[1]),
                (self.pLR[0], self.pUL[1], self.pUL[0], self.pUL[1])):
            p0 = np.array((x0, y0))
            p1 = np.array((x1, y1))
            wall = Mesh(p0, p1, Cell.WALL_MESH_RESOLUTION)
            self.walls.append(wall)

        pUL = np.array((self.pUL[0], self.pUL[1]))
        pLR = np.array((self.pLR[0], self.pLR[1]))

        self.backgroundId = self.canvas.drawRectangle(pUL, pLR, **kwargs)

        popupMenu = tk.Menu(self.canvas, tearoff=0)
        popupMenu.add_checkbutton(label="Toggle Highlight",
                                  onvalue=True, offvalue=False,
                                  command=self.onToggleHighlight)
        def onPopupBackground(evt):
            try:
                popupMenu.tk_popup(evt.x_root, evt.y_root, 0)
            finally:
                popupMenu.grab_release()

        self.canvas.tag_bind(self.backgroundId, '<Button-3>',
                             onPopupBackground)

        self.highlight = False
        self.setHighlight()

    def drawBlocking(self, pLight, wallIndex, blockingPolygon, wallTo, 
                     **kwargs):
        intercedingWall = self.walls[wallIndex]
        for i in range(Cell.WALL_MESH_RESOLUTION):
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

    def onToggleHighlight(self):
        self.highlight = not self.highlight
        self.setHighlight()

    def setHighlight(self):
        self.canvas.itemconfigure(self.backgroundId,
                                  fill= 'white' if self.highlight
                                                else 'gray80')


class VoxelizationCanvas(ApplicationCanvas):
    """an ApplicationCanvas specficially intended for this application
    """
    def __init__(self, parent, nCellRows, nCellColumns):
        self.nCellRows = nCellRows
        self.nCellColumns= nCellColumns
        self.width = (self.nCellColumns * Cell.SIZE
                      + ApplicationCanvas.CELL_ARRAY_OFFSET)
        self.height = (self.nCellRows * Cell.SIZE
                      + ApplicationCanvas.CELL_ARRAY_OFFSET)
        super().__init__(parent, 
                         width=self.width, height=self.height,
                         background='white')

        self.corners = (
            (Cell.INSET, Cell.INSET),
            (Cell.INSET, self.height - Cell.INSET),
            (self.width - Cell.INSET, self.height - Cell.INSET),
            (self.width - Cell.INSET, Cell.INSET))

        self.cells = self.drawCells()

        pLight = self.xyPosition(0.3, 0.4)
        illuminatedCell = self.cells[self.nCellRows-1, self.nCellColumns-1]
        illuminatedWall = illuminatedCell.walls[0]

        blockingPolygon = createRegularPolygon(self.xyPosition(1.2, 0.8),
                                               20, 6)
        self.drawPolygon(blockingPolygon, fill='cyan', outline='black',
                         width=1)

        blockingCell = self.cells[0, 1]
        blockingCell.drawBlocking(pLight, 0, blockingPolygon, illuminatedWall,
                                  tag=VOXELIZATION_TAG)

        targetPolygon = createRegularPolygon(self.xyPosition(2.4, 1.5), 80, 10,
                                             phaseDeg=10)
        self.drawPolygon(targetPolygon, fill='blue', outline='black', width=1)

        eyePoint=self.xyPosition(0.3, 1.25)
        self.drawRaytrace(eyePoint,
                          viewDirection=self.xyDisplacement(0.6, 0.03),
                          pLight=pLight,
                          fovDeg=30,
                          targetPolygon=targetPolygon,
                          blockingPolygon=blockingPolygon)
        self.drawCircle(pLight, 10, fill='yellow')

        # need to save image reference so it won't be garbage-collected
        self.image = tk.PhotoImage(file=EYE_FNAME)
        self.drawImage(eyePoint, self.image)

        self.setVoxelizationVisibility(True)
        self.setRaytraceVisibility(True)
            
                    
    def drawCells(self):
        cells = {}
        for i in range(self.nCellRows):
            for j in range(self.nCellColumns):
                cells[i,j] = Cell(self, i, j,
                            tag=VOXELIZATION_TAG)
        for j in range(self.nCellColumns):
            self.drawText(((j + 0.5) * Cell.SIZE + ApplicationCanvas.CELL_ARRAY_OFFSET,
                           0.5 * ApplicationCanvas.CELL_ARRAY_OFFSET),
                          "{}".format(j), font=("Times", 16),
                          tag=VOXELIZATION_TAG)
        for i in range(self.nCellColumns):
            self.drawText((0.5 * ApplicationCanvas.CELL_ARRAY_OFFSET,
                          (i + 0.5) * Cell.SIZE + ApplicationCanvas.CELL_ARRAY_OFFSET),
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

    def drawRaytrace(self, eyePoint, viewDirection, pLight, fovDeg,
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
        viewPlane = Mesh(p0, p1, 7)

        nSamples = len(viewPlane.samplePoints)
        w = 2 * halfW
        pixelWH = w / (nSamples - 1)

        for samplePoint in viewPlane.samplePoints:
            pixelPolygon = createOrientedSquare(samplePoint, pixelWH,
                                                viewDirection)
            lightRay = Ray(eyePoint, samplePoint - eyePoint)
            lightRayIntersections = lightRay.intersectsPolygon(targetPolygon)
            if lightRayIntersections:
                self.drawArrow(eyePoint, lightRayIntersections[0].p,
                               tag=RAYTRACE_TAG)
                self.drawPolygon(pixelPolygon, fill='blue',
                                 outline='black', width=1, tag=RAYTRACE_TAG)
                blockingRay = Ray(lightRayIntersections[0].p,
                                  pLight - lightRayIntersections[0].p)
                blockingRayIntersections = \
                        blockingRay.intersectsPolygon(blockingPolygon)
                if blockingRayIntersections:
                    self.drawArrow(lightRayIntersections[0].p, 
                                   blockingRayIntersections[0].p,
                                   fill='gray60',
                                   tag=RAYTRACE_TAG)
                else:
                    self.drawArrow(lightRayIntersections[0].p, pLight,
                                   tag=RAYTRACE_TAG)
            else:
                self.drawRay(lightRay, tag=RAYTRACE_TAG)
                self.drawPolygon(pixelPolygon, fill='white',
                                 outline='black', width=1, tag=RAYTRACE_TAG)

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

    def xyPosition(self, x, y):
        return np.array((x * Cell.SIZE + ApplicationCanvas.CELL_ARRAY_OFFSET,
                         y * Cell.SIZE + ApplicationCanvas.CELL_ARRAY_OFFSET),
                         dtype=float)

    def xyDisplacement(self, dx, dy):
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
