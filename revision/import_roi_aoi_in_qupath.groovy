import qupath.imagej.tools.ROIConverterIJ
import ij.plugin.filter.ThresholdToSelection
import ij.process.ByteProcessor
import ij.process.ImageProcessor
import ij.gui.ShapeRoi
import ij.gui.Roi
import javafx.application.Platform
import javafx.embed.swing.JFXPanel
import javafx.geometry.Insets
import javafx.scene.Scene
import javafx.scene.control.Button
import javafx.scene.layout.BorderPane
import javafx.stage.Stage
import loci.common.RandomAccessInputStream
import loci.common.services.ServiceFactory
import loci.formats.services.OMEXMLService
import ome.codecs.ZlibCodec
import ome.units.UNITS
import ome.units.quantity.Length
import ome.xml.meta.OMEXMLMetadata
import ome.xml.model.enums.Compression
import ome.xml.model.primitives.Color
import ome.xml.model.primitives.NonNegativeInteger
import qupath.imagej.tools.ROIConverterIJ
import qupath.imagej.tools.IJTools
import qupath.lib.common.ColorTools
import qupath.lib.geom.Point2
import qupath.lib.gui.dialogs.Dialogs
import qupath.lib.gui.dialogs.ParameterPanelFX
import qupath.lib.gui.prefs.PathPrefs
import qupath.lib.gui.scripting.QPEx
import qupath.lib.objects.PathAnnotationObject
import qupath.lib.objects.PathCellObject
import qupath.lib.objects.PathDetectionObject
import qupath.lib.objects.PathROIObject
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.plugins.parameters.ParameterChangeListener
import qupath.lib.plugins.parameters.ParameterList
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.*

// Read the OME-XML file and count the number of ROIs
def omeXml = getCurrentImageData().getServer().dumpMetadata() as String
factory = new ServiceFactory()
service = factory.getInstance(OMEXMLService.class)
omeXmlMetadata = service.createOMEXMLMetadata(omeXml)
roiCount = omeXmlMetadata.getROICount()
if (roiCount < 1) {
    print "No ROIs found to import!"
    return
} else {
    print ("Found " + roiCount + " ROIs in the OME-XML file!\n")
}

// Definitions
nameIndexes = new HashMap<String, Integer>();
newPathObjects = []
thinLineStrokeWidths = new HashSet<>()
thickLineStrokeWidths = new HashSet<>()
pathClasses = new HashSet<>()

// Unpack a NonNegativeInteger from OME-XML, using a default value of 0 if null
int getValue(NonNegativeInteger v) {
    return v == null ? 0 : v.numberValue.intValue()
}

void setPathClassAndStroke(PathROIObject path, String className, Color color, Length strokeWidth) {
    def qpColor = null
    if (color != null) {
        qpColor = ColorTools.makeRGBA(color.red, color.green, color.blue, color.alpha)
    }

    if (strokeWidth != null) {
        strokeWidthValue = strokeWidth.value(UNITS.PIXEL)
        switch (path) {
            case PathDetectionObject:
                thinLineStrokeWidths.add(strokeWidthValue)
                break
            case PathAnnotationObject:
                thickLineStrokeWidths.add(strokeWidthValue)
                break
        }
    }

    if (className == null) {
        // No class to set, so just set the color directly
        if (qpColor != null) {
            path.setColor(qpColor)
        }
    } else {
        // set the class on the object
        def qpClass = PathClassFactory.getPathClass(className, qpColor)
        path.setPathClass(qpClass)

        // update list of unique classes so that the UI can be updated later
        pathClasses.add(qpClass);
    }
}

// Extracts the name from a given Shape
String getShapeName(PathROIObject path, int roiIdx, int shapeIdx, String className) {
    def shapeType = omeXmlMetadata.getShapeType(roiIdx, shapeIdx)
    baseRoiName = "ROI #" + roiIdx
    switch (shapeType) {

        case "Ellipse":
            roiName = baseRoiName
            break

        case "Polygon":
            roiName = baseRoiName
            break

        case "Mask":
            maskText = omeXmlMetadata.getMaskText(roiIdx, shapeIdx)
            roiName = baseRoiName + " - " + "AOI: " + maskText
            break

        case "Label":
            labelText = omeXmlMetadata.getLabelText(roiIdx, shapeIdx)
            roiName = baseRoiName + " - " + "Label: " + labelText
            break

        default:
            roiName = ""
    }
    return roiName
}

// Extracts the color from a given Shape
ome.xml.model.primitives.Color getShapeColor(PathROIObject path, int roiIdx, int shapeIdx, String className) {
    def shapeType = omeXmlMetadata.getShapeType(roiIdx, shapeIdx)
    def roiColor = new ome.xml.model.primitives.Color(255, 0, 0, 0)
    if (shapeType=="Mask") {
        roiColor = omeXmlMetadata.getMaskFillColor(roiIdx, shapeIdx)
        } else {
        null
    }
    return roiColor
}

// Converts the OME-XML shape with the given indexes to a QuPath ROI
qupath.lib.roi.interfaces.ROI importShape(PathROIObject path, int roiIdx, int shapeIdx, String className) {
    def shapeType = omeXmlMetadata.getShapeType(roiIdx, shapeIdx)
    println(String.format("(Annotation %d) Shape %d has type '%s'", roiIdx, shapeIdx, shapeType))
    def locked = null
    def roi = null

    switch (shapeType) {
        case "Ellipse":
            def color = omeXmlMetadata.getEllipseStrokeColor(roiIdx, shapeIdx)
            def strokeWidth = omeXmlMetadata.getEllipseStrokeWidth(roiIdx, shapeIdx)
            setPathClassAndStroke(path, className, color, strokeWidth)

            locked = omeXmlMetadata.getEllipseLocked(roiIdx, shapeIdx)

            def c = getValue(omeXmlMetadata.getEllipseTheC(roiIdx, shapeIdx))
            def z = getValue(omeXmlMetadata.getEllipseTheZ(roiIdx, shapeIdx))
            def t = getValue(omeXmlMetadata.getEllipseTheT(roiIdx, shapeIdx))
            def plane = new ImagePlane(c, z, t)
            def centroidX = omeXmlMetadata.getEllipseX(roiIdx, shapeIdx)
            def centroidY = omeXmlMetadata.getEllipseY(roiIdx, shapeIdx)
            def radiusX = omeXmlMetadata.getEllipseRadiusX(roiIdx, shapeIdx)
            def radiusY = omeXmlMetadata.getEllipseRadiusY(roiIdx, shapeIdx)
            def x = centroidX - radiusX
            def y = centroidY - radiusY
            def width = radiusX * 2
            def height = radiusY * 2
            roi = new EllipseROI(x, y, width, height, plane)
            break

        case "Line":
            def color = omeXmlMetadata.getLineStrokeColor(roiIdx, shapeIdx)
            def strokeWidth = omeXmlMetadata.getLineStrokeWidth(roiIdx, shapeIdx)
            setPathClassAndStroke(path, className, color, strokeWidth)

            locked = omexml.getLineLocked(roiIdx, shapeIdx)

            def c = getValue(omeXmlMetadata.getLineTheC(roiIdx, shapeIdx))
            def z = getValue(omeXmlMetadata.getLineTheZ(roiIdx, shapeIdx))
            def t = getValue(omeXmlMetadata.getLineTheT(roiIdx, shapeIdx))
            def plane = new ImagePlane(c, z, t)
            def x = omeXmlMetadata.getLineX1(roiIdx, shapeIdx)
            def y = omeXmlMetadata.getLineY1(roiIdx, shapeIdx)
            def x2 = omeXmlMetadata.getLineX2(roiIdx, shapeIdx)
            def y2 = omeXmlMetadata.getLineY2(roiIdx, shapeIdx)
            roi = new LineROI(x, y, x2, y2, plane)
            break

        case "Point":
            def color = omeXmlMetadata.getPointStrokeColor(roiIdx, shapeIdx)
            def strokeWidth = omeXmlMetadata.getPointStrokeWidth(roiIdx, shapeIdx)
            setPathClassAndStroke(path, className, color, strokeWidth)

            locked = omeXmlMetadata.getPointLocked(roiIdx, shapeIdx)

            def c = getValue(omeXmlMetadata.getPointTheC(roiIdx, shapeIdx))
            def z = getValue(omeXmlMetadata.getPointTheZ(roiIdx, shapeIdx))
            def t = getValue(omeXmlMetadata.getPointTheT(roiIdx, shapeIdx))
            def plane = new ImagePlane(c, z, t)
            def x = omeXmlMetadata.getPointX(roiIdx, shapeIdx)
            def y = omeXmlMetadata.getPointY(roiIdx, shapeIdx)
            roi = new PointsROI(x, y, plane)
            break

        case "Polygon":
            def color = omeXmlMetadata.getPolygonStrokeColor(roiIdx, shapeIdx)
            def strokeWidth = omeXmlMetadata.getPolygonStrokeWidth(roiIdx, shapeIdx)
            setPathClassAndStroke(path, className, color, strokeWidth)

            locked = omeXmlMetadata.getPolygonLocked(roiIdx, shapeIdx)

            def c = getValue(omeXmlMetadata.getPolygonTheC(roiIdx, shapeIdx))
            def z = getValue(omeXmlMetadata.getPolygonTheZ(roiIdx, shapeIdx))
            def t = getValue(omeXmlMetadata.getPolygonTheT(roiIdx, shapeIdx))
            def plane = new ImagePlane(c, z, t)
            def pointsString = omeXmlMetadata.getPolygonPoints(roiIdx, shapeIdx)
            def points = pointsString.split(/ /).collect { point ->
                def (x, y) = point.split(/,/)
                new Point2(x.toDouble(), y.toDouble())
            }
            roi = new PolygonROI(points, plane)
            break

        case "Polyline":
            def color = omeXmlMetadata.getPolylineStrokeColor(roiIdx, shapeIdx)
            def strokeWidth = omeXmlMetadata.getPolylineStrokeWidth(roiIdx, shapeIdx)
            setPathClassAndStroke(path, className, color, strokeWidth)

            locked = omeXmlMetadata.getPolylineLocked(roiIdx, shapeIdx)

            def c = getValue(omeXmlMetadata.getPolylineTheC(roiIdx, shapeIdx))
            def z = getValue(omeXmlMetadata.getPolylineTheZ(roiIdx, shapeIdx))
            def t = getValue(omeXmlMetadata.getPolylineTheT(roiIdx, shapeIdx))
            def plane = new ImagePlane(c, z, t)
            def pointsString = omeXmlMetadata.getPolylinePoints(roiIdx, shapeIdx)
            def points = pointsString.split(/ /).collect { point ->
                def (x, y) = point.split(/,/)
                new Point2(x.toDouble(), y.toDouble())
            }
            roi = new PolylineROI(points, plane)
            break

        case "Rectangle":
            def color = omeXmlMetadata.getRectangleStrokeColor(roiIdx, shapeIdx)
            def strokeWidth = omeXmlMetadata.getRectangleStrokeWidth(roiIdx, shapeIdx)
            setPathClassAndStroke(path, className, color, strokeWidth)

            locked = omeXmlMetadata.getRectangleLocked(roiIdx, shapeIdx)

            def c = getValue(omeXmlMetadata.getRectangleTheC(roiIdx, shapeIdx))
            def z = getValue(omeXmlMetadata.getRectangleTheZ(roiIdx, shapeIdx))
            def t = getValue(omeXmlMetadata.getRectangleTheT(roiIdx, shapeIdx))
            def plane = new ImagePlane(c, z, t)
            def x = omexml.getRectangleX(roiIdx, shapeIdx)
            def y = omexml.getRectangleY(roiIdx, shapeIdx)
            def width = omeXmlMetadata.getRectangleWidth(roiIdx, shapeIdx)
            def height = omeXmlMetadata.getRectangleHeight(roiIdx, shapeIdx)
            roi = new RectangleROI(x, y, width, height, plane)
            break

        case "Label":
            // convert it into a rectangle
            def color = omeXmlMetadata.getLabelStrokeColor(roiIdx, shapeIdx)
            def strokeWidth = omeXmlMetadata.getLabelStrokeWidth(roiIdx, shapeIdx)
            setPathClassAndStroke(path, className, color, strokeWidth)

            locked = omeXmlMetadata.getLabelLocked(roiIdx, shapeIdx)

            def c = getValue(omeXmlMetadata.getLabelTheC(roiIdx, shapeIdx))
            def z = getValue(omeXmlMetadata.getLabelTheZ(roiIdx, shapeIdx))
            def t = getValue(omeXmlMetadata.getLabelTheT(roiIdx, shapeIdx))
            def plane = new ImagePlane(c, z, t)
            def x = omeXmlMetadata.getLabelX(roiIdx, shapeIdx)
            def y = omeXmlMetadata.getLabelY(roiIdx, shapeIdx)
            def width = 10
            def height = 10
            roi = new RectangleROI(x, y, width, height, plane)
            break

        case "Mask":
            def x = omeXmlMetadata.getMaskX(roiIdx, shapeIdx)
            def y = omeXmlMetadata.getMaskY(roiIdx, shapeIdx)
            def width = omeXmlMetadata.getMaskWidth(roiIdx, shapeIdx).intValue()
            def height = omeXmlMetadata.getMaskHeight(roiIdx, shapeIdx).intValue()
            def binData = omeXmlMetadata.getMaskBinData(roiIdx, shapeIdx)
            def compression = omeXmlMetadata.getMaskBinDataCompression(roiIdx, shapeIdx)

            bits = binData
            if (compression == Compression.ZLIB) {
                bits = new ZlibCodec().decompress(bits, null)
            }
            def stream = new RandomAccessInputStream(bits)
            def bytes = new byte[width * height]
            (0..(bytes.length - 1)).each { bitIndex ->
                bytes[bitIndex] = stream.readBits(1) * Byte.MAX_VALUE;
            }

            def bp = new ByteProcessor(width, height, bytes)
            bp.setThreshold(Byte.MAX_VALUE - 1, 255, ImageProcessor.NO_LUT_UPDATE)
            def ijROI = new ThresholdToSelection().convert(bp)

            if (ijROI != null) {
                def c = getValue(omeXmlMetadata.getMaskTheC(roiIdx, shapeIdx))
                def z = getValue(omeXmlMetadata.getMaskTheZ(roiIdx, shapeIdx))
                def t = getValue(omeXmlMetadata.getMaskTheT(roiIdx, shapeIdx))

                if (ijROI instanceof ShapeRoi) {
                    roi = ROIConverterIJ.convertToAreaROI(ijROI, -1 * x, -1 * y, 1, c, z, t)
                } else {
                    roi = ROIConverterIJ.convertToPolygonOrAreaROI(ijROI, -1 * x, -1 * y, 1, c, z, t)
                }
            }
            break

        default:
            throw new Exception(String.format("(Annotation %d) Shape %d is of unknown type: %s", roiIdx, shapeIdx, shapeType))
    }

    if (locked != null) {
        path.setLocked(locked)
    }

    return roi
}


Integer getIndex(String roiName) {
    def index = nameIndexes.get(roiName);
    if (index != null) {
        nameIndexes.put(roiName, index + 1);
        return index;
    }
    nameIndexes.put(roiName, 2);
    return 1;
}


void chooseLineWidths() {

    new JFXPanel()
    if (!Platform.isFxApplicationThread()) {
        Platform.runLater({ chooseLineWidths() })
        return
    }

    def frame = new Stage()
    frame.setTitle("Line Thickness")
    def lineParams = new ParameterList().addTitleParameter("Choose Line Thickness")
            .addEmptyParameter("More than one value for line thickness was imported.")
            .addEmptyParameter("Please choose which value you'd like to use:")
    if (thickLineStrokeWidths.size() > 1) {
        def widths = thickLineStrokeWidths.toArray().sort()
        lineParams.addChoiceParameter("annotation", "Annotation Thickness", widths[0] as Number,
                widths as ArrayList<Number>, "Numeric choice")
    }
    if (thinLineStrokeWidths.size() > 1) {
        def widths = thinLineStrokeWidths.toArray().sort()
        lineParams.addChoiceParameter("detection", "Detection Thickness", widths[0] as Number,
                widths as ArrayList<Number>, "Numeric choice")
    }

    def borderPane = new BorderPane()
    def panel = new ParameterPanelFX(lineParams)

    panel.addParameterChangeListener(new ParameterChangeListener() {
        @Override
        void parameterChanged(ParameterList params, String key, boolean isAdjusting) {
            def param = params.getParameters().get(key)
            if (key == "annotation") {
                PathPrefs.setThickStrokeThickness(param.getValue() as float)
            }
            if (key == "detection") {
                PathPrefs.setThinStrokeThickness(param.getValue() as float)
            }
        }
    })

    def button = new Button("OK")
    button.setDefaultButton(true)
    button.setOnAction({ frame.close() })
    borderPane.setCenter(panel.getPane())
    borderPane.setBottom(button)
    borderPane.setPadding(new Insets(10, 10, 10, 10))

    frame.setScene(new Scene(borderPane))
    frame.show()
}

if (thickLineStrokeWidths.size() > 1 || thinLineStrokeWidths.size() > 1) {
    chooseLineWidths()
}


// Run the extraction of ROI and AOI
for (int roiIdx = 0;roiIdx<roiCount;roiIdx++) {
    def mapAnnotations = [:]
    def annotationCount = omeXmlMetadata.getROIAnnotationRefCount(roiIdx)
    shapeCount = omeXmlMetadata.getShapeCount(roiIdx)
    println(String.format("Annotation %d contains %d shapes:", roiIdx, shapeCount))

    (0..(shapeCount - 1)).each { shapeIdx ->

        PathROIObject path
        path = new PathAnnotationObject()
        def className = mapAnnotations["qupath:class"]
        if (className == null) {
            className = mapAnnotations["class"]
            if (className == null) {
                // If there is no explicit class name, but the ROI has a name, use that as class name
                className = omeXmlMetadata.getROIName(roiIdx)
            }
        }

        roi = importShape(path, roiIdx, shapeIdx, className)
        roiName = getShapeName(path, roiIdx, shapeIdx, className)
        roiColor = getShapeColor(path, roiIdx, shapeIdx, className)

        if (roi != null) {
            path.setName(roiName)
            path.setROI(roi)
            path.setColor(roiColor.getRed(), roiColor.getGreen(), roiColor.getBlue())
            mapAnnotations.keySet().each {
                if (it.toString().startsWith("qupath:metadata:")) {
                    path.storeMetadataValue(it.toString().replaceFirst(/qupath:metadata:/, ""),
                            mapAnnotations[it].toString())
                }
            }

            // store the second shape as the cell nucleus
            if (path.isCell() && effectiveCount == 1) {
                roiIdx++
                nucleusROI = importShape(path, roiIdx, shapeIdx, className)

                // no setter for the nucleus ROI
                path = new PathCellObject(roi, nucleusROI, path.getPathClass())
            }

            newPathObjects.add(path)
        }
    }
    println(String.format("Done for Annotation %d!\n", roiIdx))
}

QPEx.getCurrentHierarchy().addObjects(newPathObjects)
