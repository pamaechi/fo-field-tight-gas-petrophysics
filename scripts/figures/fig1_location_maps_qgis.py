# -*- coding: utf-8 -*-
"""
F-O field location maps (Figure 1 of Paper 1 and Figure 1 of Paper 2).

One run exports both figures. They share the inset and block map but differ in panel (c):
    Paper 1: wells labelled with the depth of the reservoir top (IP picks, m TVDSS);
    Paper 2: wells labelled with the data used in that paper
             (P = wireline pretests, D = drill-stem test, C = core capillary pressure).

Run inside QGIS (3.28 LTR or newer):
    Plugins > Python Console > Show Editor > Open Script... > select this file > Run Script.

Panels
    (a) regional inset (Esri World Ocean Base basemap, no labels);
    (b) block map with all eight wells;
    (c) depth-structure map of 1At1/top USM, modified after PetroSA (2008; F-O8 well recommendation,
        Figure 4), georeferenced to UTM zone 34S with the well positions (residuals < 5 m), with the
        IP top-reservoir picks posted at each well.

Inputs (keep them in the same folder as this script, or edit the paths below)
    petrosa_1At1_TUSM_depth_map_UTM34S.tif   georeferenced, cleaned copy of the PetroSA map

Edit only the USER SETTINGS block.
"""
import os
from qgis.core import (
    QgsProject, QgsVectorLayer, QgsRasterLayer, QgsField, QgsFeature, QgsGeometry, QgsPointXY,
    QgsCoordinateReferenceSystem, QgsMarkerSymbol, QgsLineSymbol, QgsFillSymbol,
    QgsPalLayerSettings, QgsTextFormat, QgsTextBufferSettings, QgsVectorLayerSimpleLabeling,
    QgsCategorizedSymbolRenderer, QgsRendererCategory, QgsSingleSymbolRenderer,
    QgsPrintLayout, QgsLayoutItemMap, QgsLayoutItemLegend, QgsLayoutItemScaleBar,
    QgsLayoutItemPicture, QgsLayoutItemLabel, QgsLayoutPoint, QgsLayoutSize, QgsUnitTypes,
    QgsLayoutExporter, QgsRectangle, QgsLayoutItemMapGrid, QgsLayoutItemMapOverview,
    QgsApplication, QgsLegendStyle
)
from qgis.PyQt.QtCore import QVariant
from qgis.PyQt.QtGui import QColor, QFont

# ============================== USER SETTINGS ==============================
try:
    HERE = os.path.dirname(os.path.abspath(__file__))
except NameError:                                  # QGIS console editor
    HERE = os.path.expanduser("~/FO_map")
OUT_DIR = os.path.join(HERE, "output")
DEPTH_MAP = os.path.join(HERE, "petrosa_1At1_TUSM_depth_map_UTM34S.tif")
USE_BASEMAP = True           # inset basemap (Esri World Ocean Base, no labels); needs internet.
# Credit in the caption: "Basemap: Esri, GEBCO, NOAA, National Geographic, Garmin, HERE and other contributors".
WELL_SYMBOL_MM = 1.6         # well symbol diameter (mm)
LABEL_PT = 6.5               # label size (pt)
# ===========================================================================

os.makedirs(OUT_DIR, exist_ok=True)
project = QgsProject.instance()
CRS_WGS = QgsCoordinateReferenceSystem("EPSG:4326")
CRS_UTM = QgsCoordinateReferenceSystem("EPSG:32734")
MM = QgsUnitTypes.LayoutMillimeters

# Well positions: UTM zone 34S (WGS 84) from the Petrel well list (FO_FS_FR_wells.xlsx).
# Top reservoir: IP picks (Table_of_Tops.xlsx) = top of the uppermost reservoir package; F-R1 and F-S1: 1At1
# (Well_Tops_FS1_FR1.xlsx). Datum = drill-floor elevation above MSL (F-O6 includes the 1.45 m TVD loss).
#        well    easting      northing     top_MD   datum  group            top pick
WELLS = [("F-O1", 639328.94, 6109427.82, 3702.29, 26.0,  "F-O field",     "P4"),
         ("F-O2", 642033.26, 6112562.05, 3615.28, 25.0,  "F-O field",     "P3"),
         ("F-O3", 641979.34, 6110573.75, 3714.25, 22.0,  "F-O field",     "P2b"),
         ("F-O4", 640371.06, 6112941.81, 3709.83, 21.0,  "F-O field",     "P5"),
         ("F-O6", 640031.34, 6110776.00, 3708.00, 26.45, "F-O field",     "P5"),
         ("F-O8", 639399.18, 6109492.19, 3708.93, 25.0,  "F-O field",     "P3"),
         ("F-R1", 645729.52, 6108365.29, 3877.00, 30.0,  "Off structure", "1At1"),
         ("F-S1", 654429.32, 6108103.58, 3620.00, 22.0,  "Off structure", "1At1")]

# ---------------- text formats ----------------
def fmt(size, bold=False, buffer=True):
    f = QgsTextFormat(); fn = QFont("Arial"); fn.setPointSizeF(size); fn.setBold(bold)
    f.setFont(fn); f.setSize(size)
    if buffer:
        b = QgsTextBufferSettings(); b.setEnabled(True); b.setSize(0.7); b.setColor(QColor("white"))
        f.setBuffer(b)
    return f

# ---------------- well layers ----------------
wells = QgsVectorLayer("Point?crs=EPSG:32734", "Wells", "memory")
pr = wells.dataProvider()
pr.addAttributes([QgsField("well", QVariant.String), QgsField("group", QVariant.String),
                  QgsField("top_md", QVariant.Double), QgsField("top_tvdss", QVariant.Double),
                  QgsField("pick", QVariant.String), QgsField("label", QVariant.String),
                  QgsField("label_b", QVariant.String), QgsField("label_p2", QVariant.String)])
wells.updateFields()
# Paper 2 data codes: P = wireline pretests, D = drill-stem test, C = core capillary pressure
P2_DATA = {"F-O1": "P", "F-O2": "D, C", "F-O3": "P", "F-O4": "P, D", "F-O6": "C", "F-O8": "P",
           "F-R1": "water leg", "F-S1": ""}
for name, x, y, md, kb, grp, pick in WELLS:
    tvd = round(md - kb, 1)
    lab = f"{name} ({tvd:.0f} m)"
    lab_b = name
    lab_p2 = f"{name} ({P2_DATA[name]})" if P2_DATA[name] else name
    if name == "F-O1":                       # F-O1 and F-O8 surface locations are c. 95 m apart
        lab, lab_b = f"F-O1 ({tvd:.0f} m)\nF-O8 (3684 m)", "F-O1/F-O8"
        lab_p2 = f"F-O1 ({P2_DATA['F-O1']})\nF-O8 ({P2_DATA['F-O8']})"
    elif name == "F-O8":
        lab, lab_b, lab_p2 = "", "", ""
    f = QgsFeature(wells.fields()); f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
    f.setAttributes([name, grp, md, tvd, pick, lab, lab_b, lab_p2]); pr.addFeature(f)
wells.updateExtents()

def well_renderer():
    s1 = QgsMarkerSymbol.createSimple({"name": "circle", "color": "#d7191c", "outline_color": "black",
                                       "size": str(WELL_SYMBOL_MM), "outline_width": "0.25"})
    s2 = QgsMarkerSymbol.createSimple({"name": "circle", "color": "white", "outline_color": "black",
                                       "size": str(WELL_SYMBOL_MM), "outline_width": "0.4"})
    return QgsCategorizedSymbolRenderer("group", [
        QgsRendererCategory("F-O field", s1, "F-O field well"),
        QgsRendererCategory("Off structure", s2, "Well outside the F-O closure")])

def label_layer(layer, field, size):
    p = QgsPalLayerSettings(); p.fieldName = field
    p.placement = QgsPalLayerSettings.AroundPoint; p.dist = 0.8
    p.setFormat(fmt(size)); layer.setLabeling(QgsVectorLayerSimpleLabeling(p)); layer.setLabelsEnabled(True)

wells.setRenderer(well_renderer()); label_layer(wells, "label", LABEL_PT)
wells_b = wells.clone(); wells_b.setName("Wells (block map)")
wells_b.setRenderer(well_renderer()); label_layer(wells_b, "label_b", LABEL_PT)
wells_p2 = wells.clone(); wells_p2.setName("Wells (Paper 2)")
wells_p2.setRenderer(well_renderer()); label_layer(wells_p2, "label_p2", LABEL_PT)

# F-O field marker and towns for the inset
field_pt = QgsVectorLayer("Point?crs=EPSG:32734&field=name:string", "F-O field", "memory")
fp = QgsFeature(field_pt.fields()); fp.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(641000, 6110800)))
fp.setAttributes(["F-O field"]); field_pt.dataProvider().addFeature(fp)
field_pt.setRenderer(QgsSingleSymbolRenderer(QgsMarkerSymbol.createSimple(
    {"name": "star", "color": "#d7191c", "outline_color": "black", "size": "3.5"})))
label_layer(field_pt, "name", LABEL_PT)
towns = QgsVectorLayer("Point?crs=EPSG:4326&field=name:string", "Towns", "memory")
for nm, la, lo in [("Mossel Bay", -34.183, 22.146), ("Cape Agulhas", -34.833, 20.003),
                   ("Cape Town", -33.925, 18.424)]:
    ft = QgsFeature(towns.fields()); ft.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(lo, la)))
    ft.setAttributes([nm]); towns.dataProvider().addFeature(ft)
towns.setRenderer(QgsSingleSymbolRenderer(QgsMarkerSymbol.createSimple(
    {"name": "square", "color": "black", "size": "1.3"})))
label_layer(towns, "name", LABEL_PT)

# depth map and legend-only symbols for its content
dmap = QgsRasterLayer(DEPTH_MAP, "Depth map of 1At1/top USM (PetroSA, 2008)")
if not dmap.isValid():
    raise FileNotFoundError(f"Depth map not found: {DEPTH_MAP}")
def legend_line(name, props):
    l = QgsVectorLayer("LineString?crs=EPSG:32734", name, "memory")
    l.setRenderer(QgsSingleSymbolRenderer(QgsLineSymbol.createSimple(props))); return l
leg_contour = legend_line("Depth contour, 25 m interval (m TVDSS)", {"color": "#222222", "width": "0.2"})
leg_fault = legend_line("Fault", {"color": "black", "width": "1.0"})
leg_gwc = legend_line("Gas–water contact", {"color": "#d62728", "width": "0.4"})

layers = [wells, wells_b, wells_p2, field_pt, towns, dmap, leg_contour, leg_fault, leg_gwc]
osm = None
if USE_BASEMAP:
    osm = QgsRasterLayer("type=xyz&url=https://server.arcgisonline.com/ArcGIS/rest/services/Ocean/"
                         "World_Ocean_Base/MapServer/tile/{z}/{y}/{x}&zmax=10&zmin=0",
                         "Esri World Ocean Base", "wms")
    if osm.isValid():
        layers.append(osm)
    else:
        print("Basemap could not be loaded; the inset is drawn without it.")
        osm = None
for lyr in layers:
    project.addMapLayer(lyr)

# ---------------- layout ----------------
def build_figure(layout_name, out_stem, main_wells, code_note=None):
  lm = project.layoutManager()
  old = lm.layoutByName(layout_name)
  if old:
    lm.removeLayout(old)
  layout = QgsPrintLayout(project); layout.initializeDefaults(); layout.setName(layout_name)
  layout.pageCollection().page(0).setPageSize(QgsLayoutSize(180, 122, MM))   # double-column width

  def add_map(x, y, w, h, crs, extent, lyrs):
      m = QgsLayoutItemMap(layout); m.setCrs(crs)
      m.attemptMove(QgsLayoutPoint(x, y, MM)); m.attemptResize(QgsLayoutSize(w, h, MM))
      m.zoomToExtent(extent); m.setLayers(lyrs); m.setFrameEnabled(True)
      layout.addLayoutItem(m)
      return m

  def add_grid(m, interval):
      g = QgsLayoutItemMapGrid("Graticule", m); m.grids().addGrid(g)
      g.setCrs(CRS_WGS); g.setIntervalX(interval); g.setIntervalY(interval)
      g.setStyle(QgsLayoutItemMapGrid.FrameAnnotationsOnly)
      g.setFrameStyle(QgsLayoutItemMapGrid.Zebra); g.setFrameWidth(0.8)
      g.setAnnotationEnabled(True)
      g.setAnnotationFormat(QgsLayoutItemMapGrid.DegreeMinute); g.setAnnotationPrecision(0)
      g.setAnnotationFrameDistance(1.2)              # gap between zebra frame and coordinates
      try:
          g.setAnnotationTextFormat(fmt(5.5, buffer=False))
      except AttributeError:
          g.setAnnotationFont(QFont("Arial", 5))
      S = QgsLayoutItemMapGrid
      for side in (S.Left, S.Right, S.Bottom, S.Top):
          g.setAnnotationPosition(S.OutsideMapFrame, side)
      g.setAnnotationDisplay(S.HideAll, S.Top)
      g.setAnnotationDirection(S.Vertical, S.Left)
      g.setAnnotationDirection(S.Vertical, S.Right)
      g.setAnnotationDirection(S.Horizontal, S.Bottom)
      try:                                           # show only latitudes on the sides, longitudes below
          g.setAnnotationDisplay(S.LatitudeOnly, S.Left); g.setAnnotationDisplay(S.LatitudeOnly, S.Right)
          g.setAnnotationDisplay(S.LongitudeOnly, S.Bottom)
      except AttributeError:
          pass

  def add_label(text, x, y, size=7.5, bold=True):
      lb = QgsLayoutItemLabel(layout); lb.setText(text)
      try:
          lb.setTextFormat(fmt(size, bold, buffer=False))
      except AttributeError:
          f = QFont("Arial"); f.setPointSizeF(size); f.setBold(bold); lb.setFont(f)
      lb.attemptMove(QgsLayoutPoint(x, y, MM)); lb.attemptResize(QgsLayoutSize(8, 5, MM))
      layout.addLayoutItem(lb)

  def add_scalebar(m, x, y, seg_km):
      sb = QgsLayoutItemScaleBar(layout); sb.setLinkedMap(m); sb.setStyle("Single Box")
      sb.setUnits(QgsUnitTypes.DistanceKilometers); sb.setUnitLabel("km")
      sb.setNumberOfSegments(2); sb.setNumberOfSegmentsLeft(0); sb.setUnitsPerSegment(seg_km)
      sb.setHeight(1.2); sb.setLabelBarSpace(0.6); sb.setBoxContentSpace(0.8)   # label close to the bar
      try:
          sb.setTextFormat(fmt(5.5))
      except AttributeError:
          sb.setFont(QFont("Arial", 5))
      sb.setBackgroundEnabled(False)
      sb.attemptMove(QgsLayoutPoint(x, y, MM)); layout.addLayoutItem(sb)

  # (a) regional inset
  inset_layers = [field_pt, towns] + ([osm] if osm else [])
  inset = add_map(12, 5, 50, 36, CRS_WGS, QgsRectangle(18.0, -36.0, 24.5, -33.4), inset_layers)
  add_grid(inset, 2.0); add_label("(a)", 13, 5.5)

  # (b) block map with all wells
  bx = wells.extent(); bx.grow(2500)
  block = add_map(12, 49, 50, 40, CRS_UTM, bx, [wells_b])
  add_grid(block, 0.05); add_label("(b)", 13, 49.5)
  add_scalebar(block, 14, 82.5, 2)
  ov = QgsLayoutItemMapOverview("Block map extent", inset); ov.setLinkedMap(block)
  ov.setFrameSymbol(QgsFillSymbol.createSimple({"color": "0,0,0,0", "outline_color": "red", "outline_width": "0.5"}))
  inset.overviews().addOverview(ov)

  # (c) PetroSA depth map; well labels differ between the two papers
  cx = QgsRectangle(637400, 6107700, 646350, 6113900)          # F-O closure and F-R1
  main = add_map(74, 5, 96, 84, CRS_UTM, cx, [main_wells, dmap])
  add_grid(main, 0.02); add_label("(c)", 75, 5.5)
  add_scalebar(main, 76, 82.5, 1)
  ov2 = QgsLayoutItemMapOverview("Depth map extent", block); ov2.setLinkedMap(main)
  ov2.setFrameSymbol(QgsFillSymbol.createSimple({"color": "0,0,0,0", "outline_color": "#1f78b4", "outline_width": "0.5"}))
  block.overviews().addOverview(ov2)

  na = QgsLayoutItemPicture(layout)
  na.setPicturePath(os.path.join(QgsApplication.svgPaths()[0], "arrows", "NorthArrow_02.svg"))
  na.attemptMove(QgsLayoutPoint(161, 7, MM)); na.attemptResize(QgsLayoutSize(6, 9, MM))
  layout.addLayoutItem(na)

  # legend in a box
  leg = QgsLayoutItemLegend(layout); leg.setLinkedMap(main); leg.setTitle("")
  leg.setAutoUpdateModel(False)
  root = leg.model().rootGroup(); root.clear()
  for lyr in [wells, leg_contour, leg_fault, leg_gwc]:
      root.addLayer(lyr)
  try:                                               # hide the "Wells" group heading
      from qgis.core import QgsLegendRenderer
      QgsLegendRenderer.setNodeLegendStyle(root.findLayer(wells), QgsLegendStyle.Hidden)
  except Exception:
      pass
  leg.setColumnCount(3); leg.setSplitLayer(True); leg.setEqualColumnWidth(False)
  for st in (QgsLegendStyle.Group, QgsLegendStyle.Subgroup, QgsLegendStyle.SymbolLabel):
      try:
          f = leg.rstyle(st).textFormat(); f.setSize(6.5); leg.rstyle(st).setTextFormat(f)
      except AttributeError:
          leg.setStyleFont(st, QFont("Arial", 6))
  leg.setSymbolHeight(2.5); leg.setSymbolWidth(5); leg.setBoxSpace(2)
  leg.setFrameEnabled(True); leg.setBackgroundEnabled(True)
  leg.attemptMove(QgsLayoutPoint(12, 99, MM)); leg.setResizeToContents(True)
  layout.addLayoutItem(leg)
  if code_note:                                   # explanation of the data codes, under the legend
    nt = QgsLayoutItemLabel(layout); nt.setText(code_note)
    try:
      nt.setTextFormat(fmt(6.5, buffer=False))
    except AttributeError:
      nt.setFont(QFont("Arial", 6))
    nt.attemptMove(QgsLayoutPoint(12, 116.5, MM)); nt.attemptResize(QgsLayoutSize(160, 4.5, MM))
    layout.addLayoutItem(nt)
  lm.addLayout(layout)

  exp = QgsLayoutExporter(layout)
  ps = QgsLayoutExporter.ImageExportSettings(); ps.dpi = 300
  exp.exportToImage(os.path.join(OUT_DIR, out_stem + ".png"), ps)
  exp.exportToPdf(os.path.join(OUT_DIR, out_stem + ".pdf"), QgsLayoutExporter.PdfExportSettings())

build_figure("FO_Figure1_Paper1", "P1_Fig1_location_structure_map", wells)
build_figure("FO_Figure1_Paper2", "P2_Fig1_location_data_map", wells_p2,
             "Data used in this paper:  P = wireline pretests;  D = drill-stem test;  C = core capillary pressure.")
print("Done. Outputs in", OUT_DIR)
for f in wells.getFeatures():
    print(f["well"], f["pick"], f["top_md"], f["top_tvdss"])
