# trace generated using paraview version 5.11.0-RC2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import *

datasetIn = sys.argv[1]
print("Dataset = " + datasetIn)
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'XML Unstructured Grid Reader'
patchCXDispvtu = XMLUnstructuredGridReader(registrationName=datasetIn + '.vtu', FileName=[datasetIn + '.vtu'])
patchCXDispvtu.PointArrayStatus = ['Displacement', 'Stresses']

# Properties modified on patchCXDispvtu
patchCXDispvtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
patchCXDispvtuDisplay = Show(patchCXDispvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
patchCXDispvtuDisplay.Representation = 'Surface'
patchCXDispvtuDisplay.ColorArrayName = [None, '']
patchCXDispvtuDisplay.SelectTCoordArray = 'None'
patchCXDispvtuDisplay.SelectNormalArray = 'None'
patchCXDispvtuDisplay.SelectTangentArray = 'None'
patchCXDispvtuDisplay.OSPRayScaleArray = 'Displacement'
patchCXDispvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
patchCXDispvtuDisplay.SelectOrientationVectors = 'Displacement'
patchCXDispvtuDisplay.ScaleFactor = 0.15000000000000002
patchCXDispvtuDisplay.SelectScaleArray = 'Displacement'
patchCXDispvtuDisplay.GlyphType = 'Arrow'
patchCXDispvtuDisplay.GlyphTableIndexArray = 'Displacement'
patchCXDispvtuDisplay.GaussianRadius = 0.0075
patchCXDispvtuDisplay.SetScaleArray = ['POINTS', 'Displacement']
patchCXDispvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
patchCXDispvtuDisplay.OpacityArray = ['POINTS', 'Displacement']
patchCXDispvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
patchCXDispvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
patchCXDispvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
patchCXDispvtuDisplay.ScalarOpacityUnitDistance = 1.0542695885492728
patchCXDispvtuDisplay.OpacityArrayName = ['POINTS', 'Displacement']
patchCXDispvtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
patchCXDispvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
patchCXDispvtuDisplay.ScaleTransferFunction.Points = [0.999999999999999, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
patchCXDispvtuDisplay.OpacityTransferFunction.Points = [0.999999999999999, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.75, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.75, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# change representation type
patchCXDispvtuDisplay.SetRepresentationType('Surface With Edges')

# set scalar coloring
ColorBy(patchCXDispvtuDisplay, ('POINTS', 'Displacement', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
patchCXDispvtuDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
patchCXDispvtuDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Displacement'
displacementLUT = GetColorTransferFunction('Displacement')

# get opacity transfer function/opacity map for 'Displacement'
displacementPWF = GetOpacityTransferFunction('Displacement')

# get 2D transfer function for 'Displacement'
displacementTF2D = GetTransferFunction2D('Displacement')

# set scalar coloring
ColorBy(patchCXDispvtuDisplay, ('POINTS', 'Displacement', 'X'))

# rescale color and/or opacity maps used to exactly fit the current data range
patchCXDispvtuDisplay.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(displacementLUT, patchCXDispvtuDisplay)

# Properties modified on displacementLUT
displacementLUT.NumberOfTableValues = 10

# create a new 'Warp By Vector'
warpByVector1 = WarpByVector(registrationName='WarpByVector1', Input=patchCXDispvtu)
warpByVector1.Vectors = ['POINTS', 'Displacement']

# show data in view
warpByVector1Display = Show(warpByVector1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
warpByVector1Display.Representation = 'Surface'
warpByVector1Display.ColorArrayName = ['POINTS', 'Displacement']
warpByVector1Display.LookupTable = displacementLUT
warpByVector1Display.SelectTCoordArray = 'None'
warpByVector1Display.SelectNormalArray = 'None'
warpByVector1Display.SelectTangentArray = 'None'
warpByVector1Display.OSPRayScaleArray = 'Displacement'
warpByVector1Display.OSPRayScaleFunction = 'PiecewiseFunction'
warpByVector1Display.SelectOrientationVectors = 'Displacement'
warpByVector1Display.ScaleFactor = 0.15000000000000008
warpByVector1Display.SelectScaleArray = 'Displacement'
warpByVector1Display.GlyphType = 'Arrow'
warpByVector1Display.GlyphTableIndexArray = 'Displacement'
warpByVector1Display.GaussianRadius = 0.007500000000000003
warpByVector1Display.SetScaleArray = ['POINTS', 'Displacement']
warpByVector1Display.ScaleTransferFunction = 'PiecewiseFunction'
warpByVector1Display.OpacityArray = ['POINTS', 'Displacement']
warpByVector1Display.OpacityTransferFunction = 'PiecewiseFunction'
warpByVector1Display.DataAxesGrid = 'GridAxesRepresentation'
warpByVector1Display.PolarAxes = 'PolarAxesRepresentation'
warpByVector1Display.ScalarOpacityFunction = displacementPWF
warpByVector1Display.ScalarOpacityUnitDistance = 1.0542695885492734
warpByVector1Display.OpacityArrayName = ['POINTS', 'Displacement']
warpByVector1Display.SelectInputVectors = ['POINTS', 'Displacement']
warpByVector1Display.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
warpByVector1Display.ScaleTransferFunction.Points = [0.999999999999999, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
warpByVector1Display.OpacityTransferFunction.Points = [0.999999999999999, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

# hide data in view
Hide(patchCXDispvtu, renderView1)

# show color bar/color legend
warpByVector1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# change representation type
warpByVector1Display.SetRepresentationType('Surface With Edges')

# set active source
SetActiveSource(patchCXDispvtu)

# show data in view
patchCXDispvtuDisplay = Show(patchCXDispvtu, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
patchCXDispvtuDisplay.SetScalarBarVisibility(renderView1, True)

# set active source
SetActiveSource(warpByVector1)

# Properties modified on warpByVector1Display
warpByVector1Display.Opacity = 0.1

# Properties modified on warpByVector1Display
warpByVector1Display.Opacity = 1.0

# Properties modified on warpByVector1
warpByVector1.ScaleFactor = 0.1

# update the view to ensure updated data information
renderView1.Update()

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(974, 575)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.75, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.75, 0.0]
renderView1.CameraParallelScale = 0.9013878188659973

# save screenshot
SaveScreenshot(datasetIn + '-ux.png', renderView1, ImageResolution=[974, 575])
#SaveScreenshot('X:\Programmieren\Rust\Quad-FEM-Showcase\patchC-X-Disp-ux.png', renderView1, ImageResolution=[974, 575])

# set scalar coloring
ColorBy(warpByVector1Display, ('POINTS', 'Displacement', 'Y'))

# rescale color and/or opacity maps used to exactly fit the current data range
warpByVector1Display.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(displacementLUT, warpByVector1Display)

# layout/tab size in pixels
layout1.SetSize(974, 575)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.75, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.75, 0.0]
renderView1.CameraParallelScale = 0.9013878188659973

# save screenshot
SaveScreenshot(datasetIn + '-uy.png', renderView1, ImageResolution=[974, 575])

# set scalar coloring
ColorBy(warpByVector1Display, ('POINTS', 'Stresses', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(displacementLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
warpByVector1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
warpByVector1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Stresses'
stressesLUT = GetColorTransferFunction('Stresses')

# get opacity transfer function/opacity map for 'Stresses'
stressesPWF = GetOpacityTransferFunction('Stresses')

# get 2D transfer function for 'Stresses'
stressesTF2D = GetTransferFunction2D('Stresses')

# set scalar coloring
ColorBy(warpByVector1Display, ('POINTS', 'Stresses', 'X'))

# rescale color and/or opacity maps used to exactly fit the current data range
warpByVector1Display.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(stressesLUT, warpByVector1Display)

# hide data in view
Hide(warpByVector1, renderView1)

# set active source
SetActiveSource(patchCXDispvtu)

# set scalar coloring
ColorBy(patchCXDispvtuDisplay, ('POINTS', 'Stresses', 'X'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(displacementLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
patchCXDispvtuDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
patchCXDispvtuDisplay.SetScalarBarVisibility(renderView1, True)

# Properties modified on stressesLUT
stressesLUT.NumberOfTableValues = 10

# layout/tab size in pixels
layout1.SetSize(974, 575)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.75, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.75, 0.0]
renderView1.CameraParallelScale = 0.9013878188659973

# save screenshot
SaveScreenshot(datasetIn + '-sigmax.png', renderView1, ImageResolution=[974, 575])
#SaveScreenshot('X:\Programmieren\Rust\Quad-FEM-Showcase\patchC-X-Disp-sigmax.png', renderView1, ImageResolution=[974, 575])

# set scalar coloring
ColorBy(patchCXDispvtuDisplay, ('POINTS', 'Stresses', 'Y'))

# rescale color and/or opacity maps used to exactly fit the current data range
patchCXDispvtuDisplay.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(stressesLUT, patchCXDispvtuDisplay)

# layout/tab size in pixels
layout1.SetSize(974, 575)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.75, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.75, 0.0]
renderView1.CameraParallelScale = 0.9013878188659973

# save screenshot
SaveScreenshot(datasetIn + '-sigmay.png', renderView1, ImageResolution=[974, 575])
#SaveScreenshot('X:\Programmieren\Rust\Quad-FEM-Showcase\patchC-X-Disp-sigmay.png', renderView1, ImageResolution=[974, 575])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(974, 575)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.75, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.75, 0.0]
renderView1.CameraParallelScale = 0.9013878188659973

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).