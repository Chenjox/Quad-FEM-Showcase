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
testvtu = XMLUnstructuredGridReader(registrationName=datasetIn + '.vtu', FileName=[datasetIn + '.vtu'])
testvtu.PointArrayStatus = ['Displacement', 'Stresses']

# Properties modified on testvtu
testvtu.TimeArray = 'None'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
testvtuDisplay = Show(testvtu, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
testvtuDisplay.Representation = 'Surface'
testvtuDisplay.ColorArrayName = [None, '']
testvtuDisplay.SelectTCoordArray = 'None'
testvtuDisplay.SelectNormalArray = 'None'
testvtuDisplay.SelectTangentArray = 'None'
testvtuDisplay.OSPRayScaleArray = 'Displacement'
testvtuDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
testvtuDisplay.SelectOrientationVectors = 'Displacement'
testvtuDisplay.ScaleFactor = 0.1
testvtuDisplay.SelectScaleArray = 'Displacement'
testvtuDisplay.GlyphType = 'Arrow'
testvtuDisplay.GlyphTableIndexArray = 'Displacement'
testvtuDisplay.GaussianRadius = 0.005
testvtuDisplay.SetScaleArray = ['POINTS', 'Displacement']
testvtuDisplay.ScaleTransferFunction = 'PiecewiseFunction'
testvtuDisplay.OpacityArray = ['POINTS', 'Displacement']
testvtuDisplay.OpacityTransferFunction = 'PiecewiseFunction'
testvtuDisplay.DataAxesGrid = 'GridAxesRepresentation'
testvtuDisplay.PolarAxes = 'PolarAxesRepresentation'
testvtuDisplay.ScalarOpacityUnitDistance = 0.8908987181403394
testvtuDisplay.OpacityArrayName = ['POINTS', 'Displacement']
testvtuDisplay.SelectInputVectors = ['POINTS', 'Displacement']
testvtuDisplay.WriteLog = ''

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
testvtuDisplay.ScaleTransferFunction.Points = [-8.806071473698821e-16, 0.0, 0.5, 0.0, 0.5415384615384604, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
testvtuDisplay.OpacityTransferFunction.Points = [-8.806071473698821e-16, 0.0, 0.5, 0.0, 0.5415384615384604, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

#changing interaction mode based on data extents
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# change representation type
testvtuDisplay.SetRepresentationType('Surface With Edges')

# set scalar coloring
ColorBy(testvtuDisplay, ('POINTS', 'Displacement', 'Magnitude'))

# rescale color and/or opacity maps used to include current data range
testvtuDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
testvtuDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Displacement'
displacementLUT = GetColorTransferFunction('Displacement')

# get opacity transfer function/opacity map for 'Displacement'
displacementPWF = GetOpacityTransferFunction('Displacement')

# get 2D transfer function for 'Displacement'
displacementTF2D = GetTransferFunction2D('Displacement')

# set scalar coloring
ColorBy(testvtuDisplay, ('POINTS', 'Displacement', 'X'))

# rescale color and/or opacity maps used to exactly fit the current data range
testvtuDisplay.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(displacementLUT, testvtuDisplay)

# Properties modified on displacementLUT
displacementLUT.NumberOfTableValues = 10

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(800, 520)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# save screenshot
SaveScreenshot(datasetIn + '-ux.png', renderView1, ImageResolution=[800, 520])

# set scalar coloring
ColorBy(testvtuDisplay, ('POINTS', 'Displacement', 'Y'))

# rescale color and/or opacity maps used to exactly fit the current data range
testvtuDisplay.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(displacementLUT, testvtuDisplay)

# layout/tab size in pixels
layout1.SetSize(800, 520)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# save screenshot
SaveScreenshot(datasetIn + '-uy.png', renderView1, ImageResolution=[800, 520])

# set scalar coloring
ColorBy(testvtuDisplay, ('POINTS', 'Stresses', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(displacementLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
testvtuDisplay.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
testvtuDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Stresses'
stressesLUT = GetColorTransferFunction('Stresses')

# get opacity transfer function/opacity map for 'Stresses'
stressesPWF = GetOpacityTransferFunction('Stresses')

# get 2D transfer function for 'Stresses'
stressesTF2D = GetTransferFunction2D('Stresses')

# set scalar coloring
ColorBy(testvtuDisplay, ('POINTS', 'Stresses', 'X'))

# rescale color and/or opacity maps used to exactly fit the current data range
testvtuDisplay.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(stressesLUT, testvtuDisplay)

# layout/tab size in pixels
layout1.SetSize(800, 520)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# save screenshot
SaveScreenshot(datasetIn + '-sigmax.png', renderView1, ImageResolution=[800, 520])

# set scalar coloring
ColorBy(testvtuDisplay, ('POINTS', 'Stresses', 'Y'))

# rescale color and/or opacity maps used to exactly fit the current data range
testvtuDisplay.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(stressesLUT, testvtuDisplay)

# layout/tab size in pixels
layout1.SetSize(800, 520)

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

# save screenshot
SaveScreenshot(datasetIn + '-sigmay.png', renderView1, ImageResolution=[800, 520])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(800, 520)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.InteractionMode = '2D'
renderView1.CameraPosition = [0.5, 0.5, 10000.0]
renderView1.CameraFocalPoint = [0.5, 0.5, 0.0]
renderView1.CameraParallelScale = 0.7071067811865476

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).