# trace generated using paraview version 5.11.0-RC1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 11

import os,glob
from paraview.simple import *

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get the material library
materialLibrary1 = GetMaterialLibrary()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
renderView1.UseColorPaletteForBackground = 0
renderView1.Background = [1.0, 1.0, 1.0]

# reset view to fit data
renderView1.ResetActiveCameraToPositiveY()
renderView1.ResetCamera(False)

#vtk_des
vtk_des = input()

ss = 0
for vtkfile in glob.glob(vtk_des + '/*.vtk'):
    v = os.path.basename(vtkfile)
    vtkname = v.replace('.','_').replace('-','')
    vtkdisplayname = vtkname + 'Display'
    ss += 1

    r = (int(v.split('_')[0].split('g')[0].split('r')[-1])+1)/256.
    g = (int(v.split('_')[0].split('b')[0].split('g')[-1])+1)/256.
    b = (int(v.split('_')[0].split('b')[-1])+1)/256.
    region = v.replace(v.split('_')[0],'').replace('.vtk','')[1:]
    print(f'\nobject name: {vtkname},  r/g/b/region: {r}/{g}/{b}/{region}')

    find_str = f"{vtkname} = LegacyVTKReader(registrationName='{v}',FileNames='{vtkfile}')"
    print(f'{find_str}')# create a new 'Wavefront OBJ Reader'
    exec(find_str)

    active_str = f'SetActiveSource({vtkname})'
    print(active_str)# set active source # SetActiveSource(type89obj)
    exec(active_str)

    display_str = f'{vtkdisplayname} = Show({vtkname}, renderView1,"GeometryRepresentation")'
    print(display_str)# get display properties # type89objDisplay = GetDisplayProperties(type89obj, view=renderView1)
    exec(display_str)

    # get color transfer function/color map for 'scalars'
    scalarsLUT = GetColorTransferFunction('scalars')

    property_str = f"{vtkdisplayname}.Representation = 'Surface';\n{vtkdisplayname}.ColorArrayName = ['POINTS', 'scalars'];\n{vtkdisplayname}.LookupTable = scalarsLUT;\n{vtkdisplayname}.SelectTCoordArray = 'None';\n{vtkdisplayname}.SelectNormalArray = 'Normals';\n{vtkdisplayname}.SelectTangentArray = 'None';\n{vtkdisplayname}.OSPRayScaleArray = 'scalars';\n{vtkdisplayname}.OSPRayScaleFunction = 'PiecewiseFunction';\n{vtkdisplayname}.SelectOrientationVectors = 'None';\n{vtkdisplayname}.ScaleFactor = 52.499898004531865;\n{vtkdisplayname}.SelectScaleArray = 'scalars';\n{vtkdisplayname}.GlyphType = 'Arrow';\n{vtkdisplayname}.GlyphTableIndexArray = 'scalars';\n{vtkdisplayname}.GaussianRadius = 2.624994900226593;\n{vtkdisplayname}.SetScaleArray = ['POINTS', 'scalars'];\n{vtkdisplayname}.ScaleTransferFunction = 'PiecewiseFunction';\n{vtkdisplayname}.OpacityArray = ['POINTS', 'scalars'];\n{vtkdisplayname}.OpacityTransferFunction = 'PiecewiseFunction';\n{vtkdisplayname}.DataAxesGrid = 'GridAxesRepresentation';\n{vtkdisplayname}.PolarAxes = 'PolarAxesRepresentation';\n{vtkdisplayname}.SelectInputVectors = ['POINTS', 'Normals'];\n{vtkdisplayname}.WriteLog = '';\n{vtkdisplayname}.ScaleTransferFunction.Points = [255.0, 0.0, 0.5, 0.0, 255.03125, 1.0, 0.5, 0.0];\n{vtkdisplayname}.OpacityTransferFunction.Points = [255.0, 0.0, 0.5, 0.0, 255.03125, 1.0, 0.5, 0.0];\n{vtkdisplayname}.SetScalarBarVisibility(renderView1, True)"
    print(property_str)# trace defaults for the display properties.
    exec(property_str)

    # reset view to fit data
    renderView1.ResetCamera(False)
    # get opacity transfer function/opacity map for 'scalars'
    scalarsPWF = GetOpacityTransferFunction('scalars')
    # get 2D transfer function for 'scalars'
    scalarsTF2D = GetTransferFunction2D('scalars')

    scalar_str = f"ColorBy({vtkdisplayname}, None);\nHideScalarBarIfNotNeeded(scalarsLUT, renderView1)"
    print(scalar_str)# turn off scalar coloring# Hide the scalar bar for this color map if no visible data is colored by it.
    exec(scalar_str)    

    opacity_str = f'{vtkdisplayname}.Opacity = 0.7'
    print(opacity_str)# Properties modified on total_filterobjDisplay
    exec(opacity_str)

    color_str = f'{vtkdisplayname}.AmbientColor = [{r}, {g}, {b}];\n{vtkdisplayname}.DiffuseColor = [{r}, {g}, {b}]'
    print(color_str)# change solid color
    exec(color_str)

    if ss%10==0:
        layout1 = GetLayout()
        layout1.SetSize(1512, 838)
        renderView1.CameraPosition = [-1800, 160, 264]
        renderView1.CameraFocalPoint = [264, 160, 228]
        renderView1.CameraViewUp = [1.0, 0.0, 0.0]
        renderView1.CameraViewAngle = 15.822784810126583
        renderView1.CameraParallelScale = 360
        SaveScreenshot(f'C:/Users/12561/Desktop/paper/paper-fig1/vtk_{ss}.png', renderView1, ImageResolution=[1512, 838])

layout1 = GetLayout()
layout1.SetSize(1512, 838)
renderView1.CameraPosition = [160, -1800, 264]
renderView1.CameraFocalPoint = [264, 160, 228]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraViewAngle = 15.822784810126583
renderView1.CameraParallelScale = 360
SaveScreenshot(f'C:/Users/12561/Desktop/paper/paper-fig1/vtk_above.png', renderView1, ImageResolution=[1512, 838])

renderView1.ResetActiveCameraToPositiveZ()
renderView1.ResetCamera(False)
renderView1.AdjustRoll(-180.0)
layout1.SetSize(1512, 838)
renderView1.CameraPosition = [264, 160, -1800]
renderView1.CameraFocalPoint = [264, 160, 228]
renderView1.CameraViewUp = [0.0, -1.0, 0.0]
renderView1.CameraViewAngle = 15.822784810126583
renderView1.CameraParallelScale = 360
SaveScreenshot(f'C:/Users/12561/Desktop/paper/paper-fig1/vtk_front.png', renderView1, ImageResolution=[1512, 838])

renderView1.ResetActiveCameraToNegativeX()
renderView1.ResetCamera(False)
renderView1.AdjustRoll(-90.0)
layout1.SetSize(1512, 838)
renderView1.CameraPosition = [-1800, 160, 264]
renderView1.CameraFocalPoint = [264, 160, 228]
renderView1.CameraViewUp = [1.0, 0.0, 0.0]
renderView1.CameraViewAngle = 15.822784810126583
renderView1.CameraParallelScale = 360
SaveScreenshot(f'C:/Users/12561/Desktop/paper/paper-fig1/vtk_back.png', renderView1, ImageResolution=[1512, 838])
