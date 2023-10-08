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

ss = 0
for objfile in glob.glob(r'D:/temp_need/paper/classic_region_soma/*.obj'):
    objfile = objfile.replace('\\','/')
    s = objfile.split('/')[-1]
    objname = s.replace('.','').replace('-','')
    objdisplayname = objname + 'Display'

    r = (int(s.split('_')[0].split('g')[0].split('r')[-1])+1)/256.
    g = (int(s.split('_')[0].split('b')[0].split('g')[-1])+1)/256.
    b = (int(s.split('_')[0].split('b')[-1])+1)/256.
    print(f'object name: {objname},  r/g/b: {r}/{g}/{b}')

    find_str = f"{objname} = WavefrontOBJReader(registrationName='{s}',FileName='{objfile}')"
    print(f'{find_str}')# create a new 'Wavefront OBJ Reader'
    exec(find_str)

    active_str = f'SetActiveSource({objname})'
    print(active_str)# set active source # SetActiveSource(type89obj)
    exec(active_str)

    display_str = f'{objdisplayname} = Show({objname}, renderView1,"GeometryRepresentation")'
    print(display_str)# get display properties # type89objDisplay = GetDisplayProperties(type89obj, view=renderView1)
    exec(display_str)

    # trace defaults for the display properties.
    property_str = f"{objdisplayname}.Representation = 'Surface';\n{objdisplayname}.ColorArrayName = [None, ''];\n{objdisplayname}.SelectTCoordArray = 'None';\n{objdisplayname}.SelectNormalArray = 'None';\n{objdisplayname}.SelectTangentArray = 'None';\n{objdisplayname}.OSPRayScaleFunction = 'PiecewiseFunction';\n{objdisplayname}.SelectOrientationVectors = 'None';\n{objdisplayname}.ScaleFactor = 52.587;\n{objdisplayname}.SelectScaleArray = 'None';\n{objdisplayname}.GlyphType = 'Arrow';\n{objdisplayname}.GlyphTableIndexArray = 'None';\n{objdisplayname}.GaussianRadius = 2.62935;\n{objdisplayname}.SetScaleArray = [None, ''];\n{objdisplayname}.ScaleTransferFunction = 'PiecewiseFunction';\n{objdisplayname}.OpacityArray = [None, ''];\n{objdisplayname}.OpacityTransferFunction = 'PiecewiseFunction';\n{objdisplayname}.DataAxesGrid = 'GridAxesRepresentation';\n{objdisplayname}.PolarAxes = 'PolarAxesRepresentation';\n{objdisplayname}.SelectInputVectors = [None, ''];\n{objdisplayname}.WriteLog = ''"
    print(property_str)
    exec(property_str)

    # r = max(1-ss*0.002,0)
    # g = min(ss*0.001,1)
    # b = min(0.5+ss*0.001,1)
    color_str = f'{objdisplayname}.AmbientColor = [{r}, {g}, {b}];\n{objdisplayname}.DiffuseColor = [{r}, {g}, {b}]'
    print(color_str)# change solid color
    exec(color_str)

    opacity_str = f'{objdisplayname}.Opacity = 0.2'
    print(opacity_str)# Properties modified on total_filterobjDisplay
    exec(opacity_str)

    ss += 1
    if ss%10 == 0:
        layout1 = GetLayout()
        layout1.SetSize(1512, 838)
        renderView1.CameraPosition = [160, -1800, 264]
        renderView1.CameraFocalPoint = [264, 160, 228]
        renderView1.CameraViewUp = [0.0, 0.0, 1.0]
        renderView1.CameraViewAngle = 15.822784810126583
        renderView1.CameraParallelScale = 360
        SaveScreenshot(f'C:/Users/12561/Desktop/paper/paper-fig1/{ss}_soma.png', renderView1, ImageResolution=[1512, 838])

layout1 = GetLayout()
layout1.SetSize(1512, 838)
renderView1.CameraPosition = [160, -1800, 264]
renderView1.CameraFocalPoint = [264, 160, 228]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.CameraViewAngle = 15.822784810126583
renderView1.CameraParallelScale = 360
SaveScreenshot(f'C:/Users/12561/Desktop/paper/paper-fig1/soma_above.png', renderView1, ImageResolution=[1512, 838])

renderView1.ResetActiveCameraToPositiveZ()
renderView1.ResetCamera(False)
renderView1.AdjustRoll(-180.0)
layout1.SetSize(1512, 838)
renderView1.CameraPosition = [264, 160, -1800]
renderView1.CameraFocalPoint = [264, 160, 228]
renderView1.CameraViewUp = [0.0, -1.0, 0.0]
renderView1.CameraViewAngle = 15.822784810126583
renderView1.CameraParallelScale = 360
SaveScreenshot(f'C:/Users/12561/Desktop/paper/paper-fig1/soman_front.png', renderView1, ImageResolution=[1512, 838])

renderView1.ResetActiveCameraToNegativeX()
renderView1.ResetCamera(False)
renderView1.AdjustRoll(-90.0)
layout1.SetSize(1512, 838)
renderView1.CameraPosition = [-1800, 160, 264]
renderView1.CameraFocalPoint = [264, 160, 228]
renderView1.CameraViewUp = [1.0, 0.0, 0.0]
renderView1.CameraViewAngle = 15.822784810126583
renderView1.CameraParallelScale = 360
SaveScreenshot(f'C:/Users/12561/Desktop/paper/paper-fig1/soman_back.png', renderView1, ImageResolution=[1512, 838])
