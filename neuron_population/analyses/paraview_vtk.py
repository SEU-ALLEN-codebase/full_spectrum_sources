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


rgblist = [[112.0, 193.0, 222.0], [188.0, 160.0, 215.0], 
[246.0, 179.0, 255.0], [227.0, 106.0, 158.0], 
[227.0, 200.0, 94.0], [126.0, 233.0, 119.0], 
[223.0, 64.0, 20.0],[42.0, 145.0, 255.0],
[140.0, 255.0, 40.0],[60.0, 255.0, 0.0],
[255.0, 0.0, 255.0],[125.0, 0.0, 242.0],
[55.0, 0.0, 55.0],[0.0, 0.0, 0.0],
]


vtk_des = 'F://platform-paper//vtk//'
savefig = 'C://Users//user//Desktop//plots//'

# process 1-5(3)
para = 'targetmodules_cs'
start,end = [61,120]

# process 6
# para = 'targetmodules_is'
# start,end = [-1,1000]

#process 7
# para = 'intermodules_is'
# start,end = [-1,1000]

#process 8
# para = 'intermodules_cs'
# start,end = [-1,1000]

#process (9)
# para = 'interregionsets_cs'
# start,end = [-1,1000]

#process 10
# para = 'interregionsets_cs'
# start,end = [-1,1000]

modulefile = f'F://platform-paper//fig2modules_figs3neurites_1224//middle_files//{para}.txt'
with open(modulefile,'r') as f:
    line = f.readlines()[0]
    if line[0]=='\"':
        line=line[1:-1]
    v_strs = eval(line)


flag = False
for iv_str in v_strs:
    if para[:6]=='target':
        iv,tv,v_str,ss,s = iv_str
        if (iv>end or iv<start): 
            continue
    else:
        iv,v_str,ss,s = iv_str
        tv = iv

    figfile = f'{savefig}//{para}_{tv}.jpeg'
    if os.path.exists(figfile):
        continue

    print(iv,':',tv,v_str)

    ss = 0
    for v in v_str:
        vtkname = v.replace('-','')
        vtkfile = vtk_des + v +'.vtk'
        vtkdisplayname = vtkname + 'Display'

        r,g,b = [rgb/256. for rgb in rgblist[ss%14]]
        print(f'{ss,} object name: {vtkname},  r/g/b: {r}/{g}/{b}')

        try:
            find_str = f"{vtkname} = LegacyVTKReader(registrationName='{v}',FileNames='{vtkfile}')"
            #print(f'{find_str}')# create a new 'Wavefront OBJ Reader'
            exec(find_str)
        except:
            print(f'Error: do not find\n{find_str}')
            continue
        
        active_str = f'SetActiveSource({vtkname})'
        #print(active_str)# set active source # SetActiveSource(type89obj)
        exec(active_str)

        display_str = f'{vtkdisplayname} = Show({vtkname}, renderView1,"GeometryRepresentation")'
        #print(display_str)# get display properties # type89objDisplay = GetDisplayProperties(type89obj, view=renderView1)
        exec(display_str)

        # get color transfer function/color map for 'scalars'
        scalarsLUT = GetColorTransferFunction('scalars')

        property_str = f"{vtkdisplayname}.Representation = 'Surface';\n{vtkdisplayname}.ColorArrayName = ['POINTS', 'scalars'];\n{vtkdisplayname}.LookupTable = scalarsLUT;\n{vtkdisplayname}.SelectTCoordArray = 'None';\n{vtkdisplayname}.SelectNormalArray = 'Normals';\n{vtkdisplayname}.SelectTangentArray = 'None';\n{vtkdisplayname}.OSPRayScaleArray = 'scalars';\n{vtkdisplayname}.OSPRayScaleFunction = 'PiecewiseFunction';\n{vtkdisplayname}.SelectOrientationVectors = 'None';\n{vtkdisplayname}.ScaleFactor = 52.499898004531865;\n{vtkdisplayname}.SelectScaleArray = 'scalars';\n{vtkdisplayname}.GlyphType = 'Arrow';\n{vtkdisplayname}.GlyphTableIndexArray = 'scalars';\n{vtkdisplayname}.GaussianRadius = 2.624994900226593;\n{vtkdisplayname}.SetScaleArray = ['POINTS', 'scalars'];\n{vtkdisplayname}.ScaleTransferFunction = 'PiecewiseFunction';\n{vtkdisplayname}.OpacityArray = ['POINTS', 'scalars'];\n{vtkdisplayname}.OpacityTransferFunction = 'PiecewiseFunction';\n{vtkdisplayname}.DataAxesGrid = 'GridAxesRepresentation';\n{vtkdisplayname}.PolarAxes = 'PolarAxesRepresentation';\n{vtkdisplayname}.SelectInputVectors = ['POINTS', 'Normals'];\n{vtkdisplayname}.WriteLog = '';\n{vtkdisplayname}.ScaleTransferFunction.Points = [255.0, 0.0, 0.5, 0.0, 255.03125, 1.0, 0.5, 0.0];\n{vtkdisplayname}.OpacityTransferFunction.Points = [255.0, 0.0, 0.5, 0.0, 255.03125, 1.0, 0.5, 0.0];\n{vtkdisplayname}.SetScalarBarVisibility(renderView1, True)"
        #print(property_str)# trace defaults for the display properties.
        exec(property_str)

        # reset view to fit data
        renderView1.ResetCamera(False)
        # get opacity transfer function/opacity map for 'scalars'
        scalarsPWF = GetOpacityTransferFunction('scalars')
        # get 2D transfer function for 'scalars'
        scalarsTF2D = GetTransferFunction2D('scalars')

        scalar_str = f"ColorBy({vtkdisplayname}, None);\nHideScalarBarIfNotNeeded(scalarsLUT, renderView1)"
        #print(scalar_str)# turn off scalar coloring# Hide the scalar bar for this color map if no visible data is colored by it.
        exec(scalar_str)    

        opacity_str = f'{vtkdisplayname}.Opacity = 0.5'
        #print(opacity_str)# Properties modified on total_filterobjDisplay
        exec(opacity_str)

        color_str = f'{vtkdisplayname}.AmbientColor = [{r}, {g}, {b}];\n{vtkdisplayname}.DiffuseColor = [{r}, {g}, {b}]'
        #print(color_str)# change solid color
        exec(color_str)
        
        ss += 1

    # get layout
    layout1 = GetLayout()

    # layout/tab size in pixels
    layout1.SetSize(1932, 1098)

    # current camera placement for renderView1
    renderView1.CameraPosition = [263.5005099773407, -1255.3709390668257, 227.53300285339355]
    renderView1.CameraFocalPoint = [263.5005099773407, 153.83953976631165, 227.53300285339355]
    renderView1.CameraViewUp = [0.0, 0.0, 1.0]
    renderView1.CameraParallelScale = 364.7305104800586

    # save screenshot
    SaveScreenshot(figfile, renderView1, ImageResolution=[1932, 1098], 
        Quality=100)
    
    # set active source
    active_str = ''
    for v in v_str:
        vtkname = v.replace('-','') 
        active_str += f'SetActiveSource({vtkname});'
    #print(active_str)
    exec(active_str)
    
    # destroy aCAvtk
    destroy_str = ''
    for v in v_str:
        vtkname = v.replace('-','') 
        destroy_str += f'Delete({vtkname});del {vtkname};'
    #print(destroy_str)
    exec(destroy_str)