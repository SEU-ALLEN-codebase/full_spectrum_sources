#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Fig4b
usage:
python Fig4b_cell_viz.py --folder /Users/jiangshengdian/Desktop/Daily/PhD_project/Platform/fullmorpho/py/morpho_showdata/0

'''
import pandas as pd
import numpy as np
import os
import math
import argparse
import vtk
import vtkmodules.vtkInteractionStyle
import vtkmodules.vtkRenderingOpenGL2
import glob
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkRenderingCore import (
     vtkActor,
     vtkPolyDataMapper,
     vtkRenderWindow,
     vtkRenderWindowInteractor,
     vtkRenderer,
     vtkWindowToImageFilter
)
from vtkmodules.vtkCommonCore import(
     vtkPoints,
     vtkUnsignedCharArray
)
from vtkmodules.vtkCommonDataModel import(
     vtkCellArray,
     vtkLine,
     vtkPolyData
)
from vtkmodules.vtkIOImage import (
    vtkBMPWriter,
    vtkJPEGWriter,
    vtkPNGWriter,
    vtkPNMWriter,
    vtkPostScriptWriter,
    vtkTIFFWriter
)

def WriteImage(fileName, renWin, rgba=True):
    '''
    Write the render window view to an image file.

    Image types supported are:
     BMP, JPEG, PNM, PNG, PostScript, TIFF.
    The default parameters are used for all writers, change as needed.

    :param fileName: The file name, if no extension then PNG is assumed.
    :param renWin: The render window.
    :param rgba: Used to set the buffer type.
    :return:
    '''

    import os

    if fileName:
        # Select the writer to use.
        path, ext = os.path.splitext(fileName)
        ext = ext.lower()
        if not ext:
            ext = '.png'
            fileName = fileName + ext
        if ext == '.bmp':
            writer = vtkBMPWriter()
        elif ext == '.jpg':
            writer = vtkJPEGWriter()
        elif ext == '.pnm':
            writer = vtkPNMWriter()
        elif ext == '.ps':
            if rgba:
                rgba = False
            writer = vtkPostScriptWriter()
        elif ext == '.tiff':
            writer = vtkTIFFWriter()
        else:
            writer = vtkPNGWriter()

        windowto_image_filter = vtkWindowToImageFilter()
        windowto_image_filter.SetInput(renWin)
        windowto_image_filter.SetScale(1)  # image quality
        if rgba:
            windowto_image_filter.SetInputBufferTypeToRGBA()
        else:
            windowto_image_filter.SetInputBufferTypeToRGB()
            # Read from the front buffer.
            windowto_image_filter.ReadFrontBufferOff()
            windowto_image_filter.Update()

        writer.SetFileName(fileName)
        writer.SetInputConnection(windowto_image_filter.GetOutputPort())
        writer.Write()
    else:
        raise RuntimeError('Need a filename.')
parser = argparse.ArgumentParser()
parser.add_argument('--v', help='vtk file', type=str)
parser.add_argument('--s', help='swc file', type=str)
parser.add_argument('--folder', help='swc folder', type=str)
args = parser.parse_args()

reader=vtk.vtkPolyDataReader()
vtkrootpath='/Users/jiangshengdian/Desktop/Daily/PhD_project/Platform/ccf_vtks/root.vtk'
reader.SetFileName(vtkrootpath)
# reader.SetFileName(args.v)
reader.Update()
polydata=reader.GetOutput()
colors=vtkNamedColors()

mapper=vtkPolyDataMapper()
mapper.SetInputConnection(reader.GetOutputPort())

actor=vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetColor(colors.GetColor3d('Lightgrey'))
actor.GetProperty().SetOpacity(0.2)

renderer=vtkRenderer()
renderWindow=vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor=vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderer.AddActor(actor)
renderer.SetBackground(colors.GetColor3d('White'))
renderer.GetActiveCamera().Pitch(90)
renderer.GetActiveCamera().SetViewUp(0,0,1)
renderer.ResetCamera()

colordict={0:['Red'],1:['Blue'],2:['Dim_grey'],3:['Green'],4:['Magenta'],5:['Lime'],6:['Yellow'],7:['Blue_violet'],8:['Carrot'],9:['Orchid'],
          10:['Violet'],11:['Chocolate'],12:['Dodgerblue'],13:['Cyan'],14:['Maroon'],15:['Black']}
#colordict={0:['Maroon','Red'],1:['Blue','Dodgerblue'],
           #2:['Green','Green_pale'],3:['Chocolate','Peru'],4:['Turquoise_dark','cyan'],
           #5:['Magenta','Orchid'],6:['Gold','Yellow'],7:['Violet_dark','Violet'],8:['Black','Grey']}
### swc
N=0
for fd in sorted(glob.glob(os.path.join(args.folder,'*swc'))):
    linesPolyData=vtkPolyData()
    print(fd)
    swc=np.loadtxt(fd,comments='#')
        #name=os.path.basename(f_path)
        #n=int(name[-5:-4])
    pts=vtkPoints()
    XYZ=swc[:,2:5]
    for i in range(len(XYZ)):
        pts.InsertNextPoint(XYZ[i])
    linesPolyData.SetPoints(pts)
    lines=vtkCellArray()
    idlist=list(swc[:,0])
    for i in range(len(swc)):
        if swc[i,6]!=-1:
            if swc[i,6] not in idlist:
                continue
            p=idlist.index(swc[i,6])
            line=vtkLine()
            line.GetPointIds().SetId(0,i)
            line.GetPointIds().SetId(1,p)
            lines.InsertNextCell(line)
    linesPolyData.SetLines(lines)
    colorL=vtkUnsignedCharArray()
    mapper2=vtkPolyDataMapper()
    mapper2.SetInputData(linesPolyData)
    actor2=vtkActor()
    actor2.SetMapper(mapper2)
    actor2.GetProperty().SetLineWidth(4)
    actor2.GetProperty().SetColor(colors.GetColor3d(colordict[N][0]))
    renderer.AddActor(actor2)
    N=N+1

renderWindow.SetSize(1200,1200)
renderWindow.Render()
renderWindow.SetWindowName('Brain')
outname=args.folder.split('/')[-1]
topath='../data/Fig4/morpho_showdata'
tofile=os.path.join(topath,'cluster'+outname+'.png')
WriteImage(tofile,renWin=renderWindow,rgba=False)
renderWindowInteractor.Start()

