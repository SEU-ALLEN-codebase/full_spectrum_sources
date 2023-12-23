#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Fig4c
usage:
python main.py --p /Users/jiangshengdian/Desktop/Daily/PhD_project/Platform/fullmorpho/data/soma\ marker/soma_0_vs_1.marker

'''
import pandas as pd
import numpy as np
import os
import math
import argparse
import vtk
import vtkmodules.vtkInteractionStyle
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkFiltersSources import vtkSphereSource
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
from vtkmodules.vtkIOImage import(
     vtkPNGWriter
)

# colordict={0:'Red',1:'Blue',2:'Lime',3:'Carrot',4:'Magenta',5:'Cyan',6:'Yellow',7:'Blue_violet',8:'Dim_grey',9:'Orchid'}
colordict={0:['Red'],1:['Blue'],2:['Red'],3:['Green'],4:['Blue_violet'],5:['Lime'],6:['Yellow'],7:['Blue_violet'],8:['Carrot'],9:['Orchid'],
          10:['Violet'],11:['Chocolate'],12:['Dodgerblue'],13:['Cyan'],14:['Maroon'],15:['Black']}

parser = argparse.ArgumentParser()
parser.add_argument('--v', help='vtk file', type=str)
parser.add_argument('--p', help='point file', type=str)
args = parser.parse_args()

colors=vtkNamedColors()
#points=vtkPoints()
p_file=np.loadtxt(args.p, usecols=(0,1,2,4), skiprows=1,delimiter=',')

reader=vtk.vtkPolyDataReader()
vtkrootpath='/Users/jiangshengdian/Desktop/Daily/PhD_project/Platform/ccf_vtks/root.vtk'
reader.SetFileName(vtkrootpath)
# reader.SetFileName(args.v)
reader.Update()
polydata=reader.GetOutput()
      
mapper=vtkPolyDataMapper()
mapper.SetInputConnection(reader.GetOutputPort())
mapper.ScalarVisibilityOff()

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



#vertices=vtkCellArray()
# N=0
for i in range(p_file.shape[0]):
    #pid=points.InsertNextPoint(p_file[i])
    #vertices.InsertNextCell(1)
    #vertices.InsertCellPoint(pid)
    source = vtkSphereSource()
    source.SetRadius(4.2)
    source.SetCenter(p_file[i,0],p_file[i,1],p_file[i,2])
    source.SetPhiResolution(100)
    source.SetThetaResolution(100)
    mapper1=vtkPolyDataMapper()
    mapper1.SetInputConnection(source.GetOutputPort())
    actor1=vtkActor()
    actor1.SetMapper(mapper1)
    actor1.GetProperty().SetColor(colors.GetColor3d(colordict[int(p_file[i,3]%15)][0]))
    actor1.GetProperty().SetOpacity(0.5)
    renderer.AddActor(actor1)
#     N+=1
#     if N > 15:
#         N=0


#point=vtkPolyData()
#point.SetPoints(points)
#point.SetVerts(vertices)
print('render start....')
renderer.SetBackground(colors.GetColor3d('White'))
# renderer.GetActiveCamera().SetPosition(0,0,-1)
renderer.GetActiveCamera().Pitch(90)
# renderer.GetActiveCamera().Roll(180)
renderer.GetActiveCamera().SetViewUp(0,0,1)
renderer.ResetCamera()

renderWindow.SetSize(800,800)
renderWindow.Render()
renderWindow.SetWindowName('Soma')
renderWindowInteractor.Start()

## save
windowto_image_filter=vtkWindowToImageFilter()
windowto_image_filter.SetInput(renderWindow)
windowto_image_filter.SetScale(1)
writer=vtkPNGWriter()
topath='./soma_marker_in_vtk'
tofile=os.path.join(topath,'C'+str(int(p_file[0,3]))+'_vs_C'+str(int(p_file[p_file.shape[0]-1,3]))+'.png')
writer.SetFileName(tofile)
# writer.SetFileName('C'+str(int(p_file[0,3]))+'_vs_C'+str(int(p_file[p_file.shape[0]-1,3]))+'.png')
writer.SetInputConnection(windowto_image_filter.GetOutputPort())
writer.Write()

