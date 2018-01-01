#!/usr/bin/env python

# From https://stackoverflow.com/questions/21630987/probing-sampling-interpolating-vtk-data-using-python-tvtk-or-mayavi

# user Newfarmer

import numpy as np
from vtk.util import numpy_support as VN
from matplotlib import pyplot as plt
import vtk

def readVTK(filename):
    #read the vtk file with an unstructured grid
    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    return reader

def createLine(p1,p2,numPoints):
    # Create the line along which you want to sample
    line = vtk.vtkLineSource()
    line.SetResolution(numPoints)
    line.SetPoint1(p1)
    line.SetPoint2(p2)
    line.Update()
    return line

def probeOverLine(line,reader):
    #Interpolate the data from the VTK-file on the created line.
    data = reader.GetOutput()
    # vtkProbeFilter, the probe line is the input, and the underlying dataset is the source.
    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(line.GetOutputPort())
    probe.SetSource(data)
    probe.Update()
    #get the data from the VTK-object (probe) to an numpy array
    q=VN.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray('U'))
    numPoints = probe.GetOutput().GetNumberOfPoints() # get the number of points on the line
    #intialise the points on the line
    x = np.zeros(numPoints)
    y = np.zeros(numPoints)
    z = np.zeros(numPoints)
    points = np.zeros((numPoints , 3))
    #get the coordinates of the points on the line
    for i in range(numPoints):
        x[i],y[i],z[i] = probe.GetOutput().GetPoint(i)
        points[i,0]=x[i]
        points[i,1]=y[i]
        points[i,2]=z[i]
    return points,q

def setZeroToNaN(array):
    # In case zero-values in the data, these are set to NaN.
    array[array==0]=np.nan
    return array

#Define the filename of VTK file
filename='a-VTK-file.vtk'

#Set the points between which the line is constructed.
p1=[0.0,-0.1,0.0]
p2=[0.0,-0.1,1.0]

#Define the numer of interpolation points
numPoints=100

reader = readVTK(filename) # read the VTKfile
line=createLine(p1,p2,numPoints) # Create the line
points,U =  probeOverLine(line,reader) # interpolate the data over the line

U = setZeroToNaN(U) # Set the zero's to NaN's
plt.plot(points[:,2],U[:,0]) #plot the data
plt.show()
