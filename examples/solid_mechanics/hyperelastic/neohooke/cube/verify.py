# Script inspired by
# https://stackoverflow.com/questions/21630987/probing-sampling-interpolating
# -vtk-data-using-python-tvtk-or-mayavi

from __future__ import print_function

import sys

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

displacement = 'displacement'

def is_close(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def read_vtk_file(filename):
    # Read the vtk file with an unstructured grid
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader

def create_line(p1, p2, sample_points):
    # Create the line along which you want to sample
    line = vtk.vtkLineSource()
    line.SetResolution(sample_points)
    line.SetPoint1(p1)
    line.SetPoint2(p2)
    line.Update()
    return line

def interpolate_over_line(line, reader):

    # Interpolate the data from the VTK-file on the created line.

    # vtkProbeFilter, the probe line is the input, and the underlying dataset
    # is the source.
    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(line.GetOutputPort())
    probe.SetSourceData(reader.GetOutput())
    probe.Update()

    # Get the data from the VTK-object (probe) to an numpy array
    q = vtk_to_numpy(probe.GetOutput().GetPointData().GetArray(displacement))

    samples_on_line = probe.GetOutput().GetNumberOfPoints()

    # initialise the points on the line
    x = np.zeros(samples_on_line)
    y = np.zeros(samples_on_line)
    z = np.zeros(samples_on_line)
    points = np.zeros((samples_on_line , 3))

    # Get the coordinates of the points on the line
    for i in range(samples_on_line):
        x[i], y[i], z[i] = probe.GetOutput().GetPoint(i)
        points[i, 0] = x[i]
        points[i, 1] = y[i]
        points[i, 2] = z[i]
    return points,q

if len(sys.argv) < 3:
    print("A file name and a list of fields must be provided")
    print('Number of arguments:', len(sys.argv), 'arguments.')
    print('Argument List:', sys.argv[1::])
    quit()

file_name = sys.argv[1]
field_requests = sys.argv[2::]

start_point = [0.0, 0.0, 0.0]
end_point = [0.0, 0.0, 1.0]

interpolation_points = 10

reader = read_vtk_file(file_name)
line = create_line(start_point, end_point, interpolation_points)

if displacement in field_requests:

    points, results = interpolate_over_line(line, reader)

    if not is_close(max([max(sublist) for sublist in results]), 0.5):
        print('maximum displacement incorrect', file=sys.stderr)
        sys.exit(1)

    m, b = np.polyfit(points[:, 2], results[:, 2], 1)

    if not is_close(m, 0.5, 1.0e-7) and not is_close(b, 0.0, 1.0e-7):
        print('displacement gradient is incorrect', file=sys.stderr)
        sys.exit(1)

else:
    print('displacement not in output file', file=sys.stderr)
    sys.exit(1)

print('All checks passed')
