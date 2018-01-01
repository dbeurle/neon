#!/usr/bin/env python

#
#  ExtractVtkDataLine.py
#
#  The MIT License (MIT)
#
#  Copyright (c) 2016, Lluis Jofre Cruanyes. All rights reserved.
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in
#  all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
#
#
#  Created on Thursday November 17, 2016
#  Author: Lluis Jofre Cruanyes
#  PSAAP II, Stanford University
#

#
#  Information:
#
#    Extracts data along a given line (line) from a binary vtk-format file (input) into an ASCII file (output)
#
#    Usage: ./ExtractVtkDataLine.py input-file-path output-file-path QoI start-point-line (X, Y, Z values) end-point-line (X, Y, Z values) number-of-extraction-points
#
#    Important: (1) Load a python module, (2) load a vtk module, (3) chmod +x ExtractVtkDataLine.py
#


#---------------------------------------------------------------------------------------------------#

import sys
import numpy as np
import vtk
import vtk.util.numpy_support as VN


def readVtkFile( input_file ):

    # Read binary vtk file with rectilinear grid
    reader = vtk.vtkRectilinearGridReader()
    reader.SetFileName( input_file )
    reader.ReadAllScalarsOn()
    reader.Update()

    return reader


def writeAsciiFile( output_file, points, qoi_data ):

    # Open ASCII output file
    writer = open( output_file, 'w' )

    # Write data to output file
    num_rows = len( points )
    for row in range( 0, num_rows ):
        data_string = '%f %f %f %f' %( points[row][0], points[row][1], points[row][2], qoi_data[row] )
        writer.write( data_string )
	writer.write( '\n' )


def createExtractionLine( p1, p2, num_points ):

    # Create line along which to extract data
    line = vtk.vtkLineSource()
    line.SetResolution( num_points )
    line.SetPoint1( p1 )
    line.SetPoint2( p2 )
    line.Update()

    return line


def extractDataOverLine( line, reader, qoi_extract ):

    # Interpolate data from vtk file to extraction line
    data = reader.GetOutput()
    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection( line.GetOutputPort() )
    probe.SetSourceData( data )
    probe.Update()

    # Get data from the vtk object (probe) to a numpy array
    qoi_data = VN.vtk_to_numpy( probe.GetOutput().GetPointData().GetArray( qoi_extract ) )
    num_points = probe.GetOutput().GetNumberOfPoints()

    # Initialize line points
    x = np.zeros( num_points )
    y = np.zeros( num_points )
    z = np.zeros( num_points )
    points = np.zeros( ( num_points, 3 ) )

    # Get coordinates of line points
    for i in range( num_points ):
        x[i], y[i], z[i] = probe.GetOutput().GetPoint( i )
        points[i,0] = x[i]
        points[i,1] = y[i]
        points[i,2] = z[i]

    return points, qoi_data

#---------------------------------------------------------------------------------------------------#


# Input arguments crash message
if len( sys.argv ) < 10:
    print sys.argv[0], '<input_file> <output_file> <QoI> <point1_X> <point1_Y> <point1_Z> <point2_X> <point2_Y> <point2_Z> <num_points>'
    sys.exit()

# Manage input arguments
input_file  = sys.argv[1]		# Binary vtk input file
output_file = sys.argv[2]		# ASCII output file
qoi_extract = sys.argv[3]		# Quantity of interest to extract
point1_X    = float( sys.argv[4] )	# Point 1 x-coordinate
point1_Y    = float( sys.argv[5] )	# Point 1 y-coordinate
point1_Z    = float( sys.argv[6] )	# Point 1 z-coordinate
point2_X    = float( sys.argv[7] )	# Point 2 x-coordinate
point2_Y    = float( sys.argv[8] )	# Point 2 y-coordinate
point2_Z    = float( sys.argv[9] )	# Point 2 z-coordinate
num_points  = int( sys.argv[10] )	# Number of extraction points

# Create line start and end points
p1 = [ point1_X, point1_Y, point1_Z ]
p2 = [ point2_X, point2_Y, point2_Z ]

# Read binary vtk input file
reader = readVtkFile( input_file )

# Create extraction line
line = createExtractionLine( p1, p2, num_points )

# Extract data along line
points, qoi_data = extractDataOverLine( line, reader, qoi_extract )

# Write extracted data to output file
writeAsciiFile( output_file, points, qoi_data )
