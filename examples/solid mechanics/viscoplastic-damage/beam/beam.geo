// Gmsh project created on Thu Nov 16 22:26:03 2017
SetFactory("OpenCASCADE");
//+
SetFactory("Built-in");
//+
SetFactory("Built-in");
//+
SetFactory("Built-in");
//+
SetFactory("Built-in");
//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 1, 1, 10};
//+
Physical Surface("Clamped") = {5};
//+
Physical Surface("Zsym") = {6};
//+
Physical Volume("beam") = {1};
