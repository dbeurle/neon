
// Create a simple beam element

Point(1) = {0, 0, 0, 1.0};
Point(2) = {1.0, 0, 0, 1.0};

Line(1) = {1, 2};
//+
Physical Line("beam") = {1};
//+
Physical Point("fix") = {1};
//+
Physical Point("load") = {2};
