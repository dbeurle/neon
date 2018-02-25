
#pragma once

#include <string>

// A collection of functions that can be used to test a cube mesh
// including a mesh and input data

inline std::string cube_mesh_json()
{
    return "{\"Elements\":[{\"Indices\":[9,10,11,12,13,14,15,16,17],\"Name\":\"bottom\","
           "\"NodalConnectivity\":[[1,18,36,27],[27,36,37,26],[26,37,17,0],[18,19,38,36],"
           "[36,38,39,"
           "37],[37,39,16,17],[19,3,14,38],[38,14,15,39],[39,15,2,16]],\"Type\":3},{"
           "\"Indices\":["
           "67,80,79,78,77,76,75,74,73,72,71,70,69,68,54,66,65,64,63,62,61,60,59,58,57,"
           "56,55],"
           "\"Name\":\"cube\",\"NodalConnectivity\":[[62,60,56,58,63,61,57,59],[29,47,63,"
           "53,6,10,"
           "41,9],[28,45,62,52,29,47,63,53],[7,23,35,24,28,45,62,52],[47,46,61,63,10,11,"
           "40,41],[45,"
           "44,60,62,47,46,61,63],[23,22,34,35,45,44,60,62],[46,31,51,61,11,4,12,40],[44,"
           "30,50,60,"
           "46,31,51,61],[22,5,21,34,44,30,50,60],[53,63,59,55,9,41,43,8],[52,62,58,54,"
           "53,63,59,55]"
           ",[24,35,33,25,52,62,58,54],[63,61,57,59,41,40,42,43],[32,20,1,27,56,48,18,36]"
           ",[35,34,"
           "32,33,62,60,56,58],[61,51,49,57,40,12,13,42],[60,50,48,56,61,51,49,57],[34,"
           "21,20,32,60,"
           "50,48,56],[55,59,39,16,8,43,15,2],[54,58,37,17,55,59,39,16],[25,33,26,0,54,"
           "58,37,17],["
           "59,57,38,39,43,42,14,15],[58,56,36,37,59,57,38,39],[33,32,27,26,58,56,36,37],"
           "[57,49,19,"
           "38,42,13,3,14],[56,48,18,36,57,49,19,38]],\"Type\":5},{\"Indices\":[45,37,38,"
           "39,40,41,"
           "42,43,44,36,46,47,48,49,50,51,52,53,18,1,2,3,4,5,6,7,8,0,19,20,21,22,23,24,"
           "25,26],"
           "\"Name\":\"sides\",\"NodalConnectivity\":[[7,24,52,28],[18,48,49,19],[19,49,"
           "13,3],[20,"
           "21,50,48],[48,50,51,49],[49,51,12,13],[21,5,30,50],[50,30,31,51],[51,31,4,12]"
           ",[1,20,48,"
           "18],[28,52,53,29],[29,53,9,6],[24,25,54,52],[52,54,55,53],[53,55,8,9],[25,0,"
           "17,54],[54,"
           "17,16,55],[55,16,2,8],[4,12,40,11],[27,32,33,26],[26,33,25,0],[20,21,34,32],["
           "32,34,35,"
           "33],[33,35,24,25],[21,5,22,34],[34,22,23,35],[35,23,7,24],[1,20,32,27],[11,"
           "40,41,10],["
           "10,41,9,6],[12,13,42,40],[40,42,43,41],[41,43,8,9],[13,3,14,42],[42,14,15,43]"
           ",[43,15,2,"
           "8]],\"Type\":3},{\"Indices\":[27,28,29,30,31,32,33,34,35],\"Name\":\"top\","
           "\"NodalConnectivity\":[[5,30,44,22],[22,44,45,23],[23,45,28,7],[30,31,46,44],"
           "[44,46,47,"
           "45],[45,47,29,28],[31,4,11,46],[46,11,10,47],[47,10,6,29]],\"Type\":3}],"
           "\"Nodes\":[{"
           "\"Coordinates\":[[0,0,0],[1,0,0],[0,1,0],[1,1,0],[1,1,1],[1,0,1],[0,1,1],[0,"
           "0,1],[0,1,"
           "0.33333333333250098],[0,1,0.66666666666578744],[0.33333333333250098,1,1],[0."
           "66666666666578744,1,1],[1,1,0.66666666666759178],[1,1,0.33333333333472071],["
           "0."
           "66666666666759178,1,0],[0.33333333333472071,1,0],[0,0.66666666666759178,0],["
           "0,0."
           "33333333333472071,0],[1,0.33333333333250098,0],[1,0.66666666666578744,0],[1,"
           "0,0."
           "33333333333250098],[1,0,0.66666666666578744],[0.66666666666759178,0,1],[0."
           "33333333333472071,0,1],[0,0,0.66666666666759178],[0,0,0.33333333333472071],["
           "0."
           "33333333333250098,0,0],[0.66666666666578744,0,0],[0,0.33333333333250098,1],["
           "0,0."
           "66666666666578744,1],[1,0.33333333333250098,1],[1,0.66666666666578744,1],[0."
           "66666666666638896,0,0.33333333333324089],[0.333333333333241,0,0."
           "33333333333398069],[0."
           "66666666666699048,0,0.66666666666638896],[0.33333333333398091,0,0."
           "66666666666699026],["
           "0.66666666666638896,0.33333333333324089,0],[0.333333333333241,0."
           "33333333333398069,0],["
           "0.66666666666699048,0.66666666666638896,0],[0.33333333333398091,0."
           "66666666666699026,0],"
           "[0.66666666666638896,1,0.66666666666699026],[0.333333333333241,1,0."
           "66666666666638896],["
           "0.66666666666699048,1,0.33333333333398069],[0.33333333333398091,1,0."
           "333333333333241],["
           "0.66666666666699026,0.33333333333250098,1],[0.33333333333398091,0."
           "33333333333250098,1],"
           "[0.66666666666638874,0.66666666666578744,1],[0.33333333333324089,0."
           "66666666666578744,1]"
           ",[1,0.33333333333250098,0.33333333333324089],[1,0.66666666666578744,0."
           "33333333333398091],[1,0.33333333333250098,0.66666666666638874],[1,0."
           "66666666666578744,"
           "0.66666666666699048],[0,0.333333333333241,0.66666666666699026],[0,0."
           "66666666666638896,"
           "0.66666666666638896],[0,0.33333333333398091,0.33333333333398091],[0,0."
           "66666666666699048,0.333333333333241],[0.66666666666658936,0."
           "33333333333299442,0."
           "33333333333348758],[0.66666666666679042,0.66666666666618857,0."
           "33333333333373422],[0."
           "33333333333348752,0.33333333333348752,0.33333333333373438],[0."
           "33333333333373433,0."
           "66666666666658925,0.33333333333348758],[0.66666666666678998,0."
           "33333333333274773,0."
           "66666666666658947],[0.66666666666658991,0.66666666666598795,0."
           "66666666666679009],[0."
           "33333333333373422,0.33333333333299431,0.66666666666678975],[0."
           "33333333333348769,0."
           "66666666666618846,0.66666666666658936]],\"Indices\":[0,1,2,3,4,5,6,7,8,9,10,"
           "11,12,13,"
           "14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,"
           "39,40,41,42,"
           "43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]}]}";
}

inline std::string material_data_json()
{
    return "{\"Name\" : \"steel\",\"ElasticModulus\" : 200.0e6, \"PoissonsRatio\" "
           ": 0.45, \"Density\" : 7800.0 }";
}

inline std::string simulation_data_json()
{
    return "{ \"BoundaryConditions\" : [ "
           //
           "{\"Name\" : \"bottom\", "
           "\"Type\" : \"Displacement\","
           "\"Time\" : [0.0, 1.0],"
           "\"x\" : [0.0, 0.0], \"y\" : [0.0, 0.0], \"z\" : [0.0, 0.0]}, "
           //
           "{\"Name\" : \"top\", "
           "\"Type\" : \"Displacement\","
           "\"Time\" : [0.0, 1.0],"
           "\"z\" : [0.0, 1.0e-3]}],"
           //
           "\"ConstitutiveModel\" : {\"Name\":\"NeoHooke\"}, "
           "\"ElementOptions\" : {\"Quadrature\" : \"Full\"}, "
           "\"Name\" : \"cube\", "
           //
           "\"Visualisation\" : {\"Fields\" : [\"Displacement\", \"CauchyStress\"]},"
           //
           "\"NonlinearOptions\" : { "
           "\"DisplacementTolerance\" : 1.0e-5, "
           "\"ResidualTolerance\" : 0.001},"
           //
           "\"Time\" : {\"Period\" : 1.0, \"Increments\": { "
           "\"Initial\" : 1.0, \"Minimum\" : 0.001, \"Maximum\" : 10.0, \"Adaptive\" : true }},"
           "\"LinearSolver\" : {\"Solver\" : \"Iterative\"}}";
}

std::string const nonlinear_options_disp_broken_json("{\"NonlinearOptions\" : { "
                                                     "\"DisacementIncrementTolerance\" : 1.0e-5, "
                                                     "\"ResidualTolerance\" : 0.001}}");

std::string const nonlinear_options_force_broken_json("{\"NonlinearOptions\" : { "
                                                      "\"DisplacementTolerance\" : "
                                                      "1.0e-5, "
                                                      "\"ResidForceTolerance\" : 0.001}}");

inline std::string simulation_data_traction_json()
{
    return "{ \"BoundaryConditions\" : [ "
           //
           "{\"Name\" : \"Ysym\", "
           "\"Type\" : \"Traction\", "
           "\"Time\" : [0.0, 1.0],"
           "\"y\" : [0.0, 1.0e-3]} ], "
           //
           "\"ConstitutiveModel\" : {\"Name\":\"NeoHooke\"}, "
           "\"ElementOptions\" : {\"Quadrature\" : \"Full\"}, "
           "\"Name\" : \"cube\"}";
}

inline std::string solver_data_json()
{
    return "{\"Solver\" : \"Iterative\", \"MaxIterations\" : 1000, "
           " \"Tolerance\" : 1e-6 }";
}
