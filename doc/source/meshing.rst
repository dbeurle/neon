Meshing
=======

In order to construct the system of equations describing the deformation of a body, it must first be spatially discretised by a mesh.  This involves providing a CAD geometry file to a pre-processor (such as gmsh) to obtain the body in terms of points and element connectivities.

Neon uses a utility ``gmshreader`` to perform the conversion from Gmsh ``.msh`` format to a ``json`` representation with the file extension ``.mesh``.  This program can be found `here  <https://www.github.com/dbeurle/GmshReader>`_.  It is important to note here that `neon` expects the mesh indices as zero-based.  

The mesh file produced must be in the same directory as the executable when reading the input file.  The examples show the use of this, including the original Gmsh geometry and mesh data.


Element Options
---------------

The evaluation of the integrals in the finite element method use numerical quadrature rules specialised for each element.  For computational efficiency reasons, these could be under-integrated or fully-integrated.  However, in the current state reduced integration rules produce a rank-deficient element stiffness matrix and are therefore not recommended for use until stabilisation is applied.  For each mesh type, the element options can be specified ::

    "ElementOptions" : {
        "Quadrature" : "Full"
    }

with reduced integration selected when ``"Quadrature" : "Reduced"``.
