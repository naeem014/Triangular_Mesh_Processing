# Triangular_Mesh_Processing
A Python Implementation for Triangular Surface Mesh Processing. Mesh Coarsening is applied on the input triangle mesh via the cost minimization of an objective function that takes into account the energy of vertices to rank edges for decimation.

## Usage:
There are 4 main files in the source folder. mesh.py contains code for mesh data structure, simplification and smoothing. mesh_io.py, pqueue.py and mesh_quality.py contain implementations for minor stuff. The program can be run by navigating to the directory containing input mesh file used for the program (e.g. h0c3_mesh.vtp) The input mesh must be a triangular mesh. 

## Info:
The mesh processing program implements a mesh data structure for certain tasks e.g. mesh simplification. The mesh.py contains the mesh data structure as described below:

### Mesh:
The mesh is composed of vertices, edges and faces. Each face has 3 vertices and 3 edges. The mesh data structure consists of a list of vertices, edges and faces that make up the entire mesh. In addition the mesh data structure contains attricutes and functions to perform certain operations such as reading an input mesh, adding a vertex, edge or face to the mesh, building mesh connectivity, building normals data of the mesh, extracting pathces and features of the mesh, simplification and smoothing of the mesh.

## Connectivity:
Each building block of the mesh i.e. vertex, edge, face is connected to other building blocks in order to define the geometry of the mesh. Following is the description of the data structure for these building blocks:

### Vertex:
A vertex consists of three geometric coordinates (x,y,z) which define its position in space. In addition, a vertex also has a list of neighboring vertices, edges and faces.

### Edge:
An edge consists of two vertices that are connected to each other to build mesh faces. In addition, an edge also has a list of neighboring vertices, edges and faces.

### Face:
An face consists of three vertices and three edges to define its shape. In addition, a face also has a list of neighboring vertices, edges and faces.

# Results
![alt text](https://github.com/naeem014/Triangular_Mesh_Processing/blob/main/media/tri_result.png)

