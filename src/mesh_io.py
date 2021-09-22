import os
import vtk
import numpy as np
from vtk.util import numpy_support

vtk_readers = {
    "vtk": vtk.vtkGenericDataObjectReader,
    # "vtk": vtk.vtkPolyDataReader,
    "vtp": vtk.vtkXMLPolyDataReader,
    "obj": vtk.vtkMNIObjectReader
}

vtk_writers = {
    "vtk": vtk.vtkPolyDataWriter,
    "vtp": vtk.vtkXMLPolyDataWriter,
    "obj": vtk.vtkMNIObjectWriter
}

def get_extension(fname):
    path, ext = os.path.splitext(fname)
    return ext[1:].lower()

def read_mesh(fname):
    vtk_reader = vtk_readers[get_extension(fname)]()
    vtk_reader.SetFileName(fname)
    vtk_reader.Update()

    # if vtk_reader.IsFileUnstructuredGrid():
    #     data = vtk_reader.GetUnstructuredGridOutput()

    #     vertices = numpy_support.vtk_to_numpy(data.GetPoints().GetData())
    #     faces = numpy_support.vtk_to_numpy(data.GetCells().GetData())
    #     faces = np.reshape(faces, (int(faces.size / (faces[0]+1)), faces[0]+1))
    #     faces = faces[:, 1:]    

    data = vtk_reader.GetOutput()

    vertices = numpy_support.vtk_to_numpy(data.GetPoints().GetData())
    faces = numpy_support.vtk_to_numpy(data.GetPolys().GetData())
    faces = np.reshape(faces, (int(faces.size / (faces[0]+1)), faces[0]+1))
    faces = faces[:, 1:]

    return [vertices, faces]

def write_mesh(mesh, fname):

    featureVertices = []
    for v in mesh.Vertices:
        if v.feature:
            featureVertices.append(v.id)

    _mesh = vtk.vtkPolyData()
    num_points = len(mesh.Vertices)
    cell_dim = len(mesh.Faces[0].Vertices)
    num_cells = len(mesh.Faces)

    v_str = [' '.join(str(e) for e in list(mesh.Vertices[i].xyz())) + '\n' for i in range(num_points)]
    c_str = ['{dim} '.format(dim=cell_dim) + ' '.join(str(e) for e in list(mesh.Faces[i].Vertices)) + '\n' for i in range(num_cells)]
    type_str = ['5\n' for i in range(num_cells)]
    
    file_obj = open(fname, 'w')
    file_obj.write('# vtk DataFile Version 3.0\n{f}\nASCII\n\nDATASET UNSTRUCTURED_GRID\n'.format(f=fname))
    
    file_obj.write('POINTS {p} double\n'.format(p=num_points))
    file_obj.writelines(v_str)

    file_obj.write('CELLS {n_cells} {d_cells}\n'.format(n_cells=num_cells, d_cells=num_cells*(cell_dim+1)))
    file_obj.writelines(c_str)

    file_obj.write('CELL_TYPES {n_cells}\n'.format(n_cells=num_cells))
    file_obj.writelines(type_str)

    cell_dim = 1
    num_cells = len(featureVertices)

    v_str = [' '.join(str(e) for e in list(mesh.Vertices[i].xyz())) + '\n' for i in range(num_points)]
    c_str = ['{dim} '.format(dim=cell_dim) + ' '.join(str(e) for e in list([featureVertices[i]])) + '\n' for i in range(num_cells)]
    type_str = ['1\n' for i in range(num_cells)]
    
    file_obj = open('feature_vertices.vtk', 'w')
    file_obj.write('# vtk DataFile Version 3.0\n{f}\nASCII\n\nDATASET UNSTRUCTURED_GRID\n'.format(f=fname))
    
    file_obj.write('POINTS {p} double\n'.format(p=num_points))
    file_obj.writelines(v_str)

    file_obj.write('CELLS {n_cells} {d_cells}\n'.format(n_cells=num_cells, d_cells=num_cells*(cell_dim+1)))
    file_obj.writelines(c_str)

    file_obj.write('CELL_TYPES {n_cells}\n'.format(n_cells=num_cells))
    file_obj.writelines(type_str)

    cell_dim = len(mesh.Faces[0].Vertices)
    n = 0
    for patch in mesh.Patches:
        fname = 'patch_{p}.vtk'.format(p=n)
        num_cells = len(patch)
        v_str = [' '.join(str(e) for e in list(mesh.Vertices[i].xyz())) + '\n' for i in range(num_points)]
        c_str = ['{dim} '.format(dim=cell_dim) + ' '.join(str(e) for e in list(patch[i].Vertices)) + '\n' for i in range(num_cells)]
        type_str = ['5\n' for i in range(num_cells)]
        
        file_obj = open(fname, 'w')
        file_obj.write('# vtk DataFile Version 3.0\n{f}\nASCII\n\nDATASET UNSTRUCTURED_GRID\n'.format(f=fname))
        
        file_obj.write('POINTS {p} double\n'.format(p=num_points))
        file_obj.writelines(v_str)

        file_obj.write('CELLS {n_cells} {d_cells}\n'.format(n_cells=num_cells, d_cells=num_cells*(cell_dim+1)))
        file_obj.writelines(c_str)

        file_obj.write('CELL_TYPES {n_cells}\n'.format(n_cells=num_cells))
        file_obj.writelines(type_str)
        n += 1


    # mesh_vertices = np.array([list(mesh.Vertices[i].xyz()) for i in range(len(mesh.Vertices))])
    # mesh_cells = np.array([[len(mesh.Faces[0].Vertices)] + list(mesh.Faces[i].Vertices) for i in range(len(mesh.Faces))])
    # mesh_cells = numpy_support.numpy_to_vtkIdTypeArray(mesh_cells.flatten())
    
    # vertices = vtk.vtkPoints()
    # vertices.SetData(numpy_support.numpy_to_vtk(mesh_vertices))

    # cells = vtk.vtkCellArray()
    # cells.SetCells(len(mesh.Faces), mesh_cells)

    # _mesh.SetPoints(vertices)
    # _mesh.SetPolys(cells)
    
    # vtk_writer = vtk_writers[get_extension(fname)]()
    # vtk_writer.SetFileName(fname)
    # vtk_writer.SetInputData(_mesh)
    # vtk_writer.Write()


