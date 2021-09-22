import numpy as np
from meshio.vtk import _check_vtk_available, _get_extension_key, _get_vtk_readers 
from vtk import vtkMeshQuality, vtkGenericDataObjectReader

def getMeshQuality(mesh, q_param="Jacobian"):
    _check_vtk_available()
    
    # key = _get_extension_key(mesh.fname)
    # reader_class = _get_vtk_readers()[key]
    # reader = reader_class()
    reader = vtkGenericDataObjectReader()
    reader.SetFileName(mesh.fname)
    reader.Update()

    qualityFilter = vtkMeshQuality()
    qualityFilter.SetInputConnection(reader.GetOutputPort())
    qualityFilter.SetTriangleQualityMeasureToScaledJacobian()
    qualityFilter.Update()

    qualityArray = qualityFilter.GetOutput().GetCellData().GetArray("Quality")

    minScaledJacobian = 1
    maxScaledJacobian = -1
    avgScaledJacobian = 0
    for i in range(qualityArray.GetNumberOfTuples()):
        value = qualityArray.GetValue(i)
        if value < minScaledJacobian:
            minScaledJacobian = value
        if value > maxScaledJacobian:
            maxScaledJacobian = value
        avgScaledJacobian += value

    avgScaledJacobian /= qualityArray.GetNumberOfTuples()
    print("Mesh quality: ")  
    print("min. Scaled Jacobian: {min}".format(min=minScaledJacobian))  
    print("avg. Scaled Jacobian: {avg}".format(avg=avgScaledJacobian))  
    print("max. Scaled Jacobian: {max}".format(max=maxScaledJacobian))  