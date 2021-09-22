import numpy
import sys
import glm
import math
import copy
import h5py
import gc
from pqueue import Pqueue

from mesh_io import read_mesh, write_mesh
import mesh_quality

### Building the data structure for the mesh to be processed

# Class Vertex defined for object types vertices which carry the location information for each point in the mesh
class Vertex:
    def __init__(self, coords, Id):
        self.id = Id
        self.coords = glm.dvec3(coords[0], coords[1], coords[2])
        self.N_vids = set()
        self.N_eids = set()
        self.N_fids = set()
        self.isVisited = False
        self.normal = glm.dvec3(0, 0, 0)
        self.feature = False
        self.deleted = False

    # Returns the xyz coordinates of a vertex
    def xyz(self):
        return self.coords

    # Sets the xyz coordinates of a vertex
    def setCoords(self, coords):
        self.coords = coords

    # Adds a vertex reference (Vertex id) to the list of neighboring vertices
    def addNeighborVertex(self, n_vid):
        self.N_vids.add(n_vid)

    # Adds a Edge reference (Edge id) to the list of neighboring edges
    def addNeighborEdge(self, n_eid):
        self.N_eids.add(n_eid)
    
    # Adds a Face reference (Face id) to the list of neighboring faces
    def addNeighborFace(self, n_fid):
        self.N_fids.add(n_fid)


# Class Edge defined for object types edges in a mesh. Each edge connects two vertices in a mesh
class Edge:
    def __init__(self, Vertices):
        self.id = -1
        self.Vertices = Vertices
        self.N_vids = set()
        self.N_eids = set()
        self.N_fids = set()
        self.feature = False
        self.deleted = False
    
    # Adds a vertex reference (Vertex id) to the list of neighboring vertices
    def addNeighborVertex(self, n_vid):
        self.N_vids.add(n_vid)

    # Adds a Edge reference (Edge id) to the list of neighboring edges
    def addNeighborEdge(self, n_eid):
        self.N_eids.add(n_eid)
    
    # Adds a Face reference (Face id) to the list of neighboring faces
    def addNeighborFace(self, n_fid):
        self.N_fids.add(n_fid)
    

# Class Face defined for object types faces in a mesh. Each face has certain number of vertices (eg. 3 vertices in a triangle face)
# and a certain number of edges (eg. 3 edges in a triangle face)
class Face:
    def __init__(self, Vertices, Id):
        self.id = Id
        self.Vertices = Vertices
        self.Edges = set()
        self.N_vids = set()
        self.N_eids = set()
        self.N_fids = set()
        self.isVisited = False
        self.normal = glm.dvec3(0, 0, 0)
        self.patchID = 0
        self.feature = False
        self.deleted = False
    
    # Adds an edge that belongs to the face
    def addEdge(self, e_id):
        self.Edges.add(e_id)
    
    # Adds a vertex reference (Vertex id) to the list of neighboring vertices
    def addNeighborVertex(self, n_vid):
        self.N_vids.add(n_vid)

    # Adds a Edge reference (Edge id) to the list of neighboring edges
    def addNeighborEdge(self, n_eid):
        self.N_eids.add(n_eid)
    
    # Adds a Face reference (Face id) to the list of neighboring faces
    def addNeighborFace(self, n_fid):
        self.N_fids.add(n_fid)


# Class Mesh defined to build a data structure for mesh processing. This class contains information about Vertices, Edges, Faces
# and their relevant neighbors for ease of processing the mesh.
class Mesh:
    def __init__(self):
        self.Vertices = []
        self.Edges = []
        self.Faces = []
        
        self.Patches = []
        self.fname = ""
        self.Area = 0.0


    # Read the mesh data from a vtk mesh file and build the mesh data structure
    def readInputFile(self, fname):
        self.fname = fname
        mesh_data = read_mesh(fname)
        _vertices = mesh_data[0]
        _faces = mesh_data[1]

        self.Vertices = [0]*len(_vertices)
        for i,v in enumerate(_vertices):
            self.addVertex(Vertex(v,i))
        
        self.Faces = [0]*len(_faces)
        for i, f in enumerate(_faces):
            self.addFace(Face(list(f),i))    

    # Write the current mesh information in an output vtk mesh file 
    def writeOutputFile(self, fname):
        pass   
    
    # Add a vertex to the mesh data structure
    def addVertex(self, vertex):
        if type(vertex) is not Vertex:
            print("Error: Please provide an input of type Vertex")
            return
        else:
            self.Vertices[vertex.id] = vertex
            
    # Add an edge to the mesh data structure
    def addEdge(self, edge):
        if type(edge) is not Edge:
            print("Error: Please provide an input of type Edge")
            return
        else:
            edge.id = len(self.Edges)
            self.Edges.append(edge)
            
    # Add a face to the mesh data structure
    def addFace(self, face):
        if type(face) is not Face:
            print("Error: Please provide an input of type Face")
            return
        else:
            self.Faces[face.id] = face
    
    # Extracts the edges in the mesh and builds the connectivity information of the mesh (e.g vertex neighbors, face neighbors etc.)
    def buildMeshConnectivity(self):
        for f in self.Faces:
            for vid in f.Vertices:
                self.Vertices[vid].addNeighborFace(f.id) # vertex-face neighboring
                self.Vertices[vid].N_vids = self.Vertices[vid].N_vids.union(set(f.Vertices).difference(set({vid}))) # vertex-vertex neighboring

        for v in self.Vertices:
            v.isVisited = False

        for v in self.Vertices:
            for vid in v.N_vids:
                v2 = self.Vertices[vid]
                if v2.isVisited:
                    continue
                new_edge = Edge(list([v.id, v2.id]))
                self.addEdge(new_edge)
                v.addNeighborEdge(new_edge.id) # vertex-edge neighboring
                v2.addNeighborEdge(new_edge.id) # vertex-edge neighboring

                new_edge.N_vids = (v.N_vids.union(v2.N_vids)).difference(set({v.id, v2.id})) # edge-vertex neighboring
                new_edge.N_fids = v.N_fids.intersection(v2.N_fids) # edge-face neighboring
            v.isVisited = True
        
        for e in self.Edges:
            e.N_eids = (self.Vertices[e.Vertices[0]].N_eids.union(self.Vertices[e.Vertices[1]].N_eids)).difference(set({e.id})) # edge-edge neighboring
            for fid in e.N_fids:
                face = self.Faces[fid]
                if len(set(face.Vertices).intersection(e.Vertices)) == 2:
                    face.addEdge(e.id) # assigning an edge to a face
                else:
                    face.addNeighborEdge(e.id) # face-edge neighboring

        for f in self.Faces:
            for vid in f.Vertices:
                f.N_vids = f.N_vids.union(self.Vertices[vid].N_vids.difference(set(f.Vertices))) # face-vertex neighboring
                f.N_fids = f.N_fids.union(self.Vertices[vid].N_fids.difference(set({f.id}))) # face-face neighboring

    
    # Set Face normal of a triangular face in mesh. The Face data structure has an attribute named normal that stores it
    def setFaceNormal(self, faceId):
        face = self.Faces[faceId]
        a = glm.dvec3(self.Vertices[face.Vertices[1]].xyz() - self.Vertices[face.Vertices[0]].xyz())
        b = glm.dvec3(self.Vertices[face.Vertices[2]].xyz() - self.Vertices[face.Vertices[0]].xyz())
        face.normal = glm.normalize(glm.cross(a,b))
    
    def getFaceNormal(self, faceId):
        return self.Faces[faceId].normal

    # Set Vertex normal of a triangular face in mesh. The Vertex data structure has an attribute named normal that stores it
    def setVertexNormal(self, vertexId):
        v = self.Vertices[vertexId]
        for fid in v.N_fids:
            v.normal += self.getFaceNormal(fid)
        if len(v.N_fids) == 0:
            return
        v.normal /= len(v.N_fids)

    def getVertexNormal(self, vertexId):
        return self.Vertices[vertexId].normal

    def buildNormalData(self):
        for f in self.Faces:
            self.setFaceNormal(f.id)
        
        for v in self.Vertices:
            self.setVertexNormal(v.id)

    # Check if planes associated with two faces fall within a user defined threshold
    def getPlaneVerdict(self, f1, f2, angle_threshold):
        return (math.acos(max(0.0, min(glm.dot(f1.normal, f2.normal), 1.0))) * 180 / math.pi) < angle_threshold
    
    
    # Check if the orientation of a vertex differs from its neighbors by a user defined threshold
    def getOrientationVerdict(self, n1, n2, angle_threshold):
        verdict = glm.dot(n1, n2)
        return verdict >= 0 and verdict < math.cos(math.radians(angle_threshold))

    # Seperate triangular faces that have similar planes in the mesh
    def extractPatches(self, angle_threshold):
        patches = []
        for f in self.Faces:
            f.isVisited = False
        for f in self.Faces:
            if f.isVisited:
                continue
            patch = []
            f_stack = []
            f_stack.append(f.id)
            while len(f_stack) != 0:
                face = self.Faces[f_stack.pop()]             
                if face.isVisited:
                    continue
                if self.getPlaneVerdict(f, face, angle_threshold):
                    patch.append(face)
                    face.patchID = len(patches)
                    face.isVisited = True
                    f_stack.extend(face.N_fids)
            patches.append(patch)
        self.Patches = patches
    
    # Identify vertices that are important 
    def getFeatureVertices(self, angle_threshold):
        for v in self.Vertices:
            v.isVisited = False
        for v in self.Vertices:
            for fid in v.N_fids:
                f = self.Faces[fid]
                if self.getOrientationVerdict(v.normal, f.normal, angle_threshold):
                    v.feature = True
                    break
        
        for e in self.Edges:
            if len(e.N_fids) == 1:
                self.Vertices[e.Vertices[0]].feature = True
                self.Vertices[e.Vertices[1]].feature = True
                e.feature = True
            
    
    # Get left and right neighbors of a vertex
    def getVertexNeighborNodes(self, vid, eid):
        e = self.Edges[eid]
        v = self.Vertices[vid]
        nvid = e.Vertices[1] if vid == e.Vertices[0] else e.Vertices[0]

        f1 = self.Faces[list(e.N_fids)[0]]
        f2 = self.Faces[list(e.N_fids)[1]]
        
        t_nv1 = list(set(f1.Vertices).difference(set(e.Vertices)))[0]
        t_nv2 = list(set(f2.Vertices).difference(set(e.Vertices)))[0]

        index1 = f1.Vertices.index(nvid)
        index2 = f1.Vertices.index(t_nv1)

        if (index1 + 1) % len(f1.Vertices) == index2:
            nvid_next = t_nv1
            nvid_prev = t_nv2
        else:
            nvid_next = t_nv2
            nvid_prev = t_nv1
        
        return nvid, nvid_prev, nvid_next
    
    # Get angles associated with left and right neighbors of a vertex
    def getVertexNeighborAngles(self, vid, nvid, nvid_prev, nvid_next):
        v_n = glm.normalize(self.Vertices[vid].xyz() - self.Vertices[nvid].xyz())
        v_n_prev = glm.normalize(self.Vertices[nvid_prev].xyz() - self.Vertices[nvid].xyz())
        v_n_next = glm.normalize(self.Vertices[nvid_next].xyz() - self.Vertices[nvid].xyz())

        alpha1 = math.acos(max(0.0, min(glm.dot(v_n, v_n_next), 1.0))) * 180 / math.pi
        alpha2 = math.acos(max(0.0, min(glm.dot(v_n, v_n_prev), 1.0))) * 180 / math.pi

        return alpha1, alpha2
    
    # Get energy of a vertex 
    def getVertexEnergy(self, vid):
        v = self.Vertices[vid]
        n = len(v.N_eids)
        if n == 0:
            return 0.0
        v_ideal = 180 - (360 / n)
        E = 0.0
        for eid in v.N_eids:
            nvid, nvid_prev, nvid_next = self.getVertexNeighborNodes(vid, eid)
            alpha1, alpha2 = self.getVertexNeighborAngles(vid, nvid, nvid_prev, nvid_next)
            E += v_ideal / (alpha1 + alpha2)
        return E

    # Get cost of an edge
    def getEdgeCost(self, eid):
        e = self.Edges[eid]
        t1 = self.getVertexEnergy(e.Vertices[0]) + self.getVertexEnergy(e.Vertices[1])
        t2 = (self.getFaceArea(list(e.N_fids)[0]) + self.getFaceArea(list(e.N_fids)[1])) / self.getArea()
        cost = t1 + t2
        return cost

    # Get mesh energy
    def getEnergy(self):
        E = 0.0
        for v in self.Vertices:
            if v.feature or v.deleted:
                continue
            E += self.getVertexEnergy(v.id) - len(v.N_eids)
        return E

    # Simple Laplacian smoothing
    def smoothLaplacian(self):
        for it in range(0, 25):
            new_coords = []
            for v in self.Vertices:
                if v.feature:
                    new_coords.append(v.xyz())
                    continue
                new_coord = glm.dvec3(0, 0, 0)
                for nvid in v.N_vids:
                    new_coord += self.Vertices[nvid].xyz()
                new_coords.append(new_coord / len(v.N_vids))
            for v in self.Vertices:
                v.setCoords(new_coords[v.id])

    # Angle-based smoothing
    def smoothAngleBased(self):
        while True:
            previousE = self.getEnergy()
        
            new_coords = numpy.empty(len(self.Vertices), dtype=glm.dvec3)
            for v in self.Vertices:
                new_coords[v.id] = glm.dvec3(0.0, 0.0, 0.0)
                if len(v.N_fids) == 0:
                    continue
                feature = v.feature or v.deleted
                if feature == True:
                    continue
                else:
                    new_coord = []
                    weights = []
                    for eid in v.N_eids:   
                        nvid, nvid_prev, nvid_next = self.getVertexNeighborNodes(v.id, eid)
                        alpha1, alpha2 = self.getVertexNeighborAngles(v.id, nvid, nvid_prev, nvid_next)

                        beta = (alpha2 - alpha1) * 0.5
                        r = self.Vertices[nvid_next].xyz() - v.xyz()
                        l = math.fabs(beta / alpha1) * glm.length(r)
                        if beta > 0:
                            r = self.Vertices[nvid_prev].xyz() - v.xyz()
                            l = math.fabs(beta / alpha2) * glm.length(r)
                        new_coord.append(l * glm.normalize(r))
                        weights.append(math.fabs(beta))
                    
                    weight_agg = sum(weights)
                    vertex_weight = 0.0
                    for i in range(len(new_coord)):
                        w = (1 - weights[i] / weight_agg)
                        new_coords[v.id] += w * new_coord[i]
                        vertex_weight += w
                    new_coords[v.id] /= vertex_weight

            for v in self.Vertices:
                if len(v.N_fids) == 0:
                    continue
                feature = v.feature or v.deleted
                if feature == True:
                    continue
                new_coords[v.id] = self.remapPointToSurface(v.xyz() + new_coords[v.id], v.id)
                v.setCoords(new_coords[v.id])
            
            currentE = self.getEnergy()
            print("current energy: {ce}, previous energy: {pe}".format(ce=currentE,pe=previousE))
            if previousE - currentE < 1e-2:
                break

    def smooth(self, algorithm):
        if type(algorithm) is not str:
            print("Please provide a correct smoothing algorithm. The choices are:\nLaplacian\nAngleBased")
            return
        if algorithm == "Laplacian":
            self.smoothLaplacian()
        elif algorithm == "AngleBased":
            self.smoothAngleBased()

    def getTriangleArea(self, vertices):
        return 0.5 * glm.length(glm.cross(vertices[1] - vertices[0], vertices[2] - vertices[0]))

    def getFaceArea(self, fid):
        f = self.Faces[fid]
        return self.getTriangleArea([self.Vertices[f.Vertices[0]].xyz(), self.Vertices[f.Vertices[1]].xyz(), self.Vertices[f.Vertices[2]].xyz()])

    def setArea(self):
        for f in self.Faces:
            self.Area += self.getFaceArea(f.id)

    def getArea(self):
        return self.Area

    # Function to determine position of a point on plane most close to input mesh
    def getPointPlaneIntersection(self, p, plane, plane_normal):
        v0 = self.Vertices[plane[0]].xyz()
        v1 = self.Vertices[plane[1]].xyz()
        v2 = self.Vertices[plane[2]].xyz()

        t = (plane_normal.x * v0.x - plane_normal.x * p.x) + (plane_normal.y * v0.y - plane_normal.y * p.y) + (plane_normal.z * v0.z - plane_normal.z * p.z) 
        t /= (plane_normal.x * plane_normal.x) + (plane_normal.y * plane_normal.y) + (plane_normal.z * plane_normal.z)

        intersection = p + (t * plane_normal)

        S = self.getTriangleArea([v0, v1, v2])
        s1 = self.getTriangleArea([v1, intersection, v2])
        s2 = self.getTriangleArea([intersection, v0, v2])
        s3 = self.getTriangleArea([intersection, v1, v2])

        a = s1/S
        b = s2/S
        c = s3/S

        decision = (0 <= a and a <= 1) and (0 <= b and b <= 1) and (0 <= c and c <= 1)
        return intersection, decision

    # Remap a point to the surface of the mesh
    def remapPointToSurface(self, p, vid):
        v = self.Vertices[vid]
        remappedPoint = p
        for fid in v.N_fids:
            newPoint, decision = self.getPointPlaneIntersection(p, list(self.Faces[fid].Vertices), glm.normalize(self.Faces[fid].normal))
            if decision == True:
                remappedPoint = newPoint
                break
        
        return remappedPoint
        
    # Edge collapse operation used in simplification of mesh. Each collapse operation removes one edge, one vertex and two faces
    # The neighboring faces, edges and vertices need to be updated after a collapse operation is performed in order to avoid building
    # whole mesh again

    def collapseEdge(self, eid):
        e = self.Edges[eid]
        
        source = self.Vertices[e.Vertices[0]]
        target = self.Vertices[e.Vertices[1]]

        target.setCoords(0.5 * (source.xyz() + target.xyz()))

        canceledVertices = set({source.id})
        canceledEdges = set({e.id})
        canceledFaces = e.N_fids
        
        for fid in e.N_fids:
            f = self.Faces[fid]
            f.deleted = True
            f_edges = list(f.Edges)
            index = f_edges.index(e.id)
            e1 = self.Edges[f_edges[(index + 1) % len(f_edges)]]
            e2 = self.Edges[f_edges[(index + 2) % len(f_edges)]]

            f1_id = list(e1.N_fids)[0]
            f2_id = list(e2.N_fids)[0]
            if f1_id == f.id:
                f1_id = list(e1.N_fids)[1]
            if f2_id == f.id:
                f2_id = list(e2.N_fids)[1]
            
            f1 = self.Faces[f1_id]
            f2 = self.Faces[f2_id]

            f2.Edges = f2.Edges.difference(set({e2.id}))
            f2.Edges.add(e1.id)
            e1.N_fids = set({f1.id, f2.id})
            canceledEdges.add(e2.id)
        
        
        source.deleted = True
        for eid in canceledEdges:
            self.Edges[eid].deleted = True

        source.N_vids = source.N_vids.difference(set({target.id}))
        source.N_eids = source.N_eids.difference(canceledEdges)
        source.N_fids = source.N_fids.difference(canceledFaces)

        
        target.N_vids = target.N_vids.difference(set({source.id}))
        target.N_eids = target.N_eids.difference(canceledEdges)
        target.N_fids = target.N_fids.difference(canceledFaces)
        
        for eid in target.N_eids:
            e = self.Edges[eid]
            e.N_vids = e.N_vids.difference(set({source.id}))
            e.N_vids = e.N_vids.union(source.N_vids)
            e.N_eids = e.N_eids.difference(canceledEdges)
            e.N_eids = e.N_eids.union(source.N_eids)

        for fid in target.N_fids:
            f = self.Faces[fid]
            f.N_vids = f.N_vids.difference(set({source.id}))
            f.N_vids.union(source.N_vids.difference(set(f.Vertices)))
            f.N_eids = f.N_eids.difference(canceledEdges)
            f.N_eids = f.N_eids.union(source.N_eids)
            f.N_fids = f.N_fids.difference(canceledFaces)
            f.N_fids = f.N_fids.union(source.N_fids)
        
        for vid in source.N_vids:
            v = self.Vertices[vid]
            v.N_vids = v.N_vids.difference(set({source.id}))
            v.N_vids.add(target.id)
            v.N_eids = v.N_eids.difference(canceledEdges)
            v.N_fids = v.N_fids.difference(canceledFaces)
        
        for eid in source.N_eids:
            e = self.Edges[eid]
            e.Vertices[e.Vertices.index(source.id)] = target.id
            e.N_vids = e.N_vids.difference(set({target.id}))
            e.N_vids = e.N_vids.union(target.N_vids)
            e.N_eids = e.N_eids.difference(canceledEdges)
            e.N_eids = e.N_eids.union(target.N_eids)
        
        for fid in source.N_fids:
            f = self.Faces[fid]
            f.Vertices[f.Vertices.index(source.id)] = target.id
            f.N_vids = f.N_vids.difference(set({target.id}))
            f.N_vids = f.N_vids.union(target.N_vids.difference(set(f.Vertices)))
            f.N_eids = f.N_eids.difference(canceledEdges)
            f.N_eids = f.N_eids.union(target.N_eids.difference(set(f.Edges)))
            f.N_fids = f.N_fids.difference(canceledFaces)
            f.N_fids = f.N_fids.union(target.N_fids)

        target.N_vids = target.N_vids.union(source.N_vids)
        target.N_eids = target.N_eids.union(source.N_eids)
        target.N_fids = target.N_fids.union(source.N_fids)

        return target.id

    # Simplify the mesh using edge collapse operations where edges are sorted according to edge costs
    def simplify(self):
        p_queue = Pqueue(len(self.Edges))
        for e in self.Edges:
            if self.Vertices[e.Vertices[0]].feature or self.Vertices[e.Vertices[1]].feature:
                continue
            index = p_queue.insert([self.getEdgeCost(e.id), e.id])


        decimation_factor = 0.9
        num_elements_to_decimate = int(p_queue.size() * decimation_factor)
        print(num_elements_to_decimate)
        it = 0
        for i in range(0,num_elements_to_decimate):
            eid = p_queue.pop()[1]
            if self.Edges[eid].deleted:
                continue
            target = self.collapseEdge(eid)
            it += 1
            print("Collapsed Edge: {id} in iteration: {it_}".format(id=eid,it_=it))
            v = self.Vertices[target]
            for eid in v.N_eids:
                e = self.Edges[eid]
                if self.Vertices[e.Vertices[0]].feature or self.Vertices[e.Vertices[1]].feature:
                    continue
                c = self.getEdgeCost(eid)
                p_queue.updatePriority(p_queue.l_t[eid], c)
    
    # Remove deleted faces after mesh simplification
    def clean(self):
        print("Number of faces in input mesh: {n}".format(n=len(self.Faces)))
        deletedFaces = []
        for f in self.Faces:
            if f.deleted:
                deletedFaces.append(f.id)
        
        self.Faces = numpy.delete(self.Faces, deletedFaces)
        print("Number of faces in output mesh: {n}".format(n=len(self.Faces)))
    
    # Helper functions for debugging purposes
    def printVertexCoords(self, vertexId):
        vertex = self.Vertices[vertexId]
        print("id: {id} x: {x}, y: {y}, z: {z}".format(id=vertex.id, x=vertex.coords.x, y=vertex.coords.y, z=vertex.coords.z))

    def printEdgeData(self, edgeId):
        edge = self.Edges[edgeId]
        print("Edge id: {id}".format(id=edge.id))
        for vid in edge.Vertices:
            vertex = self.Vertices[vid]
            self.printVertexCoords(vertex)
        print("-------------------------------------")
        
    def printFaceEdges(self, faceId):
        face = self.Faces[faceId]
        res = "Edge ids in Face: {faceId} are".format(faceId=face.id)
        for eid in face.Edges:
            res += " {edgeId}".format(edgeId=eid)
        print(res)
        for eid in face.Edges:
            self.printEdgeData(eid)
        print("-------------------------------------")
    
    def printFaceVertices(self, faceId):
        face = self.Faces[faceId]
        for vid in face.Vertices:
            self.printVertexCoords(vid)
        print("-------------------------------------")

### End of classes definition for mesh data structure


if (len(sys.argv) < 2):
    print("Error: No input vtk mesh file specified")
    quit()

mesh = Mesh()
mesh.readInputFile(sys.argv[1])
mesh.buildMeshConnectivity()
mesh.buildNormalData()
mesh.setArea()
mesh.extractPatches(90)
mesh.getFeatureVertices(5)
mesh.simplify()
mesh.smooth("AngleBased")
mesh.clean()
write_mesh(mesh, "out.vtk")

# mesh_quality.getMeshQuality(mesh, "ScaledJacobian")