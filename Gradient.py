########################################################################################################################
# Compute Gradient of a scalar function on a triangular mesh
# Author : G. Vicaigne (Internship)
# Date : 2018
########################################################################################################################

from soma import aims
from numpy import *
import pickle
import sys


def TriangleGradient(mesh, texture):
    '''
    Compute gradient on a triangular mesh with a scalar function. Gradient is computed on each triangle by the function
    v := described in http://dgd.service.tu-berlin.de/wordpress/vismathws10/2012/10/17/gradient-of-scalar-functions/.
    :param mesh: Triangular mesh
    :param texture: Scalar function on Vertices
    :type texture : AimsTextureFloat
    :return: Gradient on Triangle
    :rtype: Matrix of size number of polygons x 3
    '''

    # Initialize Parameters
    vert = array(mesh.vertex())
    poly = mesh.polygon()
    l_poly = len(poly)
    n = 0
    dicgrad = zeros([l_poly, 3])

    # Calculate the Gradient
    for i in range(l_poly):

        # Percentage done
        if int(i / float(l_poly) * 100) > n:
            n += 1
            print(str(n) + ' %')
        j = []
        for jj in poly[i]:
            j.append(jj)
        eij = [vert[j[1]][0] - vert[j[0]][0], vert[j[1]][1] - vert[j[0]][1], vert[j[1]][2] - vert[j[0]][2]]
        eki = [vert[j[0]][0] - vert[j[2]][0], vert[j[0]][1] - vert[j[2]][1], vert[j[0]][2] - vert[j[2]][2]]
        ejk = [vert[j[2]][0] - vert[j[1]][0], vert[j[2]][1] - vert[j[1]][1], vert[j[2]][2] - vert[j[1]][2]]
        A = 0.5 * linalg.norm(cross(eij, ejk))
        N = 0.5 / A * cross(ejk, eki)
        dicgrad[i] = cross(0.5 * N / A, multiply(texture[0][j[0]], ejk) + multiply(texture[0][j[1]], eki) + multiply(
                    texture[0][j[2]], eij))

    return dicgrad


def Gradient(mesh, texture):
    '''
    Compute gradient on a triangular mesh with a scalar function. Gradient is computed on each triangle by the function
    v := described in http://dgd.service.tu-berlin.de/wordpress/vismathws10/2012/10/17/gradient-of-scalar-functions/.
    On each vertex, compute the mean gradient of all triangle with the vertex.
    :param mesh: Triangular mesh
    :param texture: Scalar function on Vertices
    :type texture : AimsTextureFloat
    :return: Gradient on Vertices
    :rtype: Dictionary
    '''

    # Initialize Parameters
    vert = array(mesh.vertex())
    l_vert = len(vert)
    poly = mesh.polygon()
    l_poly = len(poly)
    n = 0

    # Initialize Dictionnary
    dicgrad = dict()
    for i in range(l_vert):
        dicgrad[i] = [0, 0, 0, 0]

    # Calculate the Gradient
    for i in range(l_poly):
        # Percentage done
        if int(i / float(l_poly) * 100) > n:
            n += 1
            print(str(n) + ' %')
        j = []
        for jj in poly[i]:
            j.append(jj)
        grad = [0., 0., 0., 0.]
        eij = [vert[j[1]][0] - vert[j[0]][0], vert[j[1]][1] - vert[j[0]][1], vert[j[1]][2] - vert[j[0]][2]]
        eki = [vert[j[0]][0] - vert[j[2]][0], vert[j[0]][1] - vert[j[2]][1], vert[j[0]][2] - vert[j[2]][2]]
        ejk = [vert[j[2]][0] - vert[j[1]][0], vert[j[2]][1] - vert[j[1]][1], vert[j[2]][2] - vert[j[1]][2]]
        A = 0.5 * linalg.norm(cross(eij, ejk))
        N = 0.5 / A * cross(ejk, eki)
        grad[0:3] = cross(0.5 * N / A, multiply(
                    texture[0][j[0]], ejk) + multiply(texture[0][j[1]], eki) + multiply(texture[0][j[2]], eij))
        grad[3] = 1.
        for jj in j:
            dicgrad[jj] = add(dicgrad[jj], grad)
    for i in range(l_vert):
        dicgrad[i] = multiply(dicgrad[i][0:3], 1 / dicgrad[i][3])
    return dicgrad


def NormGradient(mesh, texture):
    '''
    Compute the norm of a vertex Gradient on vertex
    :param mesh: Triangular mesh
    :param texture: Scalar function on Vertices
    :type texture : AimsTextureFloat
    :return: Gradient's Norm
    :rtype: Texture
    '''

    # Compute the gradient of the Mesh
    Grad = Gradient(mesh, texture)

    # Initialize Parameters
    l_grad = len(Grad)
    tex = aims.TimeTexture_FLOAT(1, l_grad)
    for i in range(l_grad):
        tex[0][i] = linalg.norm(Grad[i])
    return tex


def ArrowGradient(mesh, texture):
    '''
    Create Vertex gradient arrows Mesh for all a Mesh
    :param mesh: Triangular mesh
    :param texture: Scalar function on Vertices
    :type texture : AimsTextureFloat
    :return: Gradient arrows
    :rtype: Mesh
    '''

    # Compute the gradient of the Mesh
    Grad = Gradient(mesh,texture)

    # Initialize Parameters
    surfGen = aims.SurfaceGenerator()
    surf = aims.AimsSurfaceTriangle()
    vert = array(mesh.vertex())
    l_grad = len(Grad)
    listePoints = arange(l_grad)

    # For each gradient's vertex
    for point in listePoints:

        # Start point and Final point
        Startpoint = vert[point]
        Finalpoint = Startpoint + multiply(Grad[point], 1/linalg.norm(Grad[point]))

        # Create the Gradient arrow
        arrow = surfGen.arrow(Finalpoint, Startpoint, 0.08, 0.15, 10, 0.3)

        # Stack arrow
        surf += arrow
    return surf


def ArrowGradientDPFNegative(mesh, DPF):
    '''
    Use only on a brain mesh. Reduce the number of arrow for an easiest visibility
    Create Gradient arrows Mesh for the vertex who have a negative DPF
    :param mesh: Triangular mesh
    :param DPF: Depth Potential Function
    :type DPF : AimsTextureFloat
    :return: Gradient arrows
    :rtype: Mesh
    '''

    # Compute the gradient of the Mesh
    Grad = Gradient(mesh, DPF)

    # Initialize Parameters
    surfGen = aims.SurfaceGenerator()
    surf = aims.AimsSurfaceTriangle()
    vert = array(mesh.vertex())
    l_grad = len(Grad)
    listePoints = arange(l_grad)
    for point in listePoints:
        if DPF[0][point] < 0:
            Startpoint = vert[point]
            Finalpoint = Startpoint + multiply(Grad[point], 1 / linalg.norm(Grad[point]))
            arrow = surfGen.arrow(Finalpoint, Startpoint, 0.08, 0.15, 10, 0.3)
            surf += arrow
    return surf


def ArrowGradientDPFPositive(mesh, DPF):
    '''
    Use only on a brain mesh. Reduce the number of arrow for an easiest visibility
    Create Gradient arrows Mesh for the vertex who have a positive DPF
    :param mesh: Triangular mesh
    :param DPF: Depth Potential Function
    :type DPF : texture which fit the mesh
    :return: Gradient arrows
    :rtype: Mesh
    '''

    # Compute the gradient of the Mesh
    Grad = Gradient(mesh, DPF)

    # Initialize Parameters
    surfGen = aims.SurfaceGenerator()
    surf = aims.AimsSurfaceTriangle()
    vert = array(mesh.vertex())
    l_grad = len(Grad)
    listePoints = arange(l_grad)
    for point in listePoints:
        if DPF[0][point] >= 0:
            Startpoint = vert[point]
            Finalpoint = Startpoint + multiply(Grad[point], 1 / linalg.norm(Grad[point]))
            arrow = surfGen.arrow(Finalpoint, Startpoint, 0.08, 0.15, 10, 0.3)
            surf += arrow
    return surf


def main(path_mesh, path_scalar_function, path_gradient, gradient_type='vertex'):
    '''
    :param path_mesh: Path to the triangular mesh
    :type String
    :param path_scalar_function: Path to the Scalar Function
    :type String
    :param path_gradient: Path to save the gradient
    :type string
    :param gradient_type: Gradient on the triangle ('triangles') or Gradient on the vertex ('vertex')
    :type String
    :return:
    '''
    mesh = aims.read(path_mesh)
    scalar_func = aims.read(path_scalar_function)
    if gradient_type == 'triangles':
        grad = TriangleGradient(mesh, scalar_func)
    else:
        grad = Gradient(mesh, scalar_func)
    pickle.dump(grad, path_gradient)


if __name__ == ' __main__':
    args = sys.argv
    path_mesh = sys.argv[1]
    path_scalar_function = sys.argv[2]
    path_gradient = sys.argv[3]
    gradient_type = sys.argv[4]

    main(path_mesh, path_scalar_function, path_gradient, gradient_type=gradient_type)
