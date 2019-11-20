import numpy as np
from matplotlib.tri import Triangulation
import trimesh



def quadric(K1, K2):
    """

    :param K1:
    :param K2:
    :return:
    """

    def fonction(x, y): return K1 * x ** 2 + K2 * y ** 2

    return fonction


def quadric_curv_gauss(K1, K2):
    """

    :param K1:
    :param K2:
    :return:
    """

    def curv_gauss(x, y):
        num = -4 * (K1 * K2)
        denom = (1 + 4 * K1 ** 2 * x ** 2 + 4 * K2 ** 2 * y ** 2) ** 2
        return num / denom

    return curv_gauss


def quadric_curv_mean(K1, K2):
    """

    :param K1:
    :param K2:
    :return:
    """

    def curv_mean(x, y):
        num = -(2 * K2 * (1 + 4 * K1 ** 2 * x ** 2) + 2 * K1 *
                (1 + 4 * K2 ** 2 * y ** 2))
        denom = 2 * (1 + 4 * K1 ** 2 * x ** 2 +
                     4 * K2 ** 2 * y ** 2) ** (3 / 2)

        return num / denom

    return curv_mean


def generate_quadric(Ks, nstep=50, ratio = 0, ax=1, ay=1):
    
    # Parameters
    
    xmin, xmax = [-ax, ax]
    ymin, ymax = [-ay, ay]
    randomSampling = True
    sigma = (2 / nstep) * ratio  # characteristic size of the mesh * ratio
    deformation_type = 'gaussian'

    # Coordinates
    x, stepx = np.linspace(xmin, xmax, nstep, retstep=True)
    y = np.linspace(ymin, ymax, nstep)
    X, Y = np.meshgrid(x, y)
    X[::2] += stepx / 2
    Y += np.sqrt(3) / 2
    X = X.flatten()
    Y = Y.flatten()


    if randomSampling:

        if deformation_type == ('gamma'):
            theta = np.random.rand(nstep * nstep, ) * np.pi * 2
            mean = sigma
            variance = sigma ** 2
            radius = np.random.gamma(mean ** 2 / variance, variance / mean, nstep * nstep)
            X = X + radius * np.cos(theta)
            Y = Y + radius * np.sin(theta)
        elif deformation_type == ('uniform'):
            X = X + np.random.uniform(-1, 1, 100)
            Y = Y + np.random.uniform(-1, 1, 100)
        else:
            X = X + sigma * np.random.randn(nstep * nstep, )
            Y = Y + sigma * np.random.randn(nstep * nstep, )

    # fig, ax = plt.subplots()
    # ax.plot(X, Y, '+')
    # plt.show()

    # Delaunay triangulatio
    Tri = Triangulation(X, Y)

    Zs = list()
    for K in Ks:
        Z = quadric(K[0], K[1])(X, Y)
        Zs.append(Z)

    return X, Y, Tri.triangles, Zs


def generate_ellipsiod(a, b, nstep, randomSampling):
    # Coordinates
    if randomSampling:
        THETA = (np.random.rand(nstep * nstep, 1) - 1 / 2) * np.pi
        PHI = 2 * np.pi * np.random.rand(nstep * nstep, 1)
    else:
        theta = np.linspace(-np.pi / 2, np.pi / 2, nstep)
        phi = np.linspace(0, 2 * np.pi, nstep)
        THETA, PHI = np.meshgrid(theta, phi)

    # Sphere coordinates
    X = a * np.cos(THETA) * np.cos(PHI)
    Y = b * np.cos(THETA) * np.sin(PHI)
    Z = np.sin(THETA)

    coords = np.array([X, Y, Z]).squeeze().transpose()

    return tri_from_hull(coords)


def tri_from_hull(vertices):
    """
    compute faces from vertices using trimesh convex hull
    :param vertices: (n, 3) float
    :return:
    """
    mesh = trimesh.Trimesh(vertices=vertices, process=False)
    return mesh.convex_hull


def generate_sphere(n=100):
    """

    :param n:
    :return:
    """
    coords = np.zeros((n, 3))
    for i in range(n):
        M = np.random.normal(size=(3, 3))
        Q, R = np.linalg.qr(M)
        coords[i, :] = Q[:, 0].transpose() * np.sign(R[0, 0])
    return tri_from_hull(coords)
