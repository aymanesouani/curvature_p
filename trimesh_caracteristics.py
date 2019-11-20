import numpy as np
import matplotlib.pyplot as plt

def mesh_vertex_nighbors_count(mesh):
    neig_count = []
    for neig in mesh.vertex_neighbors:
        neig_count.append(len(neig))
    return neig_count

def get_mesh_chracteristics(mesh):
    vert_neigh = mesh_vertex_nighbors_count(mesh)
    faces_area = mesh.area_faces
    faces_angles = mesh.face_angles
    edges_length = mesh.edges_unique_length
    return vert_neigh, faces_angles, faces_area, edges_length


def plot_ditributions(pop_mesh_carac, file_fig, file_fig2):
    nb_mesh = pop_mesh_carac.shape[0]
    bins = np.arange(0, 20)
    plot_vert_neigh = list()
    for p in pop_mesh_carac[:,0]:
        hist_t1, edges = np.histogram(p, bins, normed=1)
        plot_vert_neigh.append(hist_t1)
    plot_vert_neigh = np.array(plot_vert_neigh)

    X_vert_neigh = bins[:-1]
    avg_vert_neigh = np.mean(plot_vert_neigh, 0)
    std_vert_neigh = np.std(plot_vert_neigh, 0)
    Y_vert_neigh_std_sup = avg_vert_neigh + std_vert_neigh
    Y_vert_neigh_std_inf = avg_vert_neigh - std_vert_neigh
    Y_median_vert_neigh = np.median(plot_vert_neigh, 0)
    print(Y_median_vert_neigh.shape)
    print(X_vert_neigh.shape)

    bins = np.linspace(0, 3, 50)
    plot_faces_angles = list()
    for p in pop_mesh_carac[:,1]:
        hist_t1, edges = np.histogram(p.flatten(), bins, normed=1)
        plot_faces_angles.append(hist_t1)
    plot_faces_angles = np.array(plot_faces_angles)
    left, right = edges[:-1], edges[1:]
    X_angles = bins[:-1] + (bins[2] - bins[1]) / 2  # np.array([left, right]).T.flatten()
    avg_angles = np.mean(plot_faces_angles, 0)
    std_angles = np.std(plot_faces_angles, 0)
    Y_angles_std_sup = avg_angles + std_angles
    Y_angles_std_inf = avg_angles - std_angles
    Y_median_angles = np.median(plot_faces_angles, 0)
    print(Y_median_angles.shape)
    print(X_angles.shape)

    print('----------------------------------------------------')
    bins = np.linspace(0, 2, 30)
    plot_faces_area = list()
    for p in pop_mesh_carac[:,2]:
        # print(p.shape)
        hist_t1, edges = np.histogram(p.flatten(), bins, normed=1)
        plot_faces_area.append(hist_t1)
    plot_faces_area = np.array(plot_faces_area)
    left, right = edges[:-1], edges[1:]
    X_area = bins[:-1] + (bins[2] - bins[1]) / 2  # np.array([left, right]).T.flatten()
    avg_area = np.mean(plot_faces_area, 0)
    std_area = np.std(plot_faces_area, 0)
    Y_area_std_sup = avg_area + std_area
    Y_area_std_inf = avg_area - std_area
    Y_median_area = np.median(plot_faces_area, 0)
    print(Y_median_area.shape)
    print(X_area.shape)

    print('----------------------------------------------------')
    bins = np.linspace(0, 4, 50)
    plot_edge_length = list()
    for p in pop_mesh_carac[:,3]:
        # print(p.shape)
        hist_t1, edges = np.histogram(p.flatten(), bins, normed=1)
        plot_edge_length.append(hist_t1)
    plot_edge_length = np.array(plot_edge_length)
    left, right = edges[:-1], edges[1:]
    X_edge = bins[:-1] + (bins[2] - bins[1]) / 2  # np.array([left, right]).T.flatten()
    avg_edge = np.mean(plot_edge_length, 0)
    std_edge = np.std(plot_edge_length, 0)
    Y_edge_std_sup = avg_edge + std_edge
    Y_edge_std_inf = avg_edge - std_edge
    Y_median_edge = np.median(plot_edge_length, 0)
    print(Y_median_edge.shape)
    print(X_edge.shape)

    f3, axs = plt.subplots(1, 4)
    f3.suptitle('distributions across ' + str(nb_mesh) + ' meshes')
    axs[0].set_title('angles')
    axs[0].plot(np.tile(X_angles, (plot_faces_angles.shape[0], 1)).T, plot_faces_angles.T)
    # axs[0].set_xticks(X_angles)
    axs[0].grid(True)
    axs[1].set_title('area')
    axs[1].plot(np.tile(X_area, (plot_faces_area.shape[0], 1)).T, plot_faces_area.T)
    # axs[1].set_xticks(X_area)
    axs[1].grid(True)
    axs[2].set_title('neigh count')
    axs[2].plot(np.tile(X_vert_neigh, (plot_vert_neigh.shape[0], 1)).T, plot_vert_neigh.T)
    axs[2].set_xticks(X_vert_neigh)
    axs[2].grid(True)
    axs[3].set_title('edges length')
    axs[3].plot(np.tile(X_edge, (plot_edge_length.shape[0], 1)).T, plot_edge_length.T)
    # axs[3].set_xticks(X_edge)
    axs[3].grid(True)
    f3.set_size_inches(18.5, 10.5)
    plt.savefig(file_fig, bbox_inches='tight')  # , dpi=300)

    f2, axs = plt.subplots(1, 4)
    f2.suptitle('distributions across ' + str(nb_mesh) + ' meshes')
    axs[0].set_title('angles')
    axs[0].plot(X_angles, Y_angles_std_sup, '--b')
    axs[0].plot(X_angles, Y_angles_std_inf, '--b')
    axs[0].plot(X_angles, Y_median_angles, 'r')
    # axs[0].set_xticks(X_angles)
    axs[0].grid(True)
    axs[1].set_title('area')
    axs[1].plot(X_area, Y_area_std_sup, '--g')
    axs[1].plot(X_area, Y_area_std_inf, '--g')
    axs[1].plot(X_area, Y_median_area, 'r')
    # axs[1].set_xticks(X_area)
    axs[1].grid(True)
    axs[2].set_title('neigh count')
    axs[2].plot(X_vert_neigh, Y_vert_neigh_std_sup, '--g')
    axs[2].plot(X_vert_neigh, Y_vert_neigh_std_inf, '--g')
    axs[2].plot(X_vert_neigh, Y_median_vert_neigh, 'r')
    axs[2].set_xticks(X_vert_neigh)
    axs[2].grid(True)
    axs[3].set_title('egdes length')
    axs[3].plot(X_edge, Y_edge_std_sup, '--g')
    axs[3].plot(X_edge, Y_edge_std_inf, '--g')
    axs[3].plot(X_edge, Y_median_edge, 'r')
    # axs[3].set_xticks(X_edge)
    axs[3].grid(True)

    f2.set_size_inches(18.5, 10.5)
    plt.savefig(file_fig2, bbox_inches='tight')  # , dpi=300)

    # f1, axs = plt.subplots(1, 2, sharey=True)
    # f1.suptitle('FS mesh')
    # axs[0].set_title('area')
    # axs[0].boxplot(np.array(pop_faces_area), vert=False)
    # axs[0].set_yticklabels(subjects_list)
    # axs[0].grid(True)
    # axs[1].set_title('angles')
    # axs[1].boxplot(np.array(plot_faces_angles), vert=False)
    # axs[1].set_yticklabels(subjects_list)
    # axs[1].grid(True)
    # f1.set_size_inches(18.5, 10.5)
    # plt.savefig(file_fig, bbox_inches='tight')  # , dpi=300)
