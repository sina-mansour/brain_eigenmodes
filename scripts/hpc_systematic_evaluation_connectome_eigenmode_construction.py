# This script computes connectome eigenmodes according to the input parameters
# Usage:
#   python3 ./hpc_systematic_evaluation_connectome_eigenmode_construction.py \
#           <density> \
#           <binary|weighted> \
#           <with_gyral_regression|without_gyral_regression> \
#           <pang|optimized> \
#           <fwhm> \
#           <global_local_epsilon>
# 
# 
# Example: To replicate Pang et al. the following parameters can be used:
#   python3 ./hpc_systematic_evaluation_connectome_eigenmode_construction.py 0.001 binary without_gyral_regression pang 0 1
# 

# import required packages
import os
import sys
import numpy as np
import scipy.sparse as sparse
from sklearn import linear_model
import nibabel as nib
from Connectome_Spatial_Smoothing import CSS as css

# variables used in the main script
main_dir = os.path.abspath('../')
dscalar = nib.load(css._sample_cifti_dscalar)
brain_models = [x for x in dscalar.header.get_index_map(1).brain_models]
left_cortical_surface_model = brain_models[0]
node_dimension = len(left_cortical_surface_model.vertex_indices)
local_adjacency_file = f"{main_dir}/data/templates/local_adjacency.npz"
left_local_adjacency = sparse.load_npz(local_adjacency_file).tocsr()[:node_dimension, :node_dimension]

# functions used in the main script
def ensure_dir(file_name):
    os.makedirs(os.path.dirname(file_name), exist_ok=True)
    return file_name


def file_exists(file_name, path_name=os.getcwd()):
    return os.path.isfile(os.path.join(path_name, file_name))


def prepare_parameters(arguments):
    # read input arguments from command line
    density, binary, regression, tractography, fwhm, global_local_epsilon = arguments[1:]

    # change density to float
    density = float(density)

    # change binary parameter to a boolean flag
    if binary == "binary":
        binary = True
    elif binary == "weighted":
        binary = False
    else:
        raise ValueError("Incorrect option provided for the binary parameter.")

    # change regression parameter to a boolean flag
    if regression == "with_gyral_regression":
        regression = True
    elif regression == "without_gyral_regression":
        regression = False
    else:
        raise ValueError("Incorrect option provided for the regression parameter.")

    # locate the high-resolution connectome file
    if tractography == "pang":
        connectome_file = f"{main_dir}/data/connectomes/pang_group_average_connectome.npz"
    elif tractography == "ours":
        connectome_file = f"{main_dir}/data/connectomes/our_group_average_connectome.npz"
    else:
        raise ValueError("Incorrect option provided for the tractography parameter.")

    # change fwhm to float
    fwhm = float(fwhm)

    # change global_local_epsilon to float
    global_local_epsilon = float(global_local_epsilon)

    return (density, binary, regression, tractography, connectome_file, fwhm, global_local_epsilon)


def regress_gyral_bias(connectome):
    # load curvature information
    curvature = nib.load(f"{main_dir}/data/templates/S1200.curvature_MSMAll.32k_fs_LR.dscalar.nii").get_fdata()[0]

    # get left hemisphere curvature
    curvature = curvature[:node_dimension].reshape(-1, 1)

    # compute node strengths and degrees
    node_strengths = np.array(connectome.sum(0)).reshape(-1, 1)
    node_degrees = np.array((connectome > 0).sum(0)).reshape(-1, 1)

    # fit a linear regressions
    gyral_bias_predictions = linear_model.LinearRegression().fit(X=curvature, y=node_strengths).predict(curvature)

    # remove the gyral bias effect from connectome
    connectome_curvature_adjusted = connectome.tocoo()
    indices = connectome_curvature_adjusted.row
    connectome_curvature_adjusted.data -= 2 * (gyral_bias_predictions[indices, 0] / node_degrees[indices, 0])

    # add constant to ensure connection weights are not negative
    connectome_curvature_adjusted.data += 2 * (gyral_bias_predictions.max() / node_degrees[indices, 0])
    connectome_curvature_adjusted = connectome_curvature_adjusted.tocsr()

    # symmetric adjustment (to account for bias at both ends)
    connectome_curvature_adjusted = (connectome_curvature_adjusted + connectome_curvature_adjusted.T)/2

    # return adjusted connectome
    return connectome_curvature_adjusted


def threshold_connectome(connectome, threshold):
    connectome = connectome.multiply(
        connectome > threshold
    )
    connectome.eliminate_zeros()
    return connectome


def connectome_smoothing(connectome, fwhm):
    # Group-level surfaces to compute geodesic distance
    left_white_surface_file = f'{main_dir}/data/templates/S1200.L.white_MSMAll.32k_fs_LR.surf.gii'
    right_white_surface_file = f'{main_dir}/data/templates/S1200.R.white_MSMAll.32k_fs_LR.surf.gii'

    # CSS smoothing parameters
    sigma = css._fwhm2sigma(fwhm)
    kernel_radius = 6

    # Local geodesic distance computed on the cortical surfaces
    cortical_local_geodesic_distances = css._get_cortical_local_distances(
        left_white_surface_file,
        right_white_surface_file,
        kernel_radius,
        css._sample_cifti_dscalar
    )

    # keep only the left surface distances
    left_cortical_local_geodesic_distances = (
        cortical_local_geodesic_distances[:node_dimension,:][:,:node_dimension]
    )

    # convert to smoothing coefficients
    left_smoothing_kernel = css._local_distances_to_smoothing_coefficients(
        left_cortical_local_geodesic_distances,
        sigma
    )

    # apply smoothing
    connectome = threshold_connectome(
        css.smooth_high_resolution_connectome(
            connectome,
            left_smoothing_kernel
        ),
        1e-4
    )

    return connectome


def edge_count(connectome):
    connectome.eliminate_zeros()
    return (connectome.data.shape[0] - (connectome.diagonal() > 0).sum()) / 2


def find_density_threshold(connectome, density, precision=1e-5, max_iterations=50):
    # edge count to meet the desired density
    full_edge_count = (np.prod(connectome.shape) - connectome.shape[0]) * 0.5
    goal_edge_count = density * full_edge_count
    edge_count_tolerance = precision * full_edge_count

    # define a threshold range to search over
    search_range = (connectome.data.min() / 1.01, connectome.data.max())
    edge_count_range = (edge_count(connectome), 0)

    # exit if desired density is greater than the connectome density
    if goal_edge_count > edge_count_range[0]:
        print(
            f'Warning, connectome density ({edge_count_range[0] / full_edge_count}) '
            f'is lower than desired density ({density}).'
        )
        return search_range[0]

    # binary search with tolerance to find the threshold
    old_threshold, old_edge_count = search_range[1], edge_count_range[1]
    iterations = 0
    while iterations < max_iterations:
        iterations += 1
        new_threshold = np.mean(search_range)
        new_edge_count = edge_count(connectome > new_threshold)

        # check if threshold is good enough
        if abs(new_edge_count - goal_edge_count) < edge_count_tolerance:
            break

        elif new_edge_count > goal_edge_count:
            search_range = (new_threshold, search_range[1])
            edge_count_range = (new_edge_count, edge_count_range[1])

        else:
            search_range = (search_range[0], new_threshold)
            edge_count_range = (edge_count_range[0], new_edge_count)

        if iterations >= max_iterations:
            print(
                f'Warning, max iterations reached; goal density: {density}, '
                f'achieved density: {new_edge_count / full_edge_count}'
            )

    print(f"iterations: {iterations}, achieved density: {new_edge_count / full_edge_count}")
    return new_threshold


def get_connectome(parameters, main_dir):
    # expand parameters
    density, binary, regression, tractography, raw_connectome_file, fwhm, global_local_epsilon = parameters

    # file name to store/load connectome before density threshold
    connectome_file = (
        f"{main_dir}/connectomes/connectome:{tractography}_regression:{regression}_"
        f"fwhm:{fwhm}.npz"
    )

    # load the computed connectome (before density threshold)
    if file_exists(connectome_file):
        connectome = sparse.load_npz(connectome_file)
    # compute the connectome if not already existing
    else:
        # load the source connectome file
        connectome = sparse.load_npz(raw_connectome_file)

        # gyral bias regression
        if regression:
            connectome = regress_gyral_bias(connectome)

        # perform smoothing (~5mins)
        if (fwhm > 0):
            connectome = connectome_smoothing(connectome, fwhm)

        # save the computed connectome
        sparse.save_npz(ensure_dir(connectome_file), connectome)

    # apply density threshold
    if (density < 1):
        # find density threshold (~1min)
        threshold = find_density_threshold(connectome, density)

        # adjust density
        connectome = threshold_connectome(connectome, threshold)

    # combine with local
    connectome = connectome + (global_local_epsilon * left_local_adjacency)

    # binarize
    if binary:
        connectome = (connectome > 0).astype(float)

    return connectome


# function to compute the symmetric normalized adjacency
def symmetric_normalized_adj(adj):
    degree = sparse.diags(np.array(adj.sum(axis=0))[0])
    inverse_degree_square_root = degree.copy()
    inverse_degree_square_root.data **= -0.5
    SNA = (inverse_degree_square_root*adj*inverse_degree_square_root)
    return SNA


# function performing a symmetric Laplacian eigenmode decomposition over a sparse adjacency structure
def laplace_spectrum_rev(adj, k):
    # use shift invert mode
    lambdas, vectors = sparse.linalg.eigsh(
        symmetric_normalized_adj(adj),
        k=k+1,
        which="LM",
    )

    # convert to laplacian lambdas
    lambdas = 1 - lambdas
    lambda_idx = np.argsort(lambdas)
    lambdas = lambdas[lambda_idx]
    vectors = vectors[:, lambda_idx]
    return (lambdas, vectors)


def compute_eigenmodes(parameters, main_dir, overwrite=False):
    # expand parameters
    density, binary, regression, tractography, connectome_file, fwhm, global_local_epsilon = parameters

    # file name to store/load eigenmodes from
    eigenmodes_file = (
        f"{main_dir}/eigenmodes/eigenmodes_density:{density}_binary:{binary}_regression:{regression}_"
        f"connectome:{tractography}_fwhm:{fwhm}_local:{global_local_epsilon}_eigenvectors.npy"
    )
    eigenvals_file = (
        f"{main_dir}/eigenmodes/eigenmodes_density:{density}_binary:{binary}_regression:{regression}_"
        f"connectome:{tractography}_fwhm:{fwhm}_local:{global_local_epsilon}_eigenvals.npy"
    )

    # skip running if already existing and overwrite is not requested
    if ((not overwrite) and file_exists(eigenmodes_file)):
        pass
    # compute the eigenmodes if not already existing or requested to overwrite
    else:
        # get the connectome that satisfies the parameters
        connectome = get_connectome(parameters, main_dir)

        # compute the first 200 eigenmodes (~2min)
        eigenvalues, eigenvectors = laplace_spectrum_rev(connectome, k=200)

        # store the eigenmodes and eigenvectors
        np.save(ensure_dir(eigenmodes_file), eigenvectors)
        np.save(ensure_dir(eigenvals_file), eigenvalues)


if __name__ == '__main__':
    # prepare parameters from input arguments
    parameters = prepare_parameters(sys.argv)

    # compute the eigenmodes according to the parameters
    compute_eigenmodes(parameters=parameters, main_dir=main_dir)
