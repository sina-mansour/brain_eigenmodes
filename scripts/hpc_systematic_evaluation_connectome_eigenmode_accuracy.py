# import required packages
import os
import sys
import time
import numpy as np
import nibabel as nib
from Connectome_Spatial_Smoothing import CSS as css

# variables used in the main script
main_dir = os.path.abspath('../')
subjects = np.genfromtxt(f"{main_dir}/data/pang/empirical/subject_list_HCP.txt", dtype=str)
dscalar = nib.load(css._sample_cifti_dscalar)
brain_models = [x for x in dscalar.header.get_index_map(1).brain_models]
left_cortical_surface_model = brain_models[0]
node_dimension = len(left_cortical_surface_model.vertex_indices)
rest_dimension = 1200

# Load the HCPMMP1 brain atlas
atlas_dimension = 180
HCPMMP1_labels, HCPMMP1_charmat = css.parcellation_characteristic_matrix()
left_HCPMMP1_labels = HCPMMP1_labels[:atlas_dimension]
left_HCPMMP1_charmat = HCPMMP1_charmat[:atlas_dimension, :node_dimension]
left_HCPMMP1_mean = left_HCPMMP1_charmat / left_HCPMMP1_charmat.sum(1)
left_HCPMMP1_mean = left_HCPMMP1_mean.tocsr()

# load task names
task_names = np.loadtxt(f"{main_dir}/data/task_names.txt", dtype=str)

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


def load_eigenmodes(parameters, main_dir):
    # expand parameters
    density, binary, regression, tractography, connectome_file, fwhm, global_local_epsilon = parameters

    # file names to load
    eigenmodes_file = (
        f"{main_dir}/eigenmodes/eigenmodes_density:{density}_binary:{binary}_regression:{regression}_"
        f"connectome:{tractography}_fwhm:{fwhm}_local:{global_local_epsilon}_eigenvectors.npy"
    )
    eigenvals_file = (
        f"{main_dir}/eigenmodes/eigenmodes_density:{density}_binary:{binary}_regression:{regression}_"
        f"connectome:{tractography}_fwhm:{fwhm}_local:{global_local_epsilon}_eigenvals.npy"
    )

    # load and return
    eigenmodes = np.load(eigenmodes_file)
    eigenvals = np.load(eigenvals_file)
    return eigenmodes, eigenvals


def get_dot_product_coeffs(signal, eigenmodes):
    return np.dot(signal, eigenmodes)


def compute_rest_reconstruction_coefficients_and_accuracy(eigenmodes, subject):
    # compute coefficients and accuracy for all subjects
    rest_reconstruction_coeffs = np.zeros((rest_dimension, eigenmodes.shape[1]))
    rest_reconstruction_accuracy = np.zeros((rest_dimension, eigenmodes.shape[1]))
    # load reset data
    left_rest_data = np.load(f'{main_dir}/data/hcp_data/rest/{subject}_S255.npy')
    # compute reconstruction coefficients
    rest_reconstruction_coeffs = get_dot_product_coeffs(
        left_rest_data - left_rest_data.mean(),
        eigenmodes
    )
    # compute static FC
    left_static_FC = np.corrcoef(left_HCPMMP1_mean.dot(left_rest_data.T))
    left_static_FC_triu = left_static_FC[np.triu_indices_from(left_static_FC, k=1)]
    # compute reconstruction accuracy
    rest_reconstruction_accuracy = np.array([
        np.corrcoef(
            left_static_FC_triu,
            np.corrcoef(
                left_HCPMMP1_mean.dot(
                    np.dot(
                        rest_reconstruction_coeffs[:,:length],
                        eigenmodes.T[:length]
                    ).T
                )
            )[np.triu_indices_from(left_static_FC, k=1)]
        )[0, 1]
        for length in range(1, 1 + eigenmodes.shape[1])
    ])
    return rest_reconstruction_coeffs, rest_reconstruction_accuracy


def compute_reconstruction_accuracy_atlas(eigenmodes, signal, atlas):
    # dot product
    eigenmodes_normalized = eigenmodes/np.linalg.norm(eigenmodes, axis=0)
    encoded_signal = get_dot_product_coeffs(signal.reshape(1,-1), eigenmodes_normalized)
    reconstruction_loadings = np.multiply(eigenmodes_normalized, encoded_signal)
    gradually_reconstructed_signal = np.cumsum(reconstruction_loadings, axis=1)
    # downsample to atlas
    signal_atlas = np.asarray(atlas.dot(signal)).ravel()
    gradually_reconstructed_signal_atlas = np.asarray(atlas.dot(gradually_reconstructed_signal))
    # evaluate accuracy
    reconstruction_accuracy_atlas = np.corrcoef(x=signal_atlas, y=gradually_reconstructed_signal_atlas.T)[1:,0]
    return reconstruction_accuracy_atlas

def compute_task_reconstruction_accuracy(eigenmodes, task_contrast):
    # compute coefficients and accuracy for all subjects
    task_data = np.load(f'{main_dir}/data/hcp_data/task/{task_contrast}_S255.npy')
    task_reconstruction_accuracy = np.array(
        [
            compute_reconstruction_accuracy_atlas(
                eigenmodes,
                task_data[i],
                left_HCPMMP1_mean,
            ) for i in range(task_data.shape[0])
        ]
    )
    return task_reconstruction_accuracy

if __name__ == '__main__':
    # prepare parameters from input arguments
    parameters = prepare_parameters(sys.argv)
    density, binary, regression, tractography, connectome_file, fwhm, global_local_epsilon = parameters

    # load the appropriate eigenmodes
    eigenmodes, eigenvals = load_eigenmodes(parameters, main_dir)

    # compute/load reconstruction coefficients and accuracy for rest
    for subject in subjects:
        print(f"[{time.asctime()}] Rest: subject {subject} #{list(subjects).index(subject)}")
        rest_reconstruction_coeffs_file = (
            f"{main_dir}/reconstruction_coeffs/rest:{subject}_density:{density}_"
            f"binary:{binary}_regression:{regression}_"
            f"connectome:{tractography}_fwhm:{fwhm}_local:{global_local_epsilon}.npy"
        )
        rest_reconstruction_accuracy_file = (
            f"{main_dir}/reconstruction_accuracy/rest:{subject}_density:{density}_"
            f"binary:{binary}_regression:{regression}_"
            f"connectome:{tractography}_fwhm:{fwhm}_local:{global_local_epsilon}.npy"
        )
        if file_exists(rest_reconstruction_accuracy_file):
            # load existing data
            rest_reconstruction_coeffs = np.load(rest_reconstruction_coeffs_file)
            rest_reconstruction_accuracy = np.load(rest_reconstruction_accuracy_file)
        else:
            # compute reconstruction coefficients and accuracy for rest
            (
                rest_reconstruction_coeffs, rest_reconstruction_accuracy
            ) = compute_rest_reconstruction_coefficients_and_accuracy(eigenmodes, subject)

            # save the computed data
            np.save(ensure_dir(rest_reconstruction_coeffs_file), rest_reconstruction_coeffs)
            np.save(ensure_dir(rest_reconstruction_accuracy_file), rest_reconstruction_accuracy)


    # compute/load reconstruction coefficients and accuracy for task
    for task_contrast in task_names:
        print(f"[{time.asctime()}] Task: task {task_contrast} #{list(task_names).index(task_contrast)}")
        task_reconstruction_accuracy_file = (
            f"{main_dir}/reconstruction_accuracy/task:{task_contrast}_density:{density}_"
            f"binary:{binary}_regression:{regression}_connectome:{tractography}_"
            f"fwhm:{fwhm}_local:{global_local_epsilon}.npy"
        )
        if file_exists(task_reconstruction_accuracy_file):
            # load existing data
            task_reconstruction_accuracy = np.load(task_reconstruction_accuracy_file)
        else:
            # compute the reconstruction accuracy for task
            task_reconstruction_accuracy = compute_task_reconstruction_accuracy(
                eigenmodes, task_contrast
            )
    
            # save the computed data
            np.save(ensure_dir(task_reconstruction_accuracy_file), task_reconstruction_accuracy)
