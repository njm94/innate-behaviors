import argparse
import scipy.io
import numpy as np

def convert_mat_to_npy(input_file, output_file, variable_name):
    # Load the .mat file
    mat_data = scipy.io.loadmat(input_file)

    # Extract the specified variable from the .mat file
    if variable_name in mat_data:
        data = mat_data[variable_name]
    else:
        raise ValueError(f"Variable '{variable_name}' not found in the .mat file.")

    # Save the extracted data as a .npy file
    np.save(output_file, data)

def main():
    # Set up the argument parser to handle command-line arguments
    parser = argparse.ArgumentParser(description='Convert a .mat file to a .npy file.')

    # Add arguments for the input .mat file, output .npy file, and variable name
    parser.add_argument('input_file', type=str, help='Path to the input .mat file')
    parser.add_argument('output_file', type=str, help='Path to the output .npy file')
    parser.add_argument('variable_name', type=str, help='Name of the variable to extract and save')

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the conversion function with the provided arguments
    convert_mat_to_npy(args.input_file, args.output_file, args.variable_name)

# Ensure that the main function runs when this script is executed directly
if __name__ == '__main__':
    main()

# How to run this script:
# 1. Save this script as `mat_to_npy.py`.
# 2. Open your terminal or command prompt.
# 3. Navigate to the directory where this script is saved.
# 4. Run the script using the following command:
#    python mat_to_npy.py input_file.mat output_file.npy variable_name
#    - Replace `input_file.mat` with the path to your .mat file.
#    - Replace `output_file.npy` with the desired path for the output .npy file.
#    - Replace `variable_name` with the name of the variable in the .mat file you want to save as .npy.
