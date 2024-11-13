import os
import re
import numpy as np
import matplotlib.pyplot as plt
import sys
from dqmc_analysis_tools import Get_den_orb

def parse_filename(filename):
    """
    Extract relevant parameters (excluding 's' and 'mu') from the filename 
    for use in plot title and filename.
    """
    filename = filename.replace('.out', '')
    
    # Use regular expressions to find parameters excluding 's1234567' and 'mu'
    params = re.findall(r'(Ud\d+\.\d+|Up\d+\.\d+|tpd\d+\.\d+|tpp\d+\.\d+|ep\d+\.\d+|N\d+|be\d+\.\d+)', filename)
    title_str = ", ".join(params)  # Title for the plot
    filename_str = "_".join(params)  # Name for the saved file
    
    return title_str, filename_str

def extract_avg_sign(file):
    """
    Extract the average sign from the file content based on the line 
    containing 'Avg sign :' followed by the average sign value.
    """
    avg_sign = None
    with open(file, 'r') as f:
        for line in f:
            if 'Avg sign :' in line:
                try:
                    # Extract the number following 'Avg sign :'
                    avg_sign = float(line.split('Avg sign :')[1].strip().split()[0])
                except (IndexError, ValueError):
                    print("Error parsing Avg sign in file: {}".format(file))
                break
    return avg_sign

def analyze_dqmc_data(file_list, norb, Nline):
    """
    Process a list of DQMC output files to extract chemical potential values (mu), 
    density values for each orbital, and average sign values.
    """
    mu_values = []
    densities = np.zeros((len(file_list), norb))  # Initialize array to store density data
    avg_signs = []  # Initialize list to store average sign data
    print("Processing files:")
    for file in file_list:
        print(file)    
    for i, file in enumerate(file_list):
        # Extract the mu value from the filename
        try:
            mu_str = re.search(r'mu(-?\d+(\.\d+)?)', file).group(1)
            mu = float(mu_str)
            mu_values.append(mu)
        except (AttributeError, ValueError):
            print("Filename format is incorrect, skipping file: {}".format(file))
            continue

        # Use Get_den_orb function to get density
        dens = Get_den_orb(file, norb, Nline)
        
        # Store the densities for each orbital in the appropriate row
        for j in range(norb):
            densities[i, j] = dens[j, 0]  # Take the first column (density)
        
        # Extract the average sign from the file
        avg_sign = extract_avg_sign(file)
        if avg_sign is not None:
            avg_signs.append(avg_sign)
    
    # Sort by mu values
    sorted_indices = np.argsort(mu_values)
    mu_values = np.array(mu_values)[sorted_indices]
    densities = densities[sorted_indices]
    avg_signs = np.array(avg_signs)[sorted_indices]  # Sort avg_signs according to mu

    print("mu values, densities, and avg signs:")
    for mu, density, sign in zip(mu_values, densities, avg_signs):
        print("mu = {}, Density = {}, Avg sign = {}".format(mu, density, sign))

    return mu_values, densities, avg_signs

# Main function
if __name__ == "__main__":
    # Get file names from command-line arguments
    if len(sys.argv) < 2:
        print("Please provide at least one filename as an argument.")
        sys.exit(1)
    
    file_list = sys.argv[1:]  # Get the list of files from command-line arguments
    file_list = [file for file in file_list if not file.endswith(".tdm.out")]    
    # Check if files exist
    for file in file_list:
        if not file.endswith(".out") or not os.path.isfile(file):
            print("Invalid file or file does not exist: {}".format(file))
            sys.exit(1)

    norb = 3  # Number of orbitals. Cu 0 Ox 1 Oy 2
    Nline = 75  # Number of lines to read. 160 for Ncell 36;75 for Ncell 16

    # Extract title and filename strings from the first file in the list
    title_str, filename_str = parse_filename(file_list[0])

    # Initialize lists to hold the mu values, densities, and avg_signs
    mu_values, densities, avg_signs = analyze_dqmc_data(file_list, norb, Nline)

    # Plot Density vs Chemical Potential mu
    plt.figure(figsize=(8, 6))
    orbital_names = ['Cu', 'Ox', 'Oy']

    for i in range(norb):
        plt.plot(mu_values, densities[:, i], label=f'Orbital {orbital_names[i]}', marker='o')
    total_density = densities.sum(axis=1)
    plt.plot(mu_values, total_density, label='Total', marker='o')
    plt.xlabel(r'$\mu$')
    plt.ylabel('Density')
    plt.title(f'Density vs Chemical Potential μ\n({title_str})')
    plt.legend()
    plt.grid(True)

    # Ensure the results directory exists
    if not os.path.exists('./results'):
        os.makedirs('./results')

    # Save Density vs mu plot
    output_filename_density = f"./results/Density_vs_mu_{filename_str}.pdf"
    plt.savefig(output_filename_density, format='pdf')
    plt.close()

    # Plot Avg Sign vs Chemical Potential μ
    plt.figure(figsize=(8, 6))
    plt.plot(mu_values, avg_signs, label='Average Sign', marker='o', color='purple')

    plt.xlabel(r'$\mu$')
    plt.ylabel('Average Sign')
    plt.title(f'Average Sign vs Chemical Potential μ\n({title_str})')
    plt.legend()
    plt.grid(True)

    # Save Avg Sign vs mu plot
    output_filename_sign = f"./results/Avg_Sign_vs_mu_{filename_str}.pdf"
    plt.savefig(output_filename_sign, format='pdf')
    plt.close()

