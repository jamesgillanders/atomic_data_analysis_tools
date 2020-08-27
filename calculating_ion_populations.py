"""This script can pull the level information for any species the user is
interested in from NIST, and compute the partition function values for each
level. This, combined with some line transition values the user needs to supply
(wavelength, upper level g, upper level Energy, A-value) the program can
estimate the number of atoms/ions in the upper energy level needed to produce
the feature of interest. Additionally, it can then determine the total number
of atoms/ions there are across all levels (assuming LTE). If the user doesn't
want that, at its simplest, the code can be used as a simple partition function
calculator"""


"""loading in the relevant packages"""
import numpy as np
import pandas as pd
import os


def _mkdir(newdir):
    """this function is for creating the required directory for the output files
    Acquired from: http://code.activestate.com/recipes/82465-a-friendly-mkdir/"""
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError(
            "a file with the same name as the desired "
            "dir, '%s', already exists." % newdir
        )
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        if tail:
            os.mkdir(newdir)
    return newdir


"""these are the things the user has to input"""
element = "Nb"
ion = "I"
atomic_munber = 41
wavelength = 20378 * 1e-10  #  Angstroms --> m
g_upper = (1.5 * 2) + 1
E_upper = 5297.9 * 0.00012398 * 1.6e-19  # cm^-1 --> J
flux_in_feature = 3.21e32  # J s^-1
A_value = 0.031  # s^-1

"""temperature values to consider for the calculation"""
Temp_values = [500, 1000, 2000, 3000, 4000, 5000, 10000]

output_folder_name = _mkdir(f"output")

number_photons = flux_in_feature / (6.63e-34 * 3e8 / wavelength)
number_ions_in_upper_level = number_photons / A_value

"""defining lists we will populate later"""
list_partition_functions = []
list_ions_in_upper_level = []
list_ions_in_upper_level_solar_masses = []
list_total_number_of_ions = []
list_total_number_of_ions_solar_masses = []

"""reads the relevant level info from NIST into a dataframe"""
species_df = pd.read_table(
    f"https://physics.nist.gov/cgi-bin/ASD/energy1.pl?encodedlist=XXT2&de=0&spectrum={element}++{ion}&units=1&upper_limit=&parity_limit=both&conf_limit=All&conf_limit_begin=&conf_limit_end=&term_limit=All&term_limit_begin=&term_limit_end=&J_limit=&format=3&output=0&page_size=15&multiplet_ordered=0&conf_out=on&term_out=on&level_out=on&unc_out=on&g_out=on&biblio=on&temp=&submit=Retrieve+Data",
    sep="\t",
    header=0,
    names=["Configuration", "Term", "g", "Level(eV)", "Uncertainty(eV)"],
    index_col=False,
)

"""this strips all unwanted characters from the level energies column"""
for aa in range(len(species_df)):
    species_df["Level(eV)"][aa] = "".join(
        [c for c in str(species_df["Level(eV)"][aa]) if c in "1234567890."]
    )
    if species_df["Level(eV)"][aa] == "":
        species_df["Level(eV)"][aa] = "NaN"

"""casts the columns of interest as floats"""
species_df["g"] = species_df["g"].astype(float)
species_df["Level(eV)"] = species_df["Level(eV)"].astype(float)

"""adds extra column of energy levels in Joules"""
species_df["Level(J)"] = species_df["Level(eV)"] * 1.6e-19

for Temp in Temp_values:
    """calculates the partition function for each level"""
    partition_function = species_df["g"] * np.exp(
        (-1 * species_df["Level(J)"]) / (1.38e-23 * Temp)
    )
    """adds up all the partition function values to get the sum across the entire ion/atom for all levels"""
    sum_partition_function = partition_function.sum()
    """calculates the number of ions/atoms (assuming LTE)"""
    number_total_ions_in_all_levels = (
        number_ions_in_upper_level
        * (sum_partition_function / g_upper)
        * np.exp(E_upper / (1.38e-23 * Temp))
    )

    """populating the relevant lists with values"""
    list_partition_functions.append(sum_partition_function)
    list_ions_in_upper_level.append(number_ions_in_upper_level)
    list_ions_in_upper_level_solar_masses.append(
        number_ions_in_upper_level * atomic_munber * 1.67e-27 / 2e30
    )
    list_total_number_of_ions.append(number_total_ions_in_all_levels)
    list_total_number_of_ions_solar_masses.append(
        number_total_ions_in_all_levels * atomic_munber * 1.67e-27 / 2e30
    )

"""bundling all the lists into a dataframe"""
dataframe = pd.DataFrame(
    {
        "Temperature(K)": Temp_values,
        "Partition_function(Z)": list_partition_functions,
        "Wavelength(m)": wavelength,
        "Num_ions_upper_lvl": list_ions_in_upper_level,
        "Mass_ions_upper_lvl(m_sol)": list_ions_in_upper_level_solar_masses,
        "Num_of_total_ions": list_total_number_of_ions,
        "Mass_of_total_ions(m_sol)": list_total_number_of_ions_solar_masses,
    }
)

"""printing_dataframe to file"""
dataframe.to_csv(
    f"{output_folder_name}/{atomic_munber}-{element}{ion}_mass_estimates.csv",
    index=False,
)
