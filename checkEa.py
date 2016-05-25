import os

family = ['H_Abstraction']
directory = os.path.join("/home/slakman.b/Gaussian/SMD/", family)

Ea_checkfile = open(os.path.join(directory, "Ea_check.txt"), 'w')

# Go through all of the reactions for that family
for rxn_folder in os.listdir(directory):
    rxn_directory = os.path.join(directory, rxn_folder)

    if not os.path.exists(os.path.join(rxn_directory, output.txt)):
        Ea_checkfile.write("No output file for {0}".format(str(rxn_folder)))

    else:
        with open(os.path.join(rxn_directory, "output.txt"), 'r') as output_file:
            lines = output_file.read().split('\n')
            for l in lines:
                if "Activation energy" in l:
                    Ea = float(l.split()[-1])
                    if Ea < 0:
                        Ea_checkfile.write("Ea of {0} for {1}".format(Ea, rxn_folder))

Ea_checkfile.close()
