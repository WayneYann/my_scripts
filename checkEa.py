import os

family = 'H_Abstraction'
directory = os.path.join("/Users/belinda/Gaussian/SMD/", family)

Ea_checkfile = open(os.path.join(directory, "Ea_check.txt"), 'w')

total = 0
unsuccessful = 0
negative = 0

# Go through all of the reactions for that family
for rxn_folder in os.listdir(directory):
    if os.path.isdir(os.path.join(directory,rxn_folder)):
        total += 1
        rxn_directory = os.path.join(directory, rxn_folder)

        if not os.path.exists(os.path.join(rxn_directory, "output.txt")):
            Ea_checkfile.write("No output file for {0}\n".format(str(rxn_folder)))
            unsuccessful += 1
        else:
            with open(os.path.join(rxn_directory, "output.txt"), 'r') as output_file:
                lines = output_file.read().split('\n')
                for l in lines:
                    if "Activation energy" in l:
                        Ea = float(l.split()[-1])
                        if Ea < 0:
                            negative += 1
                            Ea_checkfile.write("Ea of {0} for {1}\n".format(Ea, rxn_folder))
Ea_checkfile.write("total: {0}, unsuccessful = {1}, negative = {2}".format(total, unsuccessful, negative))
Ea_checkfile.close()
