import copy
import numpy as np
import random as r
from GeneNetwork import *
from Cell import *

# -------------------------------=============------------------------------- #
#                                 REPLICATION                                 #
# -------------------------------=============------------------------------- #
def replication(environment):
    """ """
    print("\n")
    print("-" * 80)
    print("\t\t\tReplication")
    print("-" * 80)
    print("\n")

    new_environment = []



    for embryo in environment:
        print("Embryo: {}".format(embryo))
        
        embryo_network = embryo.get_network()
        embryo_network_parameters = embryo.get_network_parameters()     # FOR NOW, no mutations here. Just pass them along.

        embryo_network_copy_1 = copy.deepcopy(embryo_network)
        embryo_network_copy_2 = copy.deepcopy(embryo_network)

        embryo_network_copies = [embryo_network_copy_1, embryo_network_copy_2]


        i = 0
        for network in embryo_network_copies:
            mutate = r.sample([True, False], 1)[0]
            
            print("\n--------------")
            print("Offspring {}".format(i+1))
            if mutate:
                print("**Mutation**")
                new_network = Genetic_Network(network)
                new_network.mutate()
                traits = {"network": new_network.get_network(), "parameters": embryo_network_parameters}
                new_cell = Cell(**traits)
                new_environment.append(new_cell)

            else:
                print("No Mutation")
                traits = {"network": network, "parameters": embryo_network_parameters}
                new_cell = Cell(**traits)
                new_environment.append(new_cell)
            i += 1

            
    return new_environment

# -------------------------------=============------------------------------- #



