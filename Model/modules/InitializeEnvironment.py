from GeneNetwork import *
from Cell import *
# -------------------------------=============------------------------------- #
#                             INITIALIZE ENVIRONMENT                          #
# -------------------------------=============------------------------------- #

def initialize_environment(embryos, cells):
    """ Initializes Generation 0 with homogeneous, simple genetic networks """

    #new_network = Genetic_Network(initialize=True)

    initial_environment = []
    for i in range(0, embryos):
        new_network = Genetic_Network(initialize=True)
        
        new_traits = {"network": new_network.get_network(), "parameters": new_network.get_parameters()}
        new_cell = Cell(**new_traits)
        initial_environment.append(new_cell)

    
    return initial_environment
# -------------------------------=============------------------------------- #

