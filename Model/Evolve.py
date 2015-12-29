import numpy as np
import matplotlib.pyplot as plt
import random as r
from scipy.integrate import odeint
import copy

from modules.Cell import *
from modules.System import *
from modules.GeneNetwork import *
from modules.FitnessMap import *
from modules.Graph import *
from modules.Fitness import *
from modules.Selection import *
from modules.Replication import *
from modules.InitializeEnvironment import *




# ---------====---------====---------====---------====---------====--------- #
#                            Initializing Constants                          #
# ---------====---------====---------====---------====---------====--------- #
number_of_generations = 5                                                    #
number_of_embryos = 10                                                       #
number_of_cells = 10                                                         #
                                                                             #
generations = {}                                                             #
                                                                             #
P_log = {}                                                                   #
                                                                             #
environment = initialize_environment(number_of_embryos, number_of_cells)     #
generations[0] = environment[:]                                              #
# ---------====---------====---------====---------====---------====--------- #







# ------------------------------============------------------------------  #
#                               MAIN PROGRAM                                #
# ------------------------------============------------------------------  #
                                                                            #
for i in range(1, number_of_generations+1):                                 #        # Iterates through generations
    print("\n\n\n\n\n\n")                                                   #
    print("=" * 20)                                                         #
    print("Generation {}".format(i))                                        #                               
    print("=" * 20)                                                         #
    print("Genetic Networks: \n")                                           #
    for embryo in environment:                                              #
        print("Embryo {}".format(environment.index(embryo)))                #
        network = embryo.get_network()                                      #
        print("Gene network:\n {}\n".format(network[0]))                    #
        print("Link network:\n {}\n".format(network[1]))                    #
        print("\n" + "-"*15)                                                #       
                                                                            #                                                         
    embryo_map = assess_fitness(environment, fitness_map, number_of_cells)  #        # ASSESS FITNESS
    embryo_fitness_map = []                                                 #        # Instantiates embyro_fitness_map
    embryo_P_map = []                                                       #        # Instantiates embryo P levels
    embryo_P_max_map = []                                                   # 
                                                                            # 
    for embryo in embryo_map:                                               #        # Iterates through the f
        embryo_fitness_map.append(embryo[0])                                #        # Logs the embryo's fitness
        embryo_P_map.append(embryo[1])                                      #        # Logs the embryo's individual cellular P levels (list)
                                                                            #
    P_log[i] = embryo_P_map                                                 #        # Logs the embryo's P levels, associating it with i, the generation number
                                                                            #
    new_environment = selection(environment, embryo_fitness_map)            #        # SELECTION - creates new environment of only half of embryos from original environment
                                                                            #
    environment = replication(new_environment)                              #        # REPLICATION - replicates new environment ** Mutation happens here
                                                                            #
    generations[i] = environment[:]                                         #        # Logs the new generation before proceeding with another iteration of evolution
                                                                            #
# ------------------------------============------------------------------  #  


# GRAPH
graph(P_log, fitness_map)

