import numpy as np
import matplotlib.pyplot as plt
import random as r
from scipy.integrate import odeint
import copy

from System import *





# =-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-= #
#                                                         CELL CLASS                                                             #                                                 
# =-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-= #


class Cell():
    """ """


    def __init__(self, **traits):  
        """ Initializes the Cell object with it's embryo ID, cell ID (its position in the embryo), and its genetic network """
        
        self.genetic_network = traits.get("network", "ERROR.")
        self.network_parameters = traits.get("parameters", "ERROR.")


    def get_network(self):
        """ Returns Genetic Network """
        return self.genetic_network


    def get_network_parameters(self):
        """ Returns Genetic Network """
        return self.network_parameters


    def cell_P_log(self):
        return self.cell_P


    def get_fitness(self, fitness_map, cell_position):
        """ 
        Takes fitness map as input. 
        Takes ODE_System class object to calculate P
        calculates optimal fitness

        determines fitness
        returns fitness
        (maybe logs protein levels)

        """
        
        print("\t In get_fitness function...")


        system = ODE_System(self.genetic_network, self.network_parameters, self.morphogen_input, cell_position)         # Instantiates System
        system.construct_system()                                                                                       # Constructs System

        system.metabolize()                                                                                             # Makes system metabolize (produce P profile)

        self.cell_P = system.get_metabolism()                                                                           # Retrieve P data from system.

        optimal_fitness = fitness_map(cell_position)                                                                    # Calculate optimal fitness with fitness_map function and cell's position

        print("\n\tBack in get_fitness function...")
        print("\tCell_P: {}".format(self.cell_P))

        print("\tOptimal Fitness: {}".format(optimal_fitness))



        cell_fitness = self.cell_P * 1/(np.abs(self.cell_P - optimal_fitness))

        print("\tCell fitness: {}".format(cell_fitness))
        print("\n\tLeaving get_fitness function...")
        print("---------\n")

        return cell_fitness




    



    def morphogen_input(self, cell_position):
        G_in = np.exp(-0.05 * cell_position)
        return G_in


    def morphogen_constant(self, G, K=5, n=2):
        """ Takes in binding coefficient and n, then returns the coefficient  """

        K_G = G**n / (G**n + K**n)        

        return K_G


    def repressor_constant(self, R, K=5, n=2):
        """ 
        R(t) = 1 / (1 + (R/R_X)**n)
        """
        K_R = 1 / (1 + (R/K)**n)


        return K_R



# =-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=- END CELL CLASS -=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-= #







