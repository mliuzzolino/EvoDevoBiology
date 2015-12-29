

# -------------------------------=============------------------------------- #
#                                   FITNESS                                   #
# -------------------------------=============------------------------------- #

def assess_fitness(environment, fitness_map, number_of_cells):
    """  """
    print("\n")
    print("-" * 80)
    print("\t\t\tFitness Assessment")
    print("-" * 80)
    print("\n")

    print("Environment: {}".format(environment))

    
    

    
    embryo_fitness = []
    embryo_P_map = []
    embryo_map = []
    j = 0
    for embryo in environment:                                      # Iterates through embryo in environment
        print("--------------\n")
        print("Embryo {}".format(j))
        print("Embryo id: {}".format(id(embryo)))

        cell_fitness = []                                           # Initializes cell_fitness list to log fitness of each cell of embryo
        cell_P = []                                                 # Initializes cell_P list to log protein expression of each cell of embryo

        for i in range(0, number_of_cells):                         # Iterates through each cell in the embryo
            print("\nCell {}. Assessing Fitness...".format(i))
            cell_position = i                                       # Determines cell position in embryo
            fitness = embryo.get_fitness(fitness_map, cell_position)
            cell_fitness.append(fitness)      # calls get_fitness method of cell with fitness_map. Appends fitness to list.
            cell_P.append(embryo.cell_P_log)


        embryo_fitness = sum(cell_fitness)                          # Sums the cellular fitness to represent embryo fitness.
        print("\nEmbryo fitness: {}".format(embryo_fitness) )
        embryo_P_map.append(cell_P)                                 # Appends cellular P levels to embryo_P_map


        embryo_map.append([embryo_fitness, embryo_P_map])           # Appends fitness and P levels as list to embryo_map
        j += 1

    return embryo_map


# -------------------------------=============------------------------------- #




