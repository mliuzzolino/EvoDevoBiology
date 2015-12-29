





# -------------------------------=============------------------------------- #
#                                  SELECTION                                  #
# -------------------------------=============------------------------------- #

def selection(environment, embryo_fitness_map):
    """ """

    print("\n")
    print("-" * 80)
    print("\t\t\tSelection")
    print("-" * 80)
    print("\n")

    print("Currently Undergoing Selection...\n\n")


    print("embryo_fitness_map: {}".format(embryo_fitness_map))
    selected_environment = []
    survival_index = []

    fitness_map = embryo_fitness_map[:]
    print("fitness_map: {}".format(fitness_map))


    for i in range(0, len(fitness_map)/2):
        maximum = max(fitness_map)
        max_index = fitness_map.index(maximum)
        survival_index.append(max_index)
        fitness_map.pop(max_index)
        print("max: {}".format(maximum))

    for embryo in survival_index:
        selected_environment.append(environment[embryo])

    return selected_environment
# -------------------------------=============------------------------------- #


