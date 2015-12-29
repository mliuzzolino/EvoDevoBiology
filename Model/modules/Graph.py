
# -----------------============----------------- #
#                  GRAPH STUFF                   #
# -----------------============----------------- #



def graph(P_log, fitness_map):
    """ """
    generation_0 = P_log[1]
    generation_x = P_log[2]

    embryo_0_gen_0 = generation_0[0][0][0]
    embryo_0_gen_x = generation_x[0][0][0]
    t = np.linspace(0, 5., 50)       # time grid

   


    fitness_mapping = []

    for i in t:
        fitness_mapping.append(fitness_map(i))

    #plt.plot(t, embryo_0_gen_0, 'g', t, embryo_0_gen_x, 'b', t, fitness_mapping, 'r')
    #plt.axis([0, 50, -5, 35])
    #plt.show()
