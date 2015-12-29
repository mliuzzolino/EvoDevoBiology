
import numpy as np
import matplotlib.pyplot as plt
import random as r
from scipy.integrate import odeint
import copy

# =-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-= #
#                                                    GENETIC NETWORK CLASS                                                       #                                                 
# =-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-= #

class Genetic_Network():
    """ """

    def __init__(self, network=np.nan, parameters=np.nan, initialize=False):
        """ """

        if initialize == True:
            self.S = 100
            self.n = 2
            self.K = 1
            self.delta = 2

            self.gene_parameters  = {"S": self.S, "n": self.n, "K": self.K, "delta": self.delta}
            self.link_parameters = {"n": self.n, "K": self.K}

            # Initialize gene/link networks
            initialize_gene_network = np.zeros(10).reshape(10, 1)
            initialize_link_network = np.zeros(100).reshape(10, 10)

            initialize_gene_network[0][0] = 1

            # Initialize parameter networks
            initialize_gene_parameters = np.zeros(10)
            initialize_gene_parameter_network = list(initialize_gene_parameters)

            initialize_link_parameters = np.zeros(10)
            initialize_link_parameter_network = []
            for i in range(0, len(initialize_link_parameters)):
                initialize_link_parameter_network.append(list(initialize_link_parameters))

            for i in range(0, len(initialize_gene_parameter_network)):
                initialize_gene_parameter_network[i] = self.gene_parameters

            for i in range(0, len(initialize_link_parameter_network)):
                for j in range(0, len(initialize_link_parameter_network)):
                    initialize_link_parameter_network[i][j] = self.link_parameters


            
            self.network = [initialize_gene_network, initialize_link_network]
            self.parameters = [initialize_gene_parameter_network, initialize_link_parameter_network]
            

        else: 
            self.network = network
            self.parameters = parameters


    def get_network(self):
        return self.network

        
    def get_parameters(self):
        return self.parameters


    def mutate(self):
        """ 
        1) Genes can be turned on or off
        2) Linkages can be turned on or off
        3) Parameters (hill coeff, etc) can be mutated.
        """

        print("\n\nIn Mutation Method...\n")

        network = self.get_network()
        gene_network = network[0]
        link_network = network[1]

        print("gene network pre mutation:\n {}\n".format(gene_network))
        print("link network pre mutation:\n {}\n".format(link_network))
        #print("network id: {}".format(id(network)))



        # Determine if Genes, Linkages, or Parameters will be mutated:
        mutations = [0, 1]                                  # Key: mutations = [0, 1, 2] where 0 is genes, 1 is links, and 2 is parameters.
        mutation = r.sample(mutations, 1)[0]                # Determines mutation type

        mutate_states = [0, 1, 1, 1]                        # Mutation turns gene/link on/off: 0 for off, 1 for on
        

        # Determine which genes are currently turned on
        gene_index = []
        for i in range(0, len(gene_network)):
            if gene_network[i] == 1:
                gene_index.append(i)


        # GENES
        if mutation == 0:
            print("Mutation in Genes")

            # Determines if current gene is turned off or if a new gene is turned on.
            mutate_state = r.sample(mutate_states, 1)[0]  # Determines if turned on or off
            

            if mutate_state == 0:
                print("Existing Gene Turned Off")
                
                #try: 
                gene_mutate = r.sample(gene_index, 1)[0]            # determines which gene (amongst those currently ON) will be turned OFF
                if gene_mutate == 0:
                    gene_mutate = r.sample(gene_index, 1)[0]        # provides another chance for P to not be turned off
                    if gene_mutate == 0:
                        gene_mutate = r.sample(gene_index, 1)[0]    # provides another chance for P to not be turned off
            

                for j in range(0, len(gene_network)):               # Turns off all links associated with gene
                    gene_network[gene_mutate][j] = 0       
            
                #except:
                    #print("All genes already turned off.")
        
                     


            elif mutate_state == 1:
                print("New Gene Turned On")

                # Determine if Activator or Repressor is turned on.
                gene_turned_on = r.sample(["A", "R"], 1)[0]
                


                # Activator
                if gene_turned_on == "A":
                    if gene_network[4] == 1:
                        print("SILENT MUTATION")
                    else:
                        print("Activator gene is turned on")
                        for i in range(1, 5):                   # 1 - 4 are allocates to activators
                            if gene_network[i] == 0:
                                gene_network[i] = 1
                                link_network[i][0] = 1          # By default, linked to P
                                break
                    

                # Repressor
                elif gene_turned_on == "R":                 
                    if gene_network[9] == 1:
                        print("SILENT MUTATION")
                    else:
                        print("Repressor gene is turned on")
                        for i in range(5, len(gene_network)):            # 6 - 10 are allocated to repressors         
                            if gene_network[i] == 0:
                                gene_network[i] = 1
                                link_network[i][0] = -1          # By default, linked to P
                                break





        # LINKS
        elif mutation == 1:
            print("Mutation in Linkage")

            try:
                gene_link_mutate = r.sample(gene_index, 1)[0]   # determines which gene's link will be mutated
            
            except:
                print("Some problem with sampling from gene_index")

            try: 
                link_mutate = r.sample(gene_index, 1)[0]        # determines which link will be mutated amongst the active genes
            except:
                print("Some problem with sampling from gene_index")


            # Determine to turn on or off
            mutate_state = r.sample(mutate_states, 1)[0]  # Determines if turned on or off
            try:
                if mutate_state == 0:
                    try: 
                        link_network[gene_link_mutate][link_mutate] = 0                
                    except:
                        print("Problem")
                

                elif mutate_state == 1:
                    if gene_link_mutate < 5:                            # Activator (including P)
                        link_network[gene_link_mutate][link_mutate] = 1
                    elif gene_link_mutate >= 5:                           # Repressors
                        link_network[gene_link_mutate][link_mutate] = -1
            except:
                print("Error. Exception likely occured above.")

        # PARAMETERS
        elif mutation == 2:
            print("Mutation in Parameters")

        print("\ngene network post mutation:\n {}\n".format(gene_network))
        print("link network post mutation:\n {}\n".format(link_network))
        print("Leaving mutation...\n\n")


    




# =-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=- END GENETIC NETWORK CLASS -=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-==-=-=-= #


