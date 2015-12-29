import numpy as np
import matplotlib.pyplot as plt
import random as r
from scipy.integrate import odeint
import copy

class ODE_System():
    """ """

    def __init__(self, network, network_parameters, morphogen_input, cell_position):
        """ """

        self.gene_network = network[0]
        self.link_network = network[1]
        
        self.gene_network_parameters = network_parameters[0]
        self.link_network_parameters = network_parameters[1]

        self.G = morphogen_input(cell_position)

        self.initial_conditions = np.zeros(10)
        print("\n\t\tInitializing cell {}".format(cell_position))
        print("\t\tG_input: {}".format(self.G))



    def get_metabolism(self):
        return self.metabolism



    def morphogen_factor(self, G, n, K):
        """ """
        G_factor = G**n / (G**n + K**n)

        return G_factor


    def activator_factor(self, n, K, Ai):
        """ """
        A_factor = Ai**n / (Ai**n + K**n)

        return A_factor


    def repressor_factor(self, n, K, Ri):
        """ """
        R_factor = 1 / (1 + Ri**n / K**n)
        return R_factor



    def P_growth(self, i):

        S = self.gene_network_parameters[i]["S"]
        n = self.gene_network_parameters[i]["n"]
        K = self.gene_network_parameters[i]["K"]
        G_factor = self.morphogen_factor(self.G, n, K)
    
        P = S * G_factor

        return P



    def P_link(self, i, j):
        """ """
        
        n = self.link_network_parameters[i][j]["n"]
        K = self.link_network_parameters[i][j]["K"]
        P_factor = self.activator_factor
        
        P_link = [P_factor, n, K]

        return P_link


    def A_growth(self, i):
        S = self.gene_network_parameters[i]["S"]
        n = self.gene_network_parameters[i]["n"]
        K = self.gene_network_parameters[i]["K"]
        G_factor = self.morphogen_factor(self.G, n, K)

        parameters = [S, G_factor]
    
        A = S * G_factor

        return A



    def A_link(self, i, j):
        """ """
        n = self.link_network_parameters[i][j]["n"]
        K = self.link_network_parameters[i][j]["K"]
        A_factor = self.activator_factor
        
        A_link = [A_factor, n, K]

        return A_link


    def R_growth(self, i):
        S = self.gene_network_parameters[i]["S"]
        n = self.gene_network_parameters[i]["n"]
        K = self.gene_network_parameters[i]["K"]
        G_factor = self.morphogen_factor(self.G, n, K)

        parameters = [S, G_factor]
    
        R = S * G_factor

        return R



    def R_link(self, i, j):
        """ """
        n = self.link_network_parameters[i][j]["n"]
        K = self.link_network_parameters[i][j]["K"]
        R_factor = self.repressor_factor
        
        R_link = [R_factor, n, K]
        return R_link



    def link_factor(self, i, j):
        n = self.link_network_parameters[i][j]["n"]
        K = self.link_network_parameters[i][j]["K"]

        return [n, K]


    def decay_factor(self, decay, Xi):
        """ """

        decay_Xi = decay * Xi

        return decay_Xi


    def construct_system(self):
        print("\n\n")
        print("In Construct SYSTEM method...\n")
        print("Gene network: \n {}\n".format(self.gene_network))
        print("Link network: \n {}\n".format(self.link_network))

        self.ode_system_eqs = {}

      
        """ Determine Genes turned on """
        print("Determining Genes that are turned on...")
        for i in range(0, len(self.gene_network)):
            if self.gene_network[i] == 1:
                if i == 0:
                    # Take care of P
                    print("Adding P growth factor to ODE System...")
                    P_gene = self.P_growth(i)
                    self.ode_system_eqs[i] = {"growth": {"base": P_gene }}
                    print("Ode system:\n {}\n".format(self.ode_system_eqs))

                elif 0 < i < 5:
                    # Take care of Activators
                    print("Adding Activator growth factor to ODE System...")
                    A_gene = self.A_growth(i)
                    self.ode_system_eqs[i] = {"growth": {"base": A_gene }}
                    print("A_gene: {}".format(A_gene))
                    print("Ode system:\n {}\n".format(self.ode_system_eqs))

                elif 5 <= i < len(self.gene_network):
                    # Take care or Repressors
                    print("Adding Repressor growth factor to ODE System...")
                    R_gene = self.R_growth(i)
                    print("R_gene: {}".format(R_gene))
                    self.ode_system_eqs[i] = {"growth": {"base": R_gene }}
                    print("Ode system:\n {}\n".format(self.ode_system_eqs))
        



        genes_on = self.ode_system_eqs.keys()
        genes_on.sort()






        """ Determine Active Links """
        print("\nDetermining Active Links...")
        for gene in genes_on:

            for j in range(0, len(self.gene_network)):
                



                # ACTIVATORS

                if self.link_network[gene][j] == 1:        # P or Activator
                    print("Adding Activator!")


                    if gene == 0:      # Case of P
                        print("Currently not enabled")


                    elif gene > 0:           # Case of Activator
                        # Update Activator Linkages
                        
                        print("Updating genes that A affects...")
                        

                        print("ode system: {}".format(self.ode_system_eqs))
                        try:
                            growth_dict = self.ode_system_eqs.pop(j)
                            print("growth_dict: {}".format(growth_dict))
                            growth_values = growth_dict.values()[0]
                            print("growth_values: {}".format(growth_values))


                            new_values = self.link_factor(gene, j)
                            print("new_values: {}".format(new_values))

                            if "A" in growth_values.keys():
                                print("*****")
                                print("\n\n\n\nA already in value keys! Let's update it!!!\n")
                                # ode_system_eqs is library with gene number as key.
                                # self.ode_system_eqs.pop(j) already removed gene j and stored the growth library in growth_dict
                                # growth_values then stores the values of growth_dict, which is a new dict with the growth factors as keys.
                                # Here, we assume "A" is a key in growth_values.
                                A_growth_values = growth_values["A"]
                                A_growth_values.append(new_values)
                                # We now have new values for key "A" of growth_values.
                                growth_values["A"] = A_growth_values
                                # Since .pop(j) removed the growth library from the ode system, we now need to add
                                new_growth_values = {"growth": growth_values} 
                                self.ode_system_eqs[j] = new_growth_values


                                
                            else:
                                print("Activator link doesn't currently exist. Creating new link...")
                                growth_values["A"] = [new_values]

                                new_growth_values = {"growth": growth_values}
                                self.ode_system_eqs[j] = new_growth_values
                                print()

                            print("Ode system:\n {}\n".format(self.ode_system_eqs))

                        except:
                            print("ERROR")







                # REPRESSORS

                elif self.link_network[gene][j] == -1:    
                    print("Adding Repressor!")
                    # Update Repressor linkages
                    #try:
                    try:
                        growth_dict = self.ode_system_eqs.pop(j)
                        print("growth_dict: {}".format(growth_dict))
                        growth_values = growth_dict.values()[0]
                        print("growth_values: {}".format(growth_values))


                        new_values = self.link_factor(gene, j)
                        print("new_values: {}".format(new_values))

                        if "R" in growth_values.keys():
                            print("*****")
                            print("\n\n\n\nR already in value keys! Let's update it!!!\n")
                            # ode_system_eqs is library with gene number as key.
                            # self.ode_system_eqs.pop(j) already removed gene j and stored the growth library in growth_dict
                            # growth_values then stores the values of growth_dict, which is a new dict with the growth factors as keys.
                            # Here, we assume "A" is a key in growth_values.
                            R_growth_values = growth_values["R"]
                            R_growth_values.append(new_values)
                            # We now have new values for key "A" of growth_values.
                            growth_values["R"] = R_growth_values
                            # Since .pop(j) removed the growth library from the ode system, we now need to add
                            new_growth_values = {"growth": growth_values} 
                            self.ode_system_eqs[j] = new_growth_values
                            
                        else:
                            growth_values["R"] = [new_values]

                            new_growth_values = {"growth": growth_values}
                            self.ode_system_eqs[j] = new_growth_values
                            print()

                        print("Ode system:\n {}\n".format(self.ode_system_eqs))

                    except:
                        print("ERROR")

        """ Add decay constant to equations of ode_system """

        print("\nDetermining if Gene i has any regulators!")
        for gene in genes_on:
          
            try:
                print("\nAdding decay factor to ODE for gene {}".format(gene))
                values = self.ode_system_eqs.pop(gene)
                print("Values: {}".format(values))
                print("growth_values: {}".format(values))
            
                #values["decay"] = [self.decay_factor, self.gene_network_parameters[gene]["delta"]]
                values["decay"] = [self.gene_network_parameters[gene]["delta"]]
                print("updated_values: {}".format(values))

                self.ode_system_eqs[gene] = values
                print("Ode system:\n {}\n".format(self.ode_system_eqs))
            
            except:
                print("ERROR")  


            print("\n")


        print("\nExiting System Construction...\n")



    def metabolize(self):
    
        ode_system = self.ode_system_eqs


        print("\n\n\n")
        print("-"*20)
        print("Cell is Metabolizing!")
        print("-"*20)
        print("\n")
        print("ode system:\n")
        print(ode_system)
        print("\n")



        def repressor(values, Ri):
            
            R = 1

            idx = 0
            for value in values:

                n = value[0]
                K = value[1]
                R_x = 1 / (1 + Ri[idx]**n / K**n)
                R *= R_x

                idx += 1

            return R



        def activator(values, Ai):         

            A = []
            
            idx = 0
            for value in values:
                n = value[0]
                K = value[1]
                A_x = Ai[idx]**n / (Ai[idx]**n + K**n)
                A.append(A_x)
                idx += 1


            return max(A)
        

        # Determine which genes are turned on
        genes_on = self.ode_system_eqs.keys()
        genes_on.sort()

        if genes_on:

            def F(y, t):

                # the model equations
                solutions = []


                i = 0
                for gene in genes_on:
                 

                    # Growth
                    growth_factors = ode_system[gene]["growth"]
                    
                    growth_factor_keys = growth_factors.keys()
                    


                    # Base
                    growth_factor = ode_system[gene]["growth"]["base"]
                    
                    

                    # A  
                    if "A" in growth_factor_keys:
                        A_values = ode_system[gene]["growth"]["A"]


                        number_of_A = len(A_values)

                        A_idx = []
                        for kdx in range(1, 1 + number_of_A):       # From 1 because A starts at idx 1
                            A_idx.append(kdx)

                        y_A = []
                        for jdx in A_idx:
                            index = genes_on.index(jdx)
                            y_A.append(y[index])


                        A_factor = activator(A_values, y_A)
                    
                        #print("A_factor: {}".format(A_factor))

                        growth_factor *= A_factor
                        
                    
                    



                    # R                 
                    if "R" in growth_factor_keys:
                        # Need to find out how many repressors

                        R_values = ode_system[gene]["growth"]["R"]


                        number_of_R = len(R_values)

                        R_idx = []
                        for indx in range(5, 5 + number_of_R):              # From 5 because R starts at idx 5
                            R_idx.append(indx)

                        y_R = []
                        for idx in R_idx:
                            index = genes_on.index(idx)
                            y_R.append(y[index])

                        R_factor = repressor(R_values, y_R)
                    
                        #print("R_factor: {}".format(R_factor))

                        growth_factor *= R_factor
                        
                    

                    # Decay
                    decay_constant = ode_system[gene]["decay"][0]
                    decay = decay_constant * y[i]
                    


                    solution = growth_factor - decay

                    solutions.append(solution)
                    i += 1
               

                return solutions


            # Initial Conditions
            initial_conditions = []
            for gene in genes_on:
                initial_conditions.append(0.)




            # Time Grid
            t = np.linspace(0, 5., 100) 


            # solve the DEs
            soln = odeint(F, initial_conditions, t)
            
            P = soln[:, 0]
            P_max = soln[-1, 0]
            
            self.metabolism = P_max

            print("METABOLISM: {}".format(P_max))



        else:           # All genes are turned off.
            self.metabolism = 0


        print("\n\n******Leaving metabolism...\n\n")




