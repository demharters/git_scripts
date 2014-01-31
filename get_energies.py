#! /usr/bin/python3
import operator

class edit_files(object):

   
    def __init__(self,my_energy_file_input,my_pdb_file):

        self.my_energy_file = my_energy_file_input
        self.my_pdb_file = my_pdb_file


    def clean_file(self):

        self.my_data = {}
        i = 0

        with open(self.my_energy_file) as f:
            for line in f:
                try:
                    self.my_data[i] = float(line.rstrip())
                    i += 1

                except ValueError:
                    #print("Removed line: "+line)
                    pass

    def sort_values(self):

        self.my_sorted_keys = [] 
        self.my_sorted_values = []

        self.my_sorted_data = sorted(self.my_data.items(), key = lambda t: t[1])

        for x in self.my_sorted_data:
            self.my_sorted_keys.append(x[0])
            self.my_sorted_values.append(x[1])

    def filter_values(self,period):

        self.no_models = 0
        self.no_minima = 0

        #self.my_filtered_values = []
        
        self.my_temperature_minima = []
        self.my_temperature_minima.append(int(period) - int((period/4)))

        for value in range(len(self.my_data)):

            self.no_models += 1
        
        for i in range(1,int((self.no_models/period))):
            
            self.my_temperature_minima.append(self.my_temperature_minima[i-1] + period)
            self.no_minima += 1

    def write_models(self):
    # Takes full trajectory as input and returns trajectory filtered for
    # values in self.my_temperature_minima.
        
        my_pdb_output = open("minima.pdb","w")

        with open(self.my_pdb_file) as f:
            
            parsing = False
            i = 1

            for line in f:

                if "MODEL" in line and int(line.split()[1]) in self.my_temperature_minima:
                    
                    parsing = True

                    print(str(line.split()[0]),str(i),file= my_pdb_output)
                 
                    i += 1

                elif "MODEL" in line and int(line.split()[1]) not in self.my_temperature_minima:

                    parsing = False

                else:
                    pass



                if "MODEL" not in line and parsing:

                    print(line.rstrip("\n"),file= my_pdb_output)

                else:

                    pass


        my_pdb_output.close()

    def write_low_energy_models(self,no_low_energy_models):

        low_energy_pdb = open("low_energy.pdb","w")

        self.the_chosen_ones = self.my_sorted_keys[:no_low_energy_models]
        
        for model in self.the_chosen_ones:


            with open(self.my_pdb_file) as f:

                parsing = False
                i = 1

                for line in f:
    
                    if "MODEL" in line and int(line.split()[1]) == model:
    
                        parsing = True
   
                        print(str(line.split()[0]),str(i),file = low_energy_pdb)
    
                        i += 1
    
                    elif "MODEL" in line and int(line.split()[1]) != model:
   
                       parsing = False

                    else:
                        pass
    


                    if "MODEL" not in line and parsing:

                        print(line.rstrip("\n"),file = low_energy_pdb)

                    else:
                        pass
        

    def write_energy(self):
        
        my_energy_output = open("energy.dat","w")

        for i in self.my_temperature_minima:
            print(self.my_data[i],file = my_energy_output)

        my_energy_output.close()

    def print_values(self,no_of_values):
        
    # Print out keys and values of dictionary generated in sort_value().
    # The maximum length of the first column (i.e. the keys) is calculated
    # and used as a reference for the second column holding the values.
    # Function takes one argument: number of values to be printed
    
        self.col_width = max(len(str(row)) for row in self.my_sorted_keys) + 2

        for i in range(no_of_values):

            print(str(self.my_sorted_keys[i]).ljust(self.col_width),
                    str(self.my_sorted_values[i]).ljust(self.col_width))

energy_file = edit_files("sampled.pot_energy","sampled.pos.pdb")
my_energies = energy_file.clean_file()
energy_file.sort_values()
#energy_file.print_values(10)
energy_file.filter_values(4)
energy_file.write_models()
#energy_file.write_energy()
#energy_file.write_low_energy_models(20)
