# Selection simulations
To run this code requires Julia 1.5+ and the packages LinearAlgebra, Distributions, DelimitedFiles, Distributed, and Random

# Coalescent times
Scripts to simulate coalescent times are 4treesim.jl (for 4 lineages) and 6treesim.jl (for 6 lineages)

Step 1: Open 4treesim.jl (to simulate coalescent times for 4 lineages) or 6treesim.jl (to simulate coalescent times for 4 lineages) with any text editor

Step 2: Change line 8 to cd("filepath") where filepath is where the output will be written as csv files

Step 3: On the last line set the parameters:
        r = selection avantage,
        p_m = mutation rate,
        N = population size,
        nreps = number of coalescent times to simulate
        
Step 4: Start julia from the terminal, and run the desired script with include("script.jl")

fourtreetimes.csv and sixtreetimes.csv contain the simulated coalescent times

fourtreeh.csv and sixtreeh.csv contain the simulated distribution of h

# Tree topologies
Scripts to gene trees given a 4-taxon speices tree are ILSsym.jl (for a symmetric species tree) and ILSasym.jl (for and asymmetric species tree)

Step 1: Open ILSsym.jl (for a symmetric species tree) or ILSasym.jl (for and asymmetric species tree) with any text editor

Step 2: Change line 8 to cd("filepath") where filepath is where the output will be written as csv files

Step 3: On the last line set the parameters:
        r = selection avantage,
        p_m = mutation rate,
        tau_1, tau_2, and tau_3 are the three absolute speciation times,
        nreps = number gene trees to simulate
        
Step 4: Start julia from the terminal and run the desired script with include("script.jl")

fourtreecoal_sym.csv and fourtreecoal_asym.csv contain the coalescent times and topologies of the simulate gene trees.
