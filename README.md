# coalescent_under_selection
Code to simulate coalescent times from the coalescent under selection

To run this code requires Julia 1.5+ and the packages LinearAlgebra, Distributions, DelimitedFiles, Distributed, and Random

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

