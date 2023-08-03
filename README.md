# RandomDynamics

Welcome to the Random Dynamics package! This package was created to run simulations and trajectory tracking of Random Dynamical Systems, abbreviated as 'RDS' (or 'rds') throughout the code. Colby Shaw and Ethan Botelho created this code with the help of Dr. Jorge Gonzalez during the summer of 2023 as a research project conducted by Dr. Alex Blumenthal at the Georgia Institute of Technology.

# Functions

Here are a list of the functions that you will find in the 'rds.jl' file; please consult this file and the documentation within for specific usage.

## rds.jl
- sampleTraj (Three Required Arguments, Two Default Arguments)
    - Returns the trajectory of a given sample with respect to the given RDS and function.
- timeSeries (Two Required Arguments)
    - Returns a time series calculation for a given trajectory and observable.
- timeSeries (Three Required Arguments)
    - Returns a time series calculation for a given trajectory and random observable.
- empiricalAverage (One Required Argument)
    - Computes the empirical average of a trajectory.
- sampling (Two Required Arguments, One Default Argument)
    - Samples the indicated number of values from the given distribution.
- makeDistribution (Three Required Arguments)
    - Returns a function for estimating a change of domain from a given function to a domain with intervals modulo 1.
- makeDistributionCross (Two Required Arguments)
    - Returns a function combining given distributions with different weights.
- makeBigFloat (Two Required Arguments)
    - Creates a BigFloat approximation for a given Float64 number.

## graphics.jl
- testing (One Required Argument)
    - Plots the trajectory of a given trajectory.
- tracking (One Required Argument)
    - Creates a gif tracking the evolution of a given trajectory.
- distDisplay (Five Required Arguments)
    - Plots a surface plot for given domains and distributions up to the certain given approximation.

# Usage

For specific examples regarding these functions and defined structures, please see the 'packagePresentation.ipynb' file for some beginning examples to build off of.

# Issues

For any issues or bugs in the code found, feel free to contact us at:
- Ethan: ethanmbotelho@gmail.com
- Colby: