# Genetic Algorithm code base for python

Amlan Talukder, Chathura Jayalath

Date: April, 2020

Choosing the genetic operator types and tuning operator probabilities is an additional pre-processing step needed before running a genetic algorithm, especially since the choice varies with different problems and have significant effects on the output. To help with this issue, adaptive operator strategies were already proposed by various studies. Here we propose an extension on the adaptive operator strategy proposed by a previous study which adapts the probability of the selected multiple operators for a problem based on how the operator contributed to the fitness values so far. We also report a better result than the original study using our extension.

Branches
---------------------------------------------------------------
1. GA: Vanilla GA
2. SSGA: Steady state GA
3. Adaptive GA: SSGA proposed by Julstrom, 1995
4. Adaptive GA ext: Extension proposed by us


1. run this command 
    "python codes/search.py <problem_name>.params"
2. use "1" as additional argument to print more details
    "python codes/search.py <problem_name>.params 1"

References
---------------------------------------------------------------
Bryant A Julstrom. "What have you done forme lately? Adapting operator probabilities in a steady-state genetic algorithm". In:Proceedings of the 6th International Conference on Genetic Algorithms. 1995,pp. 81–87.
