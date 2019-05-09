# The Travelling Salesman Problem (TSP)

The Travelling Salesman Problem is a classical problem used to illustrate the benefits of implementing mathematical programming algorithms to solve transportation routing problems. In particular, this case is called the **Assignment Problem**. 

The Assignment problem is a particular case of a transportation problem that considers the number of origins to be equal to the number of destinations (*m = n*), as well as the fact that each origin has a supply of 1 unit and each destination has a demand of 1 unit. 

When solving the assignment problem, the main objective is to optimise the number of resources for a number of activities so that the cost is minimised. 

In this case, two approaches are compared:
1. The Assignment Problem Relaxation
2. Dantzig, Fulkerson and Johnson Elimination Constraints (DFJ)

Whereas the Assignment Problem Relaxation allows for subtours to be created, the DFJ algorithm constraints the creation of subtours, building a full solution to the problem.

## TO DO
- [x] Optimise, clean and refactor Matlab Code
- [ ] Add documentation
- [ ] Translate + refactor + CLI development in Python for user integration
