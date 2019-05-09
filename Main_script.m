% Main Script
%
% Author:  Group 4
% Date:    13/11/17
% Course:  Prescriptive Analytics with Mathematical Programming
% 
% 
% Description: 
%             
%% READ TABLE AND OPTIMAL SOLUTION
% Starts the solution matrix
solutions =[];
paths = [];
i = 1;
solutions_names = table({'Assignment Problem Relaxation';'Dantzig Subtour Constraints'}, [1;2]);
solutions_names.Properties.VariableNames = {'Method','id'};
% Set the optimal solution provided by the TSPLib
optimal_solution = 426;

%First we read the data from the txt file
data = readtable('data/EIL51.txt');

%% CALCULATE DISTANCE EUCLIDEAN MATRIX
% First we calculate the number of cities to connect
ncities = size(data,1);

% We convert the data to matrix format and then select the coordinate
% columns
coord_mat = table2array(data);
coord_mat = coord_mat(:,2:3);

% We compute the distance matrix based on the euclidean distance
dist_matrix = squareform(pdist(coord_mat,'Euclidean'));

% We solve the problem using the Assignment problem relaxation
[opt_sol_relax,numSubtours, cpu] = assignment_relaxation(data,dist_matrix,ncities);
solutions(i,1) = i;
solutions(i,2) = opt_sol_relax;
solutions(i,3) = numSubtours;
solutions(i,4) = 1;
solutions(i,5) = (opt_sol_relax/optimal_solution) - 1;
solutions(i,6) = cpu;
i = i + 1;

[path_dantzig,opt_sol_dantzig,iter_dant,numSubtours_dant,cpu] = subConstraint_Dantzig(data,dist_matrix,ncities);
solutions(i,1) = i;
solutions(i,2) = opt_sol_dantzig;
solutions(i,3) = numSubtours_dant;
solutions(i,4) = iter_dant;
solutions(i,5) = (opt_sol_dantzig/optimal_solution) - 1;
solutions(i,6) = cpu;

path_dantzig = cell2mat(path_dantzig);
paths(:,i-1) = path_dantzig;
i = i+1;
%% TABLES
solutions = array2table(solutions);

solutions.Properties.VariableNames = {'id', 'Solution_obtained', 'Number_of_Subtours', 'Number_of_iterations', 'Deviation_from_optimal', 'CPU_time'};
solutions_table = join(solutions_names, solutions);
solutions_table.id = [];

writetable(solutions_table,'results/solutions.xls');

paths = array2table(paths);

paths.Properties.VariableNames = {'Dantzig_Subtour_Constraints'};

                                 
writetable(paths,'results/paths.xls');

