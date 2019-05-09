function [tours,cost_path,k,numSubtours, cpu] = subConstraint_Dantzig(data,dist_matrix,ncities)

starttime = clock;
%% CALCULATE DISTANCE EUCLIDEAN MATRIX

% Then we obtain the maximum value in the distance matrix 
% This max value is then added to the zeros diagonal (+10) so that it is
% not considered as a minimum value
max_Value = max(max(dist_matrix));
dist_matrix(dist_matrix == 0) = max_Value + 10;

% The distance matrix is then reshaped as a vector so that it can be input
% in the integer linear program function
dist_matr_2D = reshape(dist_matrix, [1 ncities^2]);

% Finally, the matrix representing any possible path is computed
cityPairs = fliplr(fullfact([ncities ncities]));

% We initialise the counter of total number of subtours obtained
numSubtours = 0;
%% CONSTRUCTING EQUALITY CONSTRAINT MATRIX
% We call the function that structures the equalities constraint matrix
[Aeq] = constma(ncities);

% The RHS of the equalities constraint matrix is computed, equal to one
size_beq = size(Aeq,1);
beq = ones(size_beq,1);

%% INITIALISATION STEP

%SOLVING THE ASSIGNMENT PROBLEM RELAXATION
% We set inequality constraints (non-existent in the relaxation) to use the intlinprog
% function
A = [];
b = [];
% We set the lower bound as 0 and the upper bound as 1 for all possible
% solutions (binary integer linear problem)
ub = ones(ncities^2,1);
lb = zeros(ncities^2,1);

% We define the options of the optimization solver
opts = optimoptions('intlinprog','Display','off');

% Solution for the ASSIGNMENT PROBLEM RELAXATION
[sol_path, cost_path,exitflag,output] = intlinprog(dist_matr_2D,1:ncities^2,A,b,Aeq,beq,lb,ub,opts);
fprintf('Total distance of the optimal route: %d\n', cost_path);

% DETECTING SUBTOURS
% First the number of subtours is detected, then printed in the command
% window
tours = detectSubtours(sol_path,cityPairs);
numberOfTours = length(tours); 
fprintf('Number of subtours: %d\n',numberOfTours);

% We compute the number of subtours obtained 
numSubtours = numSubtours + numberOfTours;
%% ITERATIVE STEP

% CREATION OF INITIAL INEQUALITY CONSTRAINT MATRIX
% The function is called, with the tours created in the relaxation, number
% of subtours and number of cities as inputs
% The function provides the LHS and RHS matrices
[A_tour,b_tour] = dantzigConst(tours,numberOfTours,ncities);

% The general inequality contraints are generated;
A = [];
b = [];

% A while is used in order to continue with the iterative process until a
% single subtour is obtained
while numberOfTours > 1
    
    % A and b (LHS and RHS inequality constraint matrices) are updated for
    % each iteration, including the latest subtours obtained. The matrix
    % contained the latest subtours as well as all the previous subtours
    % obtained in previous iterations
    A = [A; A_tour];
    b = [b; b_tour];
    
    % The optimisation solver is computed, considering the subtours
    % obtained in each iteration
    [sol_path, cost_path,exitflag,output] = intlinprog(dist_matr_2D,1:ncities^2,A,b,Aeq,beq,lb,ub,opts);


    % The number of subtours obtained is computed for each iteration and
    % printed
    tours = detectSubtours(sol_path,cityPairs);
    numberOfTours = length(tours); 
    fprintf('Number of subtours: %d\n',numberOfTours);
    
    % The subtour inequality constraint matrix is updated with the latest
    % subtours obtained    
    [A_tour,b_tour] = dantzigConst(tours,numberOfTours,ncities);
    
    % The total number of subtours and the number of iterations is updated    
    numSubtours = numSubtours + numberOfTours;
    k = k + 1;
    
end

endtime = clock;
cpu = etime(endtime, starttime);

hold off

%% PRESENTATION OF RESULTS
fprintf('\nSolution Quality: %f (lesser the better)\n',output.absolutegap);
fprintf('Optimized tour route:');
celldisp(tours);
fprintf('Note: The numbers correspond to order of cities in the input file\n');
fprintf('Total distance of the optimal route: %d\n', cost_path);

end