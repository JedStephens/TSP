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
solutions_names = table({'Assignment Problem Relaxation';'Dantzig Subtour Constraints';'Connectivity Constraints';'MTZ Subtour Constraints';'Desrochers Subtour Constraints';'Single Commodity Flow-based Formulation'}, [1;2;3;4;5;6]);
solutions_names.Properties.VariableNames = {'Method','id'};
% Set the optimal solution provided by the TSPLib
optimal_solution = 426;

% To store the paths



%First we read the data from the txt file
data = readtable('EIL51.txt');

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

[path_conn,opt_sol_conn,iter_conn,numSubtours_conn,cpu] = subConstraint_Conn(data,dist_matrix,ncities);
solutions(i,1) = i;
solutions(i,2) = opt_sol_conn;
solutions(i,3) = numSubtours_conn;
solutions(i,4) = iter_conn;
solutions(i,5) = (opt_sol_conn/optimal_solution) - 1;
solutions(i,6) = cpu;

path_conn = cell2mat(path_conn);
paths(:,i-1) = path_conn;
i = i+1;

[path_MTZ,opt_sol_MTZ,iter_MTZ,numSubtours_MTZ,cpu] = subConstraint_MTZ(data,dist_matrix,ncities);
solutions(i,1) = i;
solutions(i,2) = opt_sol_MTZ;
solutions(i,3) = numSubtours_MTZ;
solutions(i,4) = iter_MTZ;
solutions(i,5) = (opt_sol_MTZ/optimal_solution) - 1;
solutions(i,6) = cpu;

path_MTZ = cell2mat(path_MTZ);
paths(:,i-1) = path_MTZ;
i = i+1;

[path_Desr,opt_sol_Desr,iter_Desr,numSubtours_Desr,cpu] = subConstraint_Desr(data,dist_matrix,ncities);
solutions(i,1) = i;
solutions(i,2) = opt_sol_Desr;
solutions(i,3) = numSubtours_Desr;
solutions(i,4) = iter_Desr;
solutions(i,5) = (opt_sol_Desr/optimal_solution) - 1;
solutions(i,6) = cpu;

path_Desr = cell2mat(path_Desr);
paths(:,i-1) = path_Desr;
i = i+1;

[path_flow,opt_sol_flow,iter_flow,numSubtours_flow,cpu] = subConstraint_flow(data,dist_matrix,ncities);
solutions(i,1) = i;
solutions(i,2) = opt_sol_flow;
solutions(i,3) = numSubtours_flow;
solutions(i,4) = iter_flow;
solutions(i,5) = (opt_sol_flow/optimal_solution) - 1;
solutions(i,6) = cpu;

path_flow = cell2mat(path_flow);
paths(:,i-1) = path_flow;
i = i+1;


%% TABLES
solutions = array2table(solutions);

solutions.Properties.VariableNames = {'id', 'Solution_obtained', 'Number_of_Subtours', 'Number_of_iterations', 'Deviation_from_optimal', 'CPU_time'};
solutions_table = join(solutions_names, solutions);
solutions_table.id = [];

writetable(solutions_table,'solutions.xls');

paths = array2table(paths);

paths.Properties.VariableNames = {'Dantzig_Subtour_Constraints', 'Connectivity_Constraints', 'MTZ_Subtour_Constraints', ...
                                     'Desrochers_Subtour_Constraints', 'Single_Commodity_Flow_based_Formulation'};

                                 
writetable(paths,'paths.xls');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% function [dist_best, cputime] = SA_algorithm(dist_matrix,connection, TempMin, B, num_neigh_to_visit)
% 
% Author1: Group 4
% Date:    02/11/17
% Course:  Business Analytics with Heuristics
% 
% Function   : SA_algorithm 
% 
% Description: The SA algorithm runs the simulated annealing method of
%              for the solutions (paths) obtained by running the 2opt or
%              3opt methods to evaluate new solutions
%
% Parameters : dist_matrix        - the Euclidean Distance Matrix of each 
%                                   instance 
%              connection         - option for running 2opt or 3opt
%              TempMin            - minimum temperature
%              B                  - cooling constant
%              num_neigh_to_visit - number of neighbours to visit
% 
% Return     : dist_best     - the new best distance
%              cputime       - running time 
%


function [dist_best, cputime] = SA_algorithm(dist_matrix,connection, TempMin, B, num_neigh_to_visit)
    tic
    
    % First step is to run the nearest neighbour algorithm to obtain the
    % initial best solution
    [path_distance, ~, path] = NN_algorithm(dist_matrix);

    % We create two dummy variables, path_better and dist_better, that will
    % be updated if the solution is accepted
    path_better = path;
    dist_better = path_distance;

    % We create two additional dummy variables that will store the best
    % solution recorded during the algorithm run, which will be the one
    % that is considered as output of this function
    path_best = path_better;
    dist_best = dist_better;
    
    % The initial temperature is calculated based on the number of moves
    % accepted
    T = initialT(path_better, dist_better, dist_matrix, connection);

    t = 1; % temperature change counter

    % We set the number of iterations to infinity so it does not impact the
    % solutions obtained
    num_iterations = inf;

    BestSolution = [];
    j = 1;

    % The first for loop runs for all iterations or until T is lower than
    % the minimum temperature
    for n = 1:num_iterations
        k = 0;
        % This for loop runs for the number of neighbours to visit,
        % providing a best solution at the end of each loop
        for i = 1:num_neigh_to_visit
            % This if statement will run either 2opt or 3opt based on the
            % input of the function
            % It will return a path to be considered as new better solution
            if strcmp(connection, '2opt')
                [path_new, dist_new] = RD_twoopt(path_better,dist_matrix);
            elseif strcmp(connection, '3opt')
                [path_new, dist_new] = RD_threeopt(path_better,dist_matrix);   
            end
            path_new = path_new';
            % The difference in cost is computed
            d = dist_new - dist_better;
            % The difference in cost and the current T are the inputs for
            % the probability function
            P = probabifunction(d,T);
            % The value from the probability function is then compared
            % against a random number between 0 and 1
            if P >= rand(1)
                % If the condition is true, then path_better is updated 
                path_better = path_new;
                dist_better = dist_new;
                BestSolution(1,j) = dist_better;
                j = j+1;
                % If path_better is better than path_best, then path_best
                % is updated
                if dist_best > dist_better
                    dist_best = dist_better;
                    path_best = path_better;
                end
            end
            k = k + 1;
        end
        % The temperature is updated, considering the current T and the
        % cooling constant
        t = t+1;
        T = B * T;
        
        % Finally, the current temperature is computed against the minimum
        % temperature
        % If it is lower, the algorithm stops
        if T <= TempMin
            break
        end
    end
    cputime = toc;

 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% function [T] = initialT(path, path_cost, dist_matrix, opt)
% 
% Author1: Group 4
% Date:    02/11/17
% Course:  Business Analytics with Heuristics
% 
% Function   : initialT
% 
% Description: This function returns the initial temperature to be
%              considered for the Simulated Annealing method, based on the
%              fact that the initial temperature is set so that 80% or more
%              of the moves generated by the 2opt or 3opt function are
%              accepted
%
% Parameters : path               - best path obtained running the nearest
%                                   neighbour method
%              path_cost          - cost of path
%              dist_matrix        - the Euclidean Distance Matrix of each 
%                                   instance 
%              opt                - 2opt or 3opt
% 
% Return     : T                  - Initial Temperature



function [T] = initialT(path, path_cost, dist_matrix, opt)

% The initial temperature is initialised at 100 degrees
T0 = 100;

while true
    
    % We gather the number of moves that are accepted
    Accepted = 0;
    NonAcc = 0;
    % 2opt or 3opt are run for 100 iterations
    for i = 1:100
        if strcmp(opt, '2opt')
            [path_new, dist_new] = RD_twoopt(path,dist_matrix);
        elseif strcmp(opt, '3opt')
            [path_new, dist_new] = RD_threeopt(path,dist_matrix);  
        end
        path_new = path_new';
        d = dist_new - path_cost;
        P = probabifunction(d,T0);
        % The number of moves is gathered (i.e. p is greater than the
        % random number generated)
        % Otherwise, the number of non-accepted solutions is increased
        if P >= rand(1)
            path_cost = dist_new;
            path = path_new;
            Accepted = Accepted + 1;
        else
            NonAcc = NonAcc + 1;
        end      
    end
    
    % Finally, if the number of accepted moves is less than 80%, the
    % temperature is doubled, otherwise the function stops and returns the
    % initial temperature
    if round((Accepted/100),1) < 0.8
        T0 = T0*2;
    else
        T = T0;
        return
        
    end
end       
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function [shortestPathLength, ctime, shortestPath] = NN_algorithm(distances)
% 
% Author1: Group 4
% Date:    06/10/17
% Course:  Business Analytics with Heuristics
% 
% Function   : NN_algorithm 
% 
% Description: The NN algorithm runs the nearest nieghbour algorithm,
%              providing a best path and the distance associated
%              to the path
%
% Parameters : distances          - the Euclidean Distance Matrix of each 
%                                   instance 
% 
% Return     : ShortestPathLength     - the new best distance
%              ctime                  - running time 
%              shortestPath           - the new shortest path
%



function [shortestPathLength, ctime, shortestPath] = NN_algorithm(distances)
% This function runs the nearest neighbour algorithm for the file selected
% in the main script

% Time at the start
tic

% Assign data from each file to the cities variables, extracting only node
% coordinates
% Get the size, i.e. number of nodes
N_cities = size(distances,1);

% We turn the 0 in the diagonal of the distance matrix to the largest
% number possible to avoid picking 0 as the minimum value
%distances(distances==0) = realmax;

for i = 1:N_cities
    distances(i,i) = realmax;
end

% We also set the initial shortest path as the largest number as well to
% make sure that the first iteration will update this value
shortestPathLength = realmax;

% This for loop is used to assign every node in the dataset as the starting
% city in each iteration
for i = 1:N_cities
    
    % Assigns starting city to the node i
    startCity = i;

    % Adds the starting city to path, which will store the city in visiting
    % order
    path = startCity;
    
    % Sets distanceTraveled to 0 and distancesNew as a copy of the distance
    % matrix that will be updated for each node added to path
    distanceTravelled = 0;
    distancesNew = distances;
    
    % Finally, it assigns the starting city as the current city for the
    % nested for loop
    currentCity = startCity; 
    
    for j = 1:N_cities-1
        
        % This function gets the minimum distance to the current node as
        % well as the index of the next city, taken from the row of the
        % current city in the distance matrix
        [minDist,nextCity] = min(distancesNew(:,currentCity));
        % Tie breaker
        if (length(nextCity) > 1)
            nextCity = nextCity(1);
        end
        
        % Path is updated with the next city 
        path(end+1,1) = nextCity;
        % distance travelled is updated with the distance to the nearest
        % node
        distanceTravelled = distanceTravelled +...
                    distances(currentCity,nextCity);
        
        % The distance matrix is updated to make the row of the current
        % city equal to realmax so it is not considered as a minimum in the
        % following iterations
        distancesNew(currentCity,:) = realmax;
        
        % Assigns next city as the current city for the next iteration
        currentCity = nextCity;
        
    end
    
    % Adds the starting city at the end of the path vector and the distance
    % from the last node to the starting city is also computed to distance
    % travelled
    path(end+1,1) = startCity;
    distanceTravelled = distanceTravelled +...
        distances(currentCity,startCity);
    
    % This if statement is used to update shortestPathLength, which
    % represents the shortest distance for the overall n iterations, being
    % n the number of nodes
    if (distanceTravelled < shortestPathLength)
        shortestPathLength = distanceTravelled;
        shortestPath = path;
        
    end 
end

% Time at the end of the run
ctime = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% function prob = probabifunction(d,T)
% 
% Author1: Group 4
% Date:    02/11/17
% Course:  Business Analytics with Heuristics
% 
% Function   : probabifunction
% 
% Description: This function returns a value to be considered in order to
%              assess whether a solution can be accepted. This probability
%              function returns different values based on the temperature
%              and the difference of the costs
%
% Parameters : d       - the difference in cost between path_better and
%                        path_new
%              T       - current temperature

% Return     : prob    - probability of being accepted
%              
%

function prob = probabifunction(d,T)

% This function will return 1 if the cost of path_new is lower than the
% cost of path_better

if d <= 0
    prob = 1;
% Otherwise, the function will return a probability based on the difference
% in cost and the current temperature
else
    prob = exp(-d/T);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% function [new_path,best_distance] = RD_threeopt(path,dist_matrix,a,b,c)
% 
% Author1: Group 4
% Date:    02/11/17
% Course:  Business Analytics with Heuristics
% 
% Function   : RD_threeopt 
% 
% Description: The RD_threeopt procedure returns the new best path and length
%              so far. It splits the path in 3 edges randomly; or four  
%              seperated parts and reconnects them in all possible ways.
%              
% 
% Parameters : path          - initially set it as NN_path   
%              dist_matrix   - the Euclidean Distance Matrix of each 
%                              instance 
%              a,b,c         - the random points where the path is going 
%                              to be broken
% 
% Return     : new_path      - the new best path
%              best_distance - the new best distance 
%
function [new_path,best_distance] = RD_threeopt(path,dist_matrix)

% The four seperated parts after breaking the three edges
n = size(path,1)-1;
a = 0; b = 0; c = 0;

while a == b || b == c || a == c

a = randi([1 n]);
b = randi([1 n]);
c = randi([1 n]);

rand_number = [a b c];
rand_number = sort(rand_number);
end

path = path';

a = rand_number(1);
b = rand_number(2);
c = rand_number(3);


A = path(1,1:a);
B = path(1,a+1:b);
C = path(1,b+1:c);
D = path(1,c+1:end);

% The possible reconnection combinations
%comb1 = [A B C D];
comb2 = [A fliplr(B) C D];
comb3 = [A fliplr(B) fliplr(C) D];
comb4 = [A B fliplr(C) D];
comb5 = [A C B D];
comb6 = [A fliplr(C) B D];
comb7 = [A fliplr(C) fliplr(B) D];
comb8 = [A C fliplr(B) D];

% Store the combinations

combinations = [ comb2; comb3; comb4; comb5; comb6; comb7; comb8];
% Store the combinations distances in a matrix
all_combinations_distances = [];

% For the 7 possible combinations the distances are calculated
for  j = 1:7
    x = combinations(j,:);
    total_distance = 0;
    for i = 1:length(path)-1
        total_distance = total_distance+dist_matrix(x(i),x(i+1));       
    end
    all_combinations_distances(j,:) = total_distance;
end

% The minimum distance and the corresponding path are computed
[best_distance, idx] = min(all_combinations_distances);
new_path = combinations(idx,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% function [new_path, best_distance] = RD_twoopt(path,dist_matrix,a,b)
% 
% Author1: Group 4
% Date:    02/11/17
% Course:  Business Analytics with Heuristics
% 
% Function   : RD_twoopt 
% 
% Description: The RD_twoopt procedure returns the new best path and length
%              so far. It splits the path in 2 edges randomly; or three  
%              seperated parts and reconnects them in all possible ways.
%              
% 
% Parameters : path          - initially set it as NN_path   
%              dist_matrix   - the Euclidean Distance Matrix of each 
%                              instance 
%              a,b           - the random points where the path is going to 
%                              be broken
% 
% Return     : new_path      - the new best path
%              best_distance - the new best distance 
%
function [new_path, total_distance] = RD_twoopt(path,dist_matrix)

n = size(path,1)-1;
a = 0;
b = 0;

path = path';
while a == b
    a = randi([1 n]);
    b = randi([1 n]);
end
% The three seperated paths after breaking the initial path
if a < b
    A = path(1,1:a);
    B = path(1,a+1:b);
    C = path(1,b+1:end);
elseif b < a
    A = path(1,1:b);
    B = path(1,b+1:a);
    C = path(1,a+1:end);
end
% The possible reconnection combinations
%comb1 = [A B C];

new_path = [A fliplr(B) C];

% Store the combinations
%combinations = [comb1; comb2;];

% Store the combinations distances in a matrix
%all_combinations_distances = [];

% For the 2 possible combinations the lengths are calculated
%for  j = 1:2
%    x = combinations(j,:);
    total_distance = 0;
    for i = 1:length(path)-1
        total_distance = total_distance+dist_matrix(new_path(i),new_path(i+1));       
    end
%    all_combinations_distances(j,:) = total_distance;
%end
% The minimum distance and the corresponding path are computed
%[best_distance, idx] = min(all_combinations_distances);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

