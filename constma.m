% Function to create the constraint matrix
function [constraint_matrix] = constma(num)

% First we create the vector constraints for both the tour visiting each
% city once and only once and the tour returning to the starting city 


% Second we create the matrix for the supply points
sup_p = zeros(num,num);

%We initialise the matrix for node 1
i = 1;
for j = 1:num
   sup_p(i,j) = 1;
end

% We then proceed to run the same process for all the nodes
% We concatenate the matrix of each node to the main matrix
for i = 1:num-1
    A = zeros(num,num);
    A(i+1,:) = 1;
    sup_p = [sup_p A];
end

% Then we create the demand points

% The demand point matrix is the identity matrix * -1
dem_p = eye(num,num);
%dem_p = eye(num-1,num) * -1;

for i = 2:num
    A = eye(num,num);
    dem_p = [dem_p A];
end


% Finally we concatenate both matrices to obtain the constraint matrix

constraint_matrix = [sup_p;dem_p];
