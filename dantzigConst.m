function [A,b] = dantzigConst(tours,numtours,n)

% First we initialise the LHS and RHS matrices with rows = number of
% subtours obtained (input) and columns equal to number of cities ^ 2
A = zeros(numtours,n^2);
b = zeros(numtours,1);

% The first for loop will take each subtour individually
for i = 1:numtours
    % We extract each subtour and calculate the size of each subtour
    tourIdx = cell2mat(tours(i));
    sizeSubTour = size(tourIdx,2);
    
    % We set a second for loop that will consider every node in the tour as
    % the initial point, i.e. for xij every i is considered
    for j = 1:sizeSubTour
        % we create a dummy tour that will be updated
        tourIdx_left = tourIdx;
        
        % each city is assigned as i in xij
        init_city = tourIdx(:,j);
        % We eliminate init_city from the tour
        tourIdx_left(:,j) = [];
        
        % The third for loop considers all the remaining cities in the
        % subtour as j in xij
        
        for z = 1:sizeSubTour-1
            sec_city = tourIdx_left(:,z);
            % We assign one to each specific xij
            A(i,(1+((init_city-1)*n))+(sec_city-1)) = 1;
        end
    end
    % For each subtour, the RHS is set as the size of the subtour - 1
    b(i,1) = sizeSubTour - 1;
end    


end
