function [A,b] = ConnectConst(tours,numtours,n)

    % First we initialise the LHS and RHS matrices with rows = number of
    % subtours obtained (input) and columns equal to number of cities ^ 2
    A = zeros(numtours,n^2);
    b = zeros(numtours,1);
    
    % A dummy variable containing the tours is created
    tours_iter = tours;

    % The first for loop will take each subtour individually
    for i = 1:numtours
        % We extract each subtour and calculate the size of each subtour
        tour_1 = cell2mat(tours(i));        
        size_tour1 = size(tour_1,2);
        
        % We eliminate the tour considered from the dummy variable 
        tours_iter(:,i) = [];
        
        % The second for loop considers all the points in the tour as the
        % initial point, i.e. i in xij        
        for j = 1:size_tour1
            % Each city in tour_1 is selected as the initial city
            tour_init = tour_1(:,j);
            
            % The third for loop goes through the remaining subtour
            for jj = 1:numtours-1
                % Each remaining subtour is selected and name as tour_2
                tour_2 = cell2mat(tours_iter(jj));
                % The size of tour_2 is computed
                size_tour2 = size(tour_2,2);
                
                % The fourth for loop goes through all the nodes in tour_2
                % and assign it as the final city, i.e. j in xij
                for z = 1:size_tour2
                    % Each node in tour_2 is considered
                    tour_end = tour_2(:,z);
                    % xij is set equal to -1
                    A(i,(1+((tour_init-1)*n))+(tour_end-1)) = -1;

                end
            end
        end
        
        % The RHS value is set -1
        b(i,:) = -1;
        
        % The dummy variable is reset with the original subtours obtained
        tours_iter = tours;
    end    


end
