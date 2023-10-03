function joinedData = mergeSpice(point, elements)
    %point structure
    
    
    %elements structure
    % 1. resistor number | 2. resistance | 3. current | 4. power 
    
    % removing non-percolating wires
    point = point( (point(:,5) == -1 & point(:,6) == -1), :);
    elements = sortrows(elements);
    
    
    joinedData = [point(:,1:6), elements(:,2:4)];
    %output structure is
    % 1. node1 | 2. node2 | 3. x_coord | 4. y_coord | 5. top_connection? |
    % 6. bottom_connection? | 7. resistance | 8. current | 9. power
    
    
end