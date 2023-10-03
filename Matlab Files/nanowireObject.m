classdef nanowireObject
    % Class of nanowires with tunable parameters
    properties
        %geometrical properties%
        x        % x coordinate of center
        y        % y coordinate of center
        x1       % x and y endpoints
        x2
        y1
        y2
        length      % length of NW
        theta       % polar angle, measured from z axis
        aRow    % this is the A row for this line in the Ax = b matrix equation
        bRow
        xMin
        xMax
        yMin
        yMax
        limsMat
        
        %percolation properties
        wireNum
        connectedWires
    end
    
    methods
        function nw = nanowireObject(lims,length) %class constructor
            %base properties
            nw.length = length;
            nw.theta = 2 * pi * rand;
            nw.x = lims(1,1) + ( lims(1,2) - lims(1,1) ) * rand;
            nw.y = lims(2,1) + ( lims(2,2) - lims(2,1) ) * rand;
            
            %calculated properties
            nw.x1 = nw.x + length/2 * cos(nw.theta);
            nw.x2 = nw.x - length/2 * cos(nw.theta);
            nw.xMin = min(nw.x1,nw.x2);
            nw.xMax = max(nw.x1,nw.x2);
            
            nw.y1 = nw.y + length/2 * sin(nw.theta);
            nw.y2 = nw.y - length * sin(nw.theta);
            nw.yMin = min(nw.y1,nw.y2);
            nw.yMax = max(nw.y1,nw.y2);
            
            nw.limsMat = [nw.xMin nw.xMax; nw.yMin nw.yMax];
            
            nw.aRow = [-tan(nw.theta) 1];
            nw.bRow = [(-nw.x * tan(nw.theta) + nw.y)];
        end
    end
end
