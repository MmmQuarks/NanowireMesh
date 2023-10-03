function [nwArray, junctions, electrodes] = makeNanowireMesh(surfaceX,surfaceY,nwLength,density, junctionResistance, junctionResistanceSD)
    %bounds     an array containing the x and y extents of the plane
    %numWires   the total number of wires we want to deposit

    import nanowireObject
    
   
    numWires = round(2 * surfaceX * 2*surfaceY * density);
% define boundaries of plane and how many wires we want on it

    nwArray = [];

    for nn = 1:numWires
        nwArray= [nwArray, nanowireObject([0,2*surfaceX; 0, 2* surfaceY],nwLength)];
        nwArray(nn).wireNum = nn;
        nwArray(nn).connectedWires = [nn];
    end
    
    junctions = [];
    
    %the below gives the bounds of the sample region (far from edges)
    surfaceXMin = 0.5 * surfaceX;
    surfaceXMax = 1.5 * surfaceX;
    surfaceYMin = 0.5 * surfaceY;
    surfaceYMax = 1.5 * surfaceY;
    
    for uu = 2:numWires %means upper bound
        for ll = 1:(uu-1) %means lower bound
            nw1 = nwArray(ll);
            nw2 = nwArray(uu);

            A = [nw1.aRow; nw2.aRow];
            b = [nw1.bRow; nw2.bRow];
            % solve matrix equation Ax = b
            X = linsolve(A,b);

            % if this point is within wire length and sample size, add to
            % junction list
            if (X(1) < nw1.xMax && X(1) < nw2.xMax)
                if( X(1) > nw1.xMin && X(1) > nw2.xMin)
                    if( X(1) > surfaceXMin && X(1) < surfaceXMax)
                        if( X(2) > surfaceYMin && X(2)< surfaceYMax)
                            %the columns of this junctions matrix are
                            %  wire1 | wire2 | xPos | yPos |
                            nwArray(ll).connectedWires = union(nw1.connectedWires, nw2.connectedWires);
                            nwArray(uu).connectedWires = nwArray(ll).connectedWires;
                            R = normrnd(junctionResistance, junctionResistanceSD);
                            while R <= 0
                                R = normrnd(junctionResistance, junctionResistanceSD);
                            end
                                
                            % assign the junction a resistance
                            junctions = [junctions; junctionObject(ll, uu, X(1), X(2), R)];
                        end
                    end
                end
            end

        end 
    end
    
    %here we determine which wires are connected to the top and bottom
    topARow = [0 1];
    topBRow = [surfaceYMax];
    bottomARow = [0 1];
    bottomBRow = [surfaceYMin];

    electrodes = [];
    %columns: wirenumber, electrodeNumber, posX, poxY
    for uu = 1:numWires

        %does it touch bottom electrode
        thisNw = nwArray(uu);
        A = [thisNw.aRow; bottomARow];
        b = [thisNw.bRow; bottomBRow];
        X = linsolve(A,b);
        if (X(1) > surfaceXMin && X(1) < surfaceXMax)
            %makes sure wires touch electrode lines in region of sample
            if  X(2) < max(thisNw.y1, thisNw.y2) && X(2) > min(thisNw.y1 , thisNw.y2) %norm(X - [thisNw.x; thisNw.y]) <= thisNw.length/2.
                %makes sure only counts if they cross within the finite extent
                %of the nanowire
                electrodes = [electrodes; uu 1 X(1) X(2)];
            end
        end

        %does it touch top electrode
        A = [thisNw.aRow; topARow];
        b = [thisNw.bRow; topBRow];
        X = linsolve(A,b);
        if (X(1) > surfaceXMin && X(1) < surfaceXMax)
            %makes sure wires touch electrode lines in region of sample
            if  X(2) < max(thisNw.y1, thisNw.y2) && X(2) > min(thisNw.y1 , thisNw.y2) %norm(X - [thisNw.x; thisNw.y]) <= thisNw.length/2.
                %makes sure only counts if they cross within the finite extent
                %of the nanowire
                electrodes = [electrodes; uu 2 X(1) X(2)];
            end
        end


    end
end
