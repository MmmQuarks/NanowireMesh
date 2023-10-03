function [totalResistance, network, points] = analyzeCircuit(network, points)
    % this analysis assumes that the top electrode is at nonzero potential
    % and bottom electrode is at ground.
    % essentially, whichever junction(s) are at the top are considered
    % connected to the potential source. whichever junction(s) are at the
    % bottom as considered connected to ground.
    
    
    %network structure
    % [ wireNum xStart yStart xFin yFin Cintercept, slope, clust1, clust2,
    % clust3, numContacts]
    
    % points structure
    % [wire1Num wire2Num xPos, yPos, wire1ClusterNum, wire2ClusterNum]
    
    % first remove any junctions that are not in the percolating cluster
    points = points( points(:,5) == -1 & points(:,6) == -1, :);
    
    % remove any wires that are not in the percolating cluster
    network = network;    %START HERE
    network = network( (network(:,9)== -1 & network(:,10)==-1),:);
    
    %sort all wire nums such that column 1 contains the lowest wire number.
    %This is necessary for some reason to ensure everything makes sense.
    wireNumsSorted = min(points(:,1:2),[],2);
    wireNumsSorted = [wireNumsSorted, max(points(:,1:2),[],2)];
    points(:,1:2) = wireNumsSorted;
    
    % re indexing wire numbers so that matrix multiplication is easier
    % (i.e. column numbers will correspond to wire numbers)
    % [ oldWireNums, newWireNums]
    newWireNums = [network(:,1) , transpose(1:length(network(:,1))) ];

    %vector of things connected to top electrode
    v = (points(:,1) == 0 & points(:,4) == max(points(:,4)) ) * 1;
    %note that wireNum = 0 in the points matrix indicates that this is
    %either the top or bottom electrode. 
    
    % create matrix to subtract potentials of connected wires across their
    % junctions. This will multiply a vector of voltages (1 for each wire
    % in the network) so A's length should be the number of junctions and
    % its width should be the number of wires. The minus one is because the
    % wirenumber 0 is electrodes and shouldn't be included here.
  
    A = zeros( length(points), length(union(points(:,1),points(:,2)))-1 );
    
    wireInd = junctions(:,1:2);s
    for( i = 1:length(A) )
        w1 = wireInd(i,1);
        w2 = wireInd(i,2);

        % only adds element to matrix if it is not on electrodes
        if( w1 >0 & w1 < max(wireNumber))
            A(i, wireInd(i,1)) = 1;
        end
        if (w2 > 0 && w2 < max(wireNumber))
            A(i, wireInd(i,2)) = -1;
        end
    end

end