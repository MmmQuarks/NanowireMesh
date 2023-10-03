import nanowireObject
import junctionObject

makePlot = false; %bool to control whether we make a plot

%write nodes data to file
jctFile = fopen('junctions.txt','w');


fprintf(jctFile, 'wire_density|wire_length|junction_density|w_dense*L^2|j_dense*L^2\n');

networkStats = []; % array with columns as above

%setting lengths and densities to be probed

densityVals = repelem(0.21:0.1:3.51, 5);
lengthVals = repelem(0.3:0.1:2, 50);

trialCounter = 1;
numTrials = length(densityVals)*length(lengthVals);

for dd = 1:length(densityVals)      
    for lengthCounter = 1:length(lengthVals) 
        surfaceX = 5; %length of sample in um
        surfaceY = 5;
        density = densityVals(dd);
        
        fprintf("Trial %d/%d, rho = %f, l = %f\n", trialCounter, numTrials, density, lengthVals(lengthCounter));
        trialCounter = trialCounter + 1;
            
        numWires = round(2* surfaceX * 2* surfaceY * density); % we generate 
        %nanowires over a larger surface than we plan to probe to
        % avoid weird edge effects. 
        nwLength = lengthVals(lengthCounter);
        % calc boundaries of plane 
        surfaceXMin = 0.5 * surfaceX;
        surfaceXMax = 1.5 * surfaceX;
        surfaceYMin = 0.5 * surfaceY;
        surfaceYMax = 1.5 * surfaceY;



        % we generate nanowires over a larger surface than we plane to probe to
        % avoid weird edge effects
        nwArray = makeNanowireMesh(2* surfaceX,2*surfaceY, nwLength,density);


        %find nodes/junctions

        junctions = []; %each row in this will be the indices for each
                        %wire connection

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
                                junctions = [junctions; junctionObject(ll, uu, X(1), X(2))];
                            end
                        end
                    end
                end

            end 
        end

        %make sure all lists of connected wires are complete
        for jj = 1:length(junctions)
           thisCluster = [];

           nw = nwArray(junctions(jj).node1);
           thisCluster = union(thisCluster, nw.connectedWires);

           for clusterCounter = 1:length(thisCluster)

           end

        end

        %remove most disconnected junctions from array
        %usefulNwIndices = unique([junctions(:,1);junctions(:,2)]);
        %nwArray = nwArray(usefulNwIndices);

        %plot to sanity check if makePlot = true

            if makePlot
                hold on
                scatter([junctions.x], [junctions.y] ,40,'MarkerFaceColor','white','MarkerEdgeColor', 'red')
                plot([nwArray.x1;nwArray.x2],[nwArray.y1;nwArray.y2],'Color','blue')

                hold off
            end


        junctionDensity = length(junctions)/ (surfaceX * surfaceY);
        networkStats = [networkStats; [density, lengthVals(lengthCounter), junctionDensity] ];


        % for nn = 1:length(junctions)
        %     fprintf(jctFile,'%d|',nn);
        %     fprintf(jctFile,'%d|', junctions(nn).node1);
        %     fprintf(jctFile,'%d|', junctions(nn).node2 );
        %     fprintf(jctFile,'%f|', junctions(nn).x );
        %     fprintf(jctFile,'%f|', junctions(nn).y );
        %     fprintf(jctFile,'%f\n',10);
        % end
    end
    
end

for nn = 1:length(networkStats)
    fprintf(jctFile, "%f|", networkStats(nn,1));
    fprintf(jctFile, "%f|", networkStats(nn,2));
    fprintf(jctFile, "%f|", networkStats(nn,3));
    fprintf(jctFile,"%f|", networkStats(nn,1).*networkStats(nn,2).*networkStats(nn,2));
    fprintf(jctFile,"%f\n", networkStats(nn,3).*networkStats(nn,2).*networkStats(nn,2));
end
    
    
    
    
