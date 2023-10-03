
% all lengths in um
Length = 200;
width = 200;


percolationMultipleRange = 5:1:20;
percolationMultipleRange = repelem(percolationMultipleRange, 10);

nwLengthRange = 10:2:40;

rMeanRange = 1;%logspace(0,6,25); %resistance Range

rSDRange = 0; % range of standard deviations of resistance



tic
%percolationMultipleRange = repelem(percolationMultipleRange,10);

junctionArray = [];
networkArray = [];
resistancesArray = [];
networkPropertiesArray = [];


parfor pc = 1:length(percolationMultipleRange) %percolation counter
    for nwc = 1:length(nwLengthRange)     
        for rc = 1:length(rMeanRange)
            for sdc = 1:length(rSDRange)
                
                %defining iteration counter for assigning values to cell arrays
                ic = (pc - 1) * length(nwLengthRange) * length(rMeanRange) * length(rSDRange) ...
                    + (nwc - 1) * length(rMeanRange) * length(rSDRange) ...
                    + (rc - 1) * length(rSDRange) ...
                    + sdc;

                
                percolationMultiple = percolationMultipleRange(pc);
                nwLength = nwLengthRange(nwc);
                rMean = rMeanRange(rc);
                rSD = rSDRange(sdc);
                
                r=width/nwLength;
                Nc=5.63726*r*r+r+5.5;
                
                
                [network, junctions] = thomasMakeNanowireMesh(Length,width,...
                nwLength, percolationMultiple);
                
                %storing the generated junctions and network in cell arrays
                junctions = [ (ic * ones(length(junctions),1)), junctions];
                junctionArray = [junctionArray ; junctions];
                
                network = [ (ic * ones(length(network),1)), network];
                networkArray = [networkArray ; network];
                
                
                %generating the junction resistance. Each row of this
                %corresponds to the same row of the junctions data
                resistances = normrnd(rMean, rSD,length(junctions),1);
                resistances = [ (ic * ones(length(resistances),1)), resistances];
                resistancesArray = [resistancesArray ; resistances];   
                
                
                
                %below are the meanings of the values in networkProperties
                %keySet = {'Iteration', 'Length', 'width' , 'nwLength', 'percolationMultiple',...
                    %'numberDensity','junctionResistanceMean', 'junctionResistanceSD'};               
                numberDensity = Nc*percolationMultiple/(Length * width);
                networkProperties = [ic, Length, width, nwLength, percolationMultiple, numberDensity, rMean, rSD];
                networkPropertiesArray = [networkPropertiesArray; networkProperties];            
            end
        end
    end
end

%storing all the data in cells at the end
for ic = 1:(length(percolationMultipleRange) * length(rMeanRange) * length(nwLengthRange) * length(rSDRange))
    % convert junctions array to cell
    thisJunctionIndices = junctionArray(:,1) == ic;
    thisJunction = junctionArray(thisJunctionIndices, :);
    junctionCell{ic} = thisJunction(:,2:length(thisJunction(1,:)));
    
    %convert network array to cell
    thisNetworkIndices = networkArray(:,1) == ic;
    thisNetwork = networkArray(thisNetworkIndices, :);
    networkCell{ic} = thisNetwork(:,2:length(thisNetwork(1,:)));
    
    %convert resistances array to cell
    thisResistancesIndices = resistancesArray(:,1) == ic;
    thisResistances = resistancesArray(thisResistancesIndices, :);
    resistancesCell{ic} = thisResistances(:,2:length(thisResistances(1,:)));

    %convert network properties array to cell (diff structure here)
    %I'm not sure what I'm doing here is necessary but I'm not sure of the
    %order things will be in in the array of properties so I want to make
    %sure they match up. Did this in different way for above. 
    thisNetworkPropertiesIndices = networkPropertiesArray(:,1) == ic;
    thisNetworkProperties = networkPropertiesArray(thisNetworkPropertiesIndices, :);
    networkPropertiesCell{ic} = thisNetworkProperties(:,2:length(thisNetworkProperties(1,:)));    
end


toc
%structure of network columns is
%1 wireNumber before non-percolating wires are removed
%2 wireNumber after non percolating wires are removed
%3 center x
%4 center y
%5 angle in degrees
%6 length

%structure of junctions columns is (I think)
%1 junction number before non-percolating wires removed
%2 junction number after non-percolating wires removed
%3 x coord
%4 y coord
%5 distance to zero


%save(strcat('/home/amwt/TPV/3D_nanowires/data/','nw',datestr(datetime)),'networkCell','junctionCell','resistancesCell','networkPropertiesCell')
save(strcat('/Users/adamtrebach/Dropbox (MIT)/MIT/Research/Data/nw',datestr(datetime)),'junctionCell','networkCell','resistancesCell','networkPropertiesCell')


