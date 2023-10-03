import nanowireObject
import junctionObject


tic


nwLength = 2;
lengthRatio = 5; % ratio of sample length to nanowire length


%critical percolation threshold is 5.64/l^2
nc = 5.64/nwLength^2;

densityValues = 1:.5:2.5 * nc;
densityValues = repelem(densityValues,5);



junctionResistance = 10;
junctionResistanceSD = 1:.25:2.5;
junctionResistanceSD = repelem(junctionResistanceSD, 8);

surfaceX = lengthRatio * nwLength; %length of sample in um
surfaceY = lengthRatio * nwLength;

junctionsData = [];
electrodesData = [];
trialNum = 1;

for densityCounter = densityValues
    for sdCounter = junctionResistanceSD
        [nwArray, junctions, electrodes] = makeNanowireMesh( surfaceX,surfaceY, nwLength,densityCounter, junctionResistance, sdCounter);

        electrodesData = [electrodesData ; [trialNum*ones(length(electrodes),1) electrodes]];
        numJunctions = length(junctions);
        % junction matrix has this structure
        % trialNum | junctionResistanceMean | junctionResistanceSD |
        % junctionResistanceValue | node1 | node2 | xPos | yPos
        thisJunctionData = [trialNum * ones(numJunctions,1) , junctionResistance * ones(numJunctions,1) , sdCounter * ones(numJunctions,1)];
        thisJunctionData = [thisJunctionData , transpose([junctions.R; junctions.node1; junctions.node2; junctions.x; junctions.y])];
        junctionsData = [junctionsData ; thisJunctionData];

        trialNum = trialNum + 1;
    end
end

electrodesData = array2table(electrodesData,'VariableNames',{'trialNum' 'wireNum' 'electrodeNum' 'xPos' 'yPos'});
%electrodes matrix has these properties
% wireNumber| electrodeNumber | xPos | yPos

%electrodes output file will have following data
%trialNumber | wireNumber | electrodeNumber | xPos | yPos

junctionsData = array2table(junctionsData, 'VariableNames',{'trialNum' 'junctionResistanceMean' 'junctionResistanceSD' 'junctionResistanceValue' 'node1' 'node2' 'xPos' 'yPos'});

save(strcat('nw',date),'electrodesData','junctionsData');

toc