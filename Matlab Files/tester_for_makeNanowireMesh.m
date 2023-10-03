import nanowireObject
import junctionObject


density = 5.64*2.5/4;

nwLength = 2;

lengthRatio = 5; % ratio of sample length to nanowire length

junctionResistance = 10;
junctionResistanceSD = 2;

makePlot = true; %bool to control whether we make a plot
surfaceX = lengthRatio * nwLength; %length of sample in um
surfaceY = lengthRatio * nwLength;
surfaceXMin = 0.5 * surfaceX;
surfaceXMax = 1.5 * surfaceX;
surfaceYMin = 0.5 * surfaceY;
surfaceYMax = 1.5 * surfaceY;
    
tic

% we generate nanowires over a larger surface than we plane to probe to
% avoid weird edge effects
[nwArray, junctions, electrodes] = makeNanowireMesh( surfaceX,surfaceY, nwLength,density, junctionResistance, junctionResistanceSD);


if makePlot
    hold on
    scatter([junctions.x], [junctions.y] ,50,'MarkerFaceColor','red','MarkerEdgeColor', 'red')
    plot([nwArray.x1;nwArray.x2],[nwArray.y1;nwArray.y2],'Color','blue')
    edgesX = [surfaceXMin surfaceXMin surfaceXMax surfaceXMax; surfaceXMin surfaceXMax surfaceXMax surfaceXMin];
    edgesY = [surfaceYMin surfaceYMax surfaceYMax surfaceYMin; surfaceYMax surfaceYMax surfaceYMin surfaceYMin];
    plot(edgesX,edgesY,'Color','black','LineWidth',1)
    scatter(electrodes(:,3), electrodes(:,4),50,'MarkerFaceColor','blue','MarkerEdgeColor','green')
    hold off
end

toc