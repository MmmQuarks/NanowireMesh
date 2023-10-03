%this function creates a mesh of equipotentials and calculates tortuosity
%wireVolt is the csvread output of the spice voltages, presumed to be done
%in the main script that calls each function
%function [averageTort, C, h, volts]= tortuosityCalculation(points_internalres, outfile, meshSize)
function [averageTort, C, h, volts]= tortuosityCalculation(wireVolt, meshSize)
    %wireVolt = csvread(strcat(outfile, '_voltages.csv'));
    %wireVolt = horzcat(wireVolt,zeros(length(wireVolt(:,1)),2));
    wireVolt = sortrows(wireVolt,1); %order by wire number
    wireVolt = wireVolt( wireVolt(:,3) >= 0, :);
    wireVolt = wireVolt( wireVolt(:,3) <= 700, :);
    wireVolt = wireVolt( wireVolt(:,4) >= 0, :);
    wireVolt = wireVolt( wireVolt(:,4) <= 700, :);
    maxVals = max(wireVolt(:,3:4)); %x,y bounds
    voltages = zeros(meshSize);
    %find the corresponding midpoint to each wire/voltage pair generated in
    %spice
    % correct format for wirevolt should be
    % wireNum | voltage | x | y
%     maxnum = max(wireVolt(:,1));
%     for i = 1:maxnum
%         temporary = points_internalres(points_internalres(:,1) == i | points_internalres(:,2) == i,:);
%         temporary(temporary(:,9) == 0,:) = [];
%         wireVolt(i,3:4) = mean(temporary(:,3:4),1);
%         clear temporary
%     end
%     wireVolt(:,3) = discretize(wireVolt(:,3),0:maxVals(1)/meshSize:maxVals(1));
%     wireVolt(:,4) = discretize(wireVolt(:,4),0:maxVals(2)/meshSize:maxVals(2));
%     wireVolt = array2table(wireVolt);
%     heatmap(wireVolt,'wireVolt3','wireVolt4','ColorVariable','wireVolt2');

    for i = 1:meshSize %average all voltages within a half-square about each mesh point
       for j = 1:meshSize
          clustermat = wireVolt( ( abs((wireVolt(:,3)-(i*maxVals(1)/meshSize))) <= maxVals(1)/(2*meshSize)) & (abs((wireVolt(:,4)-(j*maxVals(2)/meshSize))) <= maxVals(2)/(2*meshSize)),2);
          voltages(j,i) = mean(clustermat);
          clear clustermat
       end
    end
    
    surf(voltages);
    savefig('voltagesurf');
    
    [C, h] = contour(voltages, meshSize); %create contour plot and call calculation function
    averageTort = ContourCalculations(C, h);
    volts = voltages;
savefig('contourplotfigure');