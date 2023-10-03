function toSpiceTxt(pointsSpice, props, inFile, topElectrodeNum)

% % % % % % CHANGED POINTS --> POINTSSPICE BECAUSE OF OTHER CHANGES BELOW

    % remove any junctions which are not on the percolating cluster (or are
    % otherwise useless)
    
    %note: infile is the name of the spice netlist we're writing.
    %outfile is the name of the file THAT NETLIST will generate
    
% % % % % %     %COMMENTED OUT THIS LINE VVV
% % % % % %     %pointsSpice = points( (points(:,5) == -1 & points(:,6) == -1), :);
   
% % % % % % %     %REMOVED +1 AT END BC TAKEN CARE OF IN InternalRes.m VVV
    %topElectrodeNum = max( max(pointsSpice(:,1)), max(pointsSpice(:,2)));
function toSpiceTxt(pointsSpice, props, inFile, outFile)
    %note: infile is the name of the spice netlist we're writing.
    %outfile is the name of the file THAT NETLIST will generate
    topElectrodeNum = max( max(pointsSpice(:,1)), max(pointsSpice(:,2)));
    
    % currently both top and bottom electrodes have node number 0. the
    % bottom electrode is at y = 0. Now we set all node numbers that are
    % currently zero but whose position is not y = 0 to have nodeNum =
    % topElectrodeNum    
    
    %write to file
    fid = fopen(inFile, 'wt+');
    
    
    %print structure of network properties array to file in comment
    %{'Length', 'width' , 'nwLength', 'percolationMultiple',...
                    %'numberDensity','junctionResistanceMean', 'junctionResistanceSD'};
    fprintf(fid,'%.5e_%.5e_%.5e_%.5e_%.5e_%.5e_%.5e\n', props);
    fprintf(fid, '* length width nwLength percolationMultiple numberDensity junctionResistanceMean junctionResistanceSD\n');
    
    %insert electrodes
    fprintf(fid,'vin %d 0 DC 1\n', topElectrodeNum);
    
    %make version of points matrix with resistor numbers added to it and
    %with only the columns spice needs in the order it wants them
    % note that this is transposed because matlab reads things in fprintf
    % by going down columns like an idiot
    pointsSpice = [ transpose(1:length(pointsSpice)), pointsSpice(:,[1,2,7])]';
    
    %resistor_id node1 node2 resistance
    %fprintf(fid, 'r%d\n %d\n %d\n %f\n',pointsSpice);
    fprintf(fid, 'r%d %d %d %.5e\n',pointsSpice);
    fprintf(fid,'.op\n.end');
    fclose(fid);
end
