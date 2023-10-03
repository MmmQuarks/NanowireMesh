% to spice main

randomize = true;
for nwind = 1:20
    clear junctionCell networkCell networkPropertiesCell 
    mfile = strcat('/Users/adamtrebach/Documents/Research/TPV/2D_nanowires/3D_nanowires/Data/nw',num2str(nwind),'.mat');
    load(mfile)
    for i = 1:length(junctionCell)
        if(randomize)
            rsdRange = 0.05:0.05:0.5; %the range of relative standard deviations
            for sdc = 1:5
                % editing the file names so we don't overwrite previous
                % stuff (and so we can tell these are randomized files)
                inFile = strcat('/Users/adamtrebach/Documents/Research/TPV/2D_nanowires/3D_nanowires/Data/spice_files/nw',num2str(nwind),'_spice_',num2str(i));
                outFile = strcat('/home/amwt/TPV/3D_nanowires/data/nw1_',num2str(i));
                inFile = strcat(inFile,'_normrnd_',num2str(sdc));
                outFile = strcat(outFile,'_normrnd_', num2str(sdc));
                
                % drawing a random standard deviation from the list of
                % them we've allowed
                sdIndex = round(1 + (length(rsdRange) - 1) * rand());
                sd = rsdRange(sdIndex) * networkPropertiesCell{i}(6);
                
                %writing this value back to the network properties cell to
                %be used when writing to spice
                networkPropertiesCell{i}(7) = sd;
                
                %generating a list of resistances according to the proper
                %distribution
                mu = networkPropertiesCell{i}(6);
                
                rlist = normrnd(mu, sd, length(junctionCell{i}), 1);
                
                % make sure no resistances are negative
                while min(rlist) <0 
                    rlist( rlist<0) = normrnd(mu, sd, length(rlist(rlist<0)), 1);
                end
                
                junctionCell{i}(:,7) = rlist;
                toSpiceTxt(junctionCell{i}, networkPropertiesCell{i}, inFile,outFile);
                
            end
            
        else
    %         inFile = '/Users/adamtrebach/Documents/Research/TPV/2D_nanowires/3D_nanowires/Data/nw1_spice/';
    %         inFile = strcat(inFile, 'nw1_',num2str(i));
            inFile = strcat('/Users/adamtrebach/Documents/Research/TPV/2D_nanowires/3D_nanowires/Data/nw',num2str(nwind),'_spice_',num2str(i));
            outFile = strcat('/home/amwt/TPV/3D_nanowires/data/nw1_',num2str(i));
            toSpiceTxt(junctionCell{i}, networkPropertiesCell{i}, inFile,outFile);
        end
    end
end