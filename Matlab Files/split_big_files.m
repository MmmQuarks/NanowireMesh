
for nwind = 1:10 %nanowire index counter
    clear junctionCell networkCell networkPropertiesCell 
    fileName = strcat('/Users/adamtrebach/Documents/Research/TPV/2D_nanowires/3D_nanowires/Data/nw',num2str(nwind),'.mat');
    load(fileName)

    % add resistances to junctions data
%     for i = 1:length(junctionCell)
%         junctionCell{i} = [junctionCell{i}, resistancesCell{i}];
%     end

    %splitting files
    jWhole = junctionCell;
    nWhole = networkCell;
    npWhole = networkPropertiesCell;

    junctionCell = junctionCell(1:2:length(junctionCell));
    networkCell = networkCell(1:2:length(networkCell));
    networkPropertiesCell = networkPropertiesCell(1:2:length(networkPropertiesCell));


    save(strcat('/Users/adamtrebach/Dropbox (MIT)/MIT/Research/Data/new/nw',num2str(2*nwind-1) ),'junctionCell','networkCell','networkPropertiesCell')

    junctionCell = jWhole(2:2:length(jWhole));
    networkCell = nWhole(2:2:length(nWhole));
    networkPropertiesCell = npWhole(2:2:length(npWhole));

    save(strcat('/Users/adamtrebach/Dropbox (MIT)/MIT/Research/Data/new/nw',num2str(2*nwind) ),'junctionCell','networkCell','networkPropertiesCell')

    clear all
    
end