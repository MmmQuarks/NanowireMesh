%this function calculations tortuosity of a network
function tortuosity = ContourCalculations(C, h) %inputs are a (contour plot data, contour object)
    bins = [];
    maxCon = 1;
    size = length(h.LevelList);
    for i = 1:size  %sort the contour data by level, extra logic added in the case of multiple lines of the same level
        temp = find(C(1,:) == h.LevelList(i));
        if length(temp) > maxCon        %Matrix dimensional agreement logic
            maxCon = length(temp);
            if isempty(bins)
                bins = temp;
            else
                bins = horzcat(bins,zeros(length(bins(:,1)),maxCon-length(bins(1,:))));
                bins = vertcat(bins, temp);
            end
        else
            temp = horzcat(temp,zeros(1,maxCon-length(temp)));
            bins = vertcat(bins, temp);
        end
        clear temp
    end
    [real, numpts] = contourlengths(C); 
    bins = sortrows(bins(bins ~=0)); %remove zeros and sort rows linear by index
    sepVals = zeros(length(bins),3); 
    for i = 1:length(bins)
        if i == length(bins)
        sepVals(i,:) = [sqrt(sum((C(:,length(C)) - C(:,bins(i)+1)).^2)), real(i), numpts(i)];
        else
        sepVals(i,:) = [sqrt(sum((C(:,bins(i+1) - 1) - C(:,bins(i)+1)).^2)), real(i), numpts(i)];
        end
    end
    sepVals = sepVals((sepVals(:,1) ~= 0 )&(sepVals(:,3) > 2),:);
    sepVals = sepVals(sepVals(:,2)./sepVals(:,1) < 1000,:);
    tortuosity = sum((sepVals(:,2)./sepVals(:,1)).*sepVals(:,3))/(sum(sepVals(:,3)));