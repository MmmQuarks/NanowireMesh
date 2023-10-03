function [network point] = thomasMakeNanowireMesh(Length, width, NW_length, percolationMultiple)
    %
    % This file is used to generate the percolation network.
    % The wires and junctions(single junctions and the composited junctions) is ready to be regenerated in Comsol later
    % with another code called "Mat_2_Comsol_RandomJunction_solved"

    % To change the systemsize, please change the Length and width value in line 14 and 15
    % To change the network denstiy, please change the factor in front of Nc in line 46
    % The step_size in line 19 can also be adjusted depending on how many wires there are in total

    % After having a satisfying network, please save the work space with command: save('file name')
    rng('shuffle');
    dbstop if error
    close all
    %clear all
    %Length=[50]; % the unit is um
    %width=[50];
    %NW_length=10
    angulardist=1;
    Ldist=1;
    step_size=2;
    Numsims=1;
    NW_numt=0;

    r=width/NW_length;
    Nc=5.63726*r*r+r+5.5 % the number of nanowire needed to have 50% chance of reaching percolation
    N99=5.86*r*r+r+9.9; % the number of nanowire needed to have 99% chance of reaching percolation

    disp('loop start');
    tic % start to count the elapsed time     
    for o=1:Numsims
        percolation = 0; % Trigger to end a simulation once percolation occurs this value becomes 0
        ic = 0; % iteration counter used to determine the total number of nanowires in a given network
        num = 1; % Position counter to define the location in the array which the next nanowires data should go
        clusterc = 0; % Intitial cluster numbering
        network =	[0 0 0 width 0 0.0000000001 0 -1 -1 0 0;...
            0 0 Length width Length -0.00000000001 Length -1 0 -1 0];
        % network = [wire_number xstart ystart xfin yfin intercept slope cluster_1 cluster_2 cluster_3 intercept_number]
        % this first two rows will be covered later
        delta2 = 1;
        Xint = 0;
        Yint = 0;
        num1 = 1;
        point=0; % Point is the interception point matrix, A=[wire1 wire2 Xint Yint cluster1 cluster2]
        point_p=0;
        int_N=zeros(1000*floor(N99),1); % count how many interception point that each wire has

    while NW_numt<percolationMultiple*Nc
        percolation = 0; % Trigger to end a simulation once percolation occurs this value becomes 0

    while percolation == 0;
        ic = ic + 1;  % Count current iteration
        NW_numt = ic * step_size; % Calculate total nanowire number

        if rem(ic,200)==0
            disp(ic);
        end
        NW_numt = ic * step_size; % Check wether this line is useless...

    % --------------------Generate attributes for the nanowires about to be added to the network.-----------------
        L = NW_length + (NW_length * Ldist)/4 * randn(step_size,1); % lengths
        angle = -90 + angulardist*(360 * rand(step_size,1)); % angle of each of the new wires
        originX = width * rand(step_size,1); % start of nanowire x position
        originY = Length * rand(step_size,1); % start of nanowire y position
        endX = originX + L .* cosd(angle); % End of nanowire x position
        endY = originY + L .* sind(angle); % End of nanowire y position

        cenX = originX + 0.5 .* L .* cosd(angle); % center of nanowire x position
        cenY = originY + 0.5 .* L .* sind(angle); % center of nanowire y position
        clear angle L;
        slope = (endY - originY) ./ (endX - originX);
        Cintercept = endY - slope .* endX;
        % Non iterative generation of cluster numbers for new nanowires
        clustercon(1:step_size,1) = 1:step_size; % cluster number
        clustercon(1:step_size,2) = 1:step_size; % connection to bottom
        clustercon(1:step_size,3) = 1:step_size; % connection to top
        cluster =clustercon + clusterc;
        clusterc = clusterc + step_size;

    % --------------------------Check for NNWs extending out from edges ----------------------------------------

        % Bottom and left edges
        for i = 1:step_size
            if (endX(i,1) < 0)
                endX(i,1) = 0;
                endY(i,1) = Cintercept(i,1);
            end
            if (endY(i,1) < 0)
                endY(i,1) = 0;
                endX(i,1) = (-1 * Cintercept(i,1)) / slope(i,1);
                if endX(i,1)<width || endX(i,1)==width 
                    cluster(i,2) = -1; % Sets nanowire that cross the bottom of the network to the same cluster number
                    point(delta2,3) = -Cintercept(i)/slope(i); %Xint
                    point(delta2,4) = 0; %Yint
                    point(delta2,1)=0;   % intercept with x axis
                    point(delta2,2)=i+(ic-1)*step_size;   % wire_i+(ic-1)*step_size intercept
                    delta2 = delta2+1;
                    int_N(i+(ic-1)*step_size,1)=int_N(i+(ic-1)*step_size,1)+1; % the relative wire has 
                end                                                            % one more intercept point
            end
        end
        clear i;

        % Top and right edges
        for i = 1:step_size
            if (endX(i,1) > width)
                endX(i,1) = width;
                endY(i,1) = slope(i,1) * endX(i,1) + Cintercept(i,1);
            end
            if (endY(i,1) > Length)
                endY(i,1) = Length;
                endX(i,1) = (endY(i,1) - Cintercept(i,1)) / slope(i,1);
                cluster(i,3) = -1; % Sets nanowire that cross the top of the network to the same cluster number
                point(delta2,3) = (Length-Cintercept(i,1))/slope(i,1); % Xint
                point(delta2,4) = Length; % Yint
                point(delta2,1) = i+(ic-1)*step_size; % wire_i+(ic-1)*step_size intercept
                point(delta2,2) = 0; % intercept with the top
                delta2 = delta2+1;
                int_N(i+(ic-1)*step_size,1)=int_N(i+(ic-1)*step_size,1)+1;
            end
        end
        clear i;

        %---------------------Calculate any intersections between the new nanowires x1=(c2-c1)/(m1-m2)----------
        for i=1:step_size
            for g=1:step_size-i
                XintTEMP(delta2,1) = -(Cintercept(i)-Cintercept(step_size+1-g))/(slope(i)-slope(step_size+1-g));
                if ...
                        XintTEMP(delta2,1)>=min(originX(i,1),endX(i,1))...
                        && XintTEMP(delta2,1)<=max(originX(i,1),endX(i,1))...
                        && XintTEMP(delta2,1)>=min(originX(step_size+1-g,1),endX(step_size+1-g,1))...
                        && XintTEMP(delta2,1)<=max(originX(step_size+1-g,1),endX(step_size+1-g,1));
                    YintTEMP(delta2,1) = slope(i,1)*XintTEMP(delta2,1)+Cintercept(i,1);
                    if ...
                            YintTEMP(delta2,1)>=min(originY(i,1),endY(i,1))...
                            && YintTEMP(delta2,1)<=max(originY(i,1),endY(i,1))...
                            && YintTEMP(delta2,1)>=min(originY(step_size+1-g,1),endY(step_size+1-g,1))...
                            && YintTEMP(delta2,1)<=max(originY(step_size+1-g,1),endY(step_size+1-g,1));

                        point(delta2,3) = XintTEMP(delta2,1); % Xint
                        point(delta2,4) = YintTEMP(delta2,1); % Yint
                        point(delta2,1) = i+(ic-1)*step_size; % wire_i intercept with wire_(step_size+1-g)
                        point(delta2,2) = step_size+1-g+(ic-1)*step_size;
                        delta2 = delta2+1;

                        int_N(i+(ic-1)*step_size,1)=int_N(i+(ic-1)*step_size,1)+1;
                        int_N(ic*step_size-g+1,1)=int_N(ic*step_size-g+1,1)+1; %ic*step_size-g+1 corresponds to the
                        %wire number in the network matrix

                        new_cluster=min(cluster(i,1),cluster(step_size+1-g,1));
                        old_cluster_1=cluster(i,1);
                        old_cluster_2=cluster(step_size+1-g,1);
                        for d=1:step_size
                            if cluster(d,1)==old_cluster_1 || cluster(d,1)==old_cluster_2
                                cluster(d,1)=new_cluster;
                            end
                        end
                    end
                end
            end
        end
        clear belta i m n g
        %------------------------Calculate any intersections between the new nanowires and the network-----------------
        for i=1:step_size
            [m,n]=size(network);
            for g = 1:m
                XintTEMP(delta2,1) = -(Cintercept(i)-network(g,6))/(slope(i)-network(g,7));
                if ...
                        XintTEMP(delta2,1)>=min(originX(i,1),endX(i,1))...
                        && XintTEMP(delta2,1)<=max(originX(i,1),endX(i,1))...
                        && XintTEMP(delta2,1)>=min(network(g,2),network(g,4))...
                        && XintTEMP(delta2,1)<=max(network(g,2),network(g,4));
                    YintTEMP(delta2,1) = slope(i)*XintTEMP(delta2,1)+Cintercept(i);
                    if ...
                            YintTEMP(delta2,1)>=min(originY(i,1),endY(i,1))...
                            && YintTEMP(delta2,1)<=max(originY(i,1),endY(i,1))...
                            && YintTEMP(delta2,1)>=min(network(g,3),network(g,5))...
                            && YintTEMP(delta2,1)<=max(network(g,3),network(g,5));

                        point(delta2,3) = XintTEMP(delta2,1); % Xint
                        point(delta2,4) = YintTEMP(delta2,1); % Yint
                        point(delta2,1) = i+(ic-1)*step_size; % wire_i intercept with wire_(step_size+1-g)
                        point(delta2,2) = g;
                        delta2 = delta2+1;

                        int_N(i+(ic-1)*step_size,1)=int_N(i+(ic-1)*step_size,1)+1;
                        int_N(g,1)=int_N(g,1)+1;

                        new_cluster = min(cluster(i,1),network(g,8));
                        old_cluster_1=network(g,8);
                        old_cluster_2=cluster(i,1);

                        for d = 1:m
                            if network(d,8) == old_cluster_1 || network(d,8) == old_cluster_2
                                network(d,8) = new_cluster;
                            end
                        end
                        for h=1:step_size
                            if cluster(h,1)==old_cluster_1 || cluster(h,1) == old_cluster_2
                                cluster(h,1)=new_cluster;
                            end
                        end
                    end
                end
            end
        end
        clear belta i g m n

        % Combine nanowire attributes into the network

        network (num1:num1+step_size-1,1) = num1:num1+step_size-1;
        network (num1:num1+step_size-1,2) = originX;
        network (num1:num1+step_size-1,3) = originY;
        network (num1:num1+step_size-1,4) = endX;
        network (num1:num1+step_size-1,5) = endY;
        network (num1:num1+step_size-1,6) = Cintercept;
        network (num1:num1+step_size-1,7) = slope;
        network (num1:num1+step_size-1,8) = cluster(:,1);
        network (num1:num1+step_size-1,9) = cluster(:,2);
        network (num1:num1+step_size-1,10) = cluster(:,3);
        num1 = num1 + step_size;

        %------For each group of cluster, decide whether the cluster is percolating -------------------------------------------

        network=sortrows(network,8); % rearrange the matrix according to the value of cluster_1
        cluster_1=unique(network(:,8)); % get the different values in the column
        [f v]=hist(network(:,8),cluster_1); % value v in matrix cluster_1 appears f times in column network(:,7)
        group_cluster=[f' v]; % we have f' here because f is a row vector, and v is a column vector

        % Attribute the proper cluster_2(connect to the bottom) and cluster_3(connect to the top) value to the wires
        [m dummy]=size(group_cluster);
        c=1;
        for i=1:m
            d=min(network(c:c+group_cluster(i,1)-1,9));
            e=min(network(c:c+group_cluster(i,1)-1,10));
            network(c:c+group_cluster(i,1)-1,9)=d;
            network(c:c+group_cluster(i,1)-1,10)=e;
            c=c+group_cluster(i,1);
            if d==-1 && e==-1
                percolation=1;
            end
        end
        clear i m c dummy

        network=sortrows(network,1);
        num = num + 1;
    end
    end

    network (:,11) = int_N(1:ic*step_size,1);

    %-------------------------------------Remove some useless wires from the percolation network--------------------
    [m dummy]=size(group_cluster);
    network=sortrows(network,8);
    c=1;
    for i=1:m
        for g=c:c+group_cluster(i,1)-1
            if any(network(c:c+group_cluster(i,1)-1,9)==-1)...
               && any(network(c:c+group_cluster(i,1)-1,10)==-1) % make sure the cluster is percolating
                n=0; % n is the wire has only one interception in the current cluster; initialize n
                n1=0;
                int_N=network(:,11);
                n1=find(int_N(c:c+group_cluster(i,1)-1)==1);
                n=n1+c-1; % the wire has only one interception in the current cluster
                if ~isempty(n)
                    network(n,9)=-0.5;
                    network(n,10)=-0.5;
                    network(n,11)=network(n,11)-1;
                    [a dummy]=size(n);
                    for f=1:a
                        [m2 n2]=find(point(:,1:2)==network(n(f,1),1)); % help to find out which wire intercepts with wire_n
                        network=sortrows(network,1);
                        [m3 dummy]=size(m2);
                        for b=1:m3
                            if n2==1
                                w=point(m2(b,1),2); % w is the wire intercept with wire_n
                            else
                                w=point(m2(b,1),1);
                            end
                            network(w,11)=network(w,11)-1;
                        end
                        network=sortrows(network,8);
                    end
                end
            end
        end
        c=c+group_cluster(i,1);
    end
    clear m c g n n1 a m2 n2 b

    network=sortrows(network,1);
    %-------------------------------------- Interception points part -----------------------------------------------------
    % Decide the cluster number of the interception points
    point_p=zeros(1,6);
    num3=1;
    [a b]=size(point);
    for i=1:a
        if all(point(i,1:2)~=0)  % when the wire connect to the bottom or top, the corresponding cluster number should be -1
            % and cannot be be attributed by [cluster] matrix
            point(i,5)=max(network(point(i,1),9),network(point(i,2),9));
            point(i,6)=max(network(point(i,1),10),network(point(i,2),10));
        elseif point(i,1)==0
            point(i,5)=-1;
            point(i,6)=network(point(i,2),10);
        elseif point(i,2)==0
            point(i,6)=-1;
            point(i,5)=network(point(i,1),9);
        end

        if point(i,5)==-1 && point(i,6)==-1
            point_p(num3,:)=point(i,:);
            num3=num3+1;
        end
    end
    clear belta i m n a b

    end
    toc
    %disp('loop stops');

    %--------------------------------------------------- plot --------------------------------------------------------------
%     hold all
%     plot([0 width],[0 0],'b','LineWidth',3)
%     plot([0 width],[Length Length],'b','LineWidth',3)
    [a b]=size(network);
    network_p=zeros(1,11); % percolation network
    delta3=1;
    for i=1:a
        if network(i,9)==-1 && network(i,10)==-1
%            plot([network(i,2) network(i,4)],[network(i,3) network(i,5)],'b')
           network_p(delta3,:)=network(i,:); %[wire_number xstart ystart xfin yfin intercept slope cluster_1 cluster_2 cluster_3 intercept_number]
           delta3=delta3+1;
        else
            if network(i,9)==-0.5 && network(i,10)==-0.5
%              plot([network(i,2) network(i,4)],[network(i,3) network(i,5)],'m')
            else
%              plot([network(i,2) network(i,4)],[network(i,3) network(i,5)],'k')
            end
        end
    end
    %scatter(point_p(:,3),point_p(:,4),'.','r')
%     axis('equal')
%     axis([0,width,0,Length])
%     hold off

    %-----------------------------------------------------------------------------------------------------------------------
    %------------------------------------------------- Comsol part ---------------------------------------------------------
    %-------------------- calculate the network_comsol matrix--------------------
    [M N] = size(network_p);
    clear i
    i = sqrt(-1);
    L=network_p(:,2)+network_p(:,3)*i-(network_p(:,4)+network_p(:,5)*i);
    network_comsol(:,1) = network_p(:,1);                    % wire number in the network matrix
    network_comsol(:,2) = 1:M;                               % wire number in the network_comsol matrix
    network_comsol(:,3) = (network_p(:,2)+network_p(:,4))/2; % center_x
    network_comsol(:,4) = (network_p(:,3)+network_p(:,5))/2; % center_y
    network_comsol(:,5) = atan(network_p(:,7))*180/pi;       % angle in degree
    network_comsol(:,6) = abs(L);                            % wire length


    %--------------- Combine the junctions nearby -----------------
%     %disp('Combine the junctions nearby')
%     a = all(point_p,2);               % find in A_p the rows that don't contain 0
%     int_tempt = point_p(a,1:4);       % int_tempt=[w1_i w1_j Xint Yint dist_2_zero]; not include the junctions in the electrode
%     dist_0 = sqrt(int_tempt(:,3).^2+int_tempt(:,4).^2);
%     int_tempt(:,5) = dist_0;
% 
%     [m dummy] = size(int_tempt);
%     for i = 1:m                   % this loop is to transfer w1_i, w1_j to w2_i, w2_j
%         w1_1 = find(network_comsol(:,1)==int_tempt(i,1));
%         w1_2 = find(network_comsol(:,1)==int_tempt(i,2));
%         int_tempt(i,1) = network_comsol(w1_1,2);
%         int_tempt(i,2) = network_comsol(w1_2,2);
%     end                           % now,int_tempt=[w2_i w2_j Xint Yint distance_2_zero]
% 
%     int_tempt = sortrows(int_tempt,5);
%     int_remove = []; % store the junctions that are going to be combined and removed from the junction list
%     int_combine = {}; % in each row of the int_comsol cell, we have a matrix = [int_number w1 w2 Xint Yint]
%     range = 0.07;
%     for i = 1:m
%         int_remove_tempt = [];
%         limit = int_tempt(i,5)+range;
%     for g = i+1:m
%         if int_tempt(g,5)>limit
%             break
%         end
% 
%         dist_r =  sqrt((int_tempt(i,3)-int_tempt(g,3)).^2+(int_tempt(i,4)-int_tempt(g,4)).^2); % calculate the
%         % distance between two points
%         if dist_r<range
%             int_remove_tempt = union(int_remove_tempt,[i;g]);
%         end
%     end
% 
%     if isempty(int_remove_tempt)~=1
%         int_remove = union(int_remove,int_remove_tempt);
%         [m2 dummy] = size(int_combine);
%         p2 = 0; % a variable to tell whether int_remove_tempt is unioned in the some row of int_combine
%         for k = 1:m2
%             if any(ismember(int_remove_tempt,int_combine{k,1}(:,1))) % int_remove_tempt and int_combine{k,1} have the same element
%                 int_combine{k,1} = union(int_remove_tempt,int_combine{k,1}(:,1));
%                 int_combine{k,1}(:,2:5) = int_tempt(int_combine{k,1}(:,1),1:4); % in each row of the int_comsol cell, we have a matrix
%                                                             % the matrix = [int_number w1 w2 Xint Yint]
%                 p2 = 1;
%                 break
%             end
%         end
% 
%         if p2 == 0
%             int_combine{m2+1,1} = int_remove_tempt;
%             int_combine{m2+1,1}(:,2:5) = int_tempt(int_combine{m2+1,1}(:,1),1:4);
%             if m2>0
%                 int_com_i = int_combine{m2,1};
%                 [m3 dummy] = size(int_com_i);
%                 W = setxor(int_com_i(1,2:3),int_com_i(2,2:3));
%                 if m3 == 2 && length(W) == 2 
%                     int_com_i(3,2:3) = W;
%                     [a b] = find(int_tempt(:,1:2) == int_com_i(3,2));
%                     for j =1:length(a)
%                         if int_tempt(a(j),3-b(j)) == int_com_i(3,3)
%                             int_remove = [int_remove;a(j)];
%                             int_com_i(3,1) = a(j);
%                             int_com_i(3,4:5) = int_tempt(a(j),3:4);
%                             int_combine{m2,1} = int_com_i;
%                             break
%                         end
%                     end
%                 end
%             end
%         end
%     end
% 
%     end
%     clear i g k
% 
%     int_Comsol = int_tempt; % int_Comsol=[w2_i w2_j Xint Yint distance_2_zero]=[wire number updated] INTERSECTION MATRIX
%     int_Comsol(int_remove,:) = []; % remove the junctions that are already combined together 


    % ------ plot to examine the combined junctions ------
    % figure(2)
    % hold on
    % [m2 n2] = size(int_combine);
    % for i=1:m2
    %     scatter(int_combine{i,1}(:,4),int_combine{i,1}(:,5),'b');    
    % end
    % hold off

    %save(['sys_' num2str(length)])
end
