
%network matrix format
%[wire_number xstart ystart xfin yfin intercept slope cluster_1 cluster_2 cluster_3 intercept_number]

%junction matrix format
% [wireNum1, wireNum2, junction X, Junction Y,  clusterNumWire1, clusterNumWire2, junctionResistance]


wireNumber = transpose(0:9);
R = transpose(1:12);


% note: this will fail if you swap the order of the wires if one of the
% wires in the electrode. otherwise order doesn't matter.
junctions = [0 1 0 1 0 0 R(1);
    0 2 0 0 0 0 R(2);
    3 1 .2 1 0 0 R(3);
    2 5 .3 0 0 0 R(4);
    3 4 .4 .6 0 0 R(5);
    3 5 .3 .4 0 0 R(6); 
    3 6 .5 .5 0 0 R(7);
    6 7 .7 .5 0 0 R(8);
    4 6 .5 .7 0 0 R(9);
    7 9 1 .4 0 0 R(10);
    7 8 .8 .3 0 0 R(11);
    9 8 1 0 0 0 R(12)];


v = ( (junctions(:,1) == 0) | (junctions(:,2) == 0)); % vector where 1 means connected to left electrode

% the minus one at the end is so we don't include the wires touching the
% left or right electrodes
A = zeros( length(junctions), length(union(junctions(:,1),junctions(:,2)))-2 );

wireInd = junctions(:,1:2);

for( i = 1:length(A) )
    w1 = wireInd(i,1);
    w2 = wireInd(i,2);
    
    % only adds element to matrix if it is not on electrodes
    if( w1 >0 & w1 < max(wireNumber))
        A(i, wireInd(i,1)) = 1;
    end
    if (w2 > 0 && w2 < max(wireNumber))
        A(i, wireInd(i,2)) = -1;
    end
end

% matrix with 1/resistances on main diagonal
K = diag( ones(length(junctions),1)./junctions(:,7) );


x = linsolve(transpose(A) * K * A, transpose(A) * K * v);


% find junctions toiuching right (ground) electrode
rightNum = max(wireNumber);
groundResistors = junctions( junctions(:,1) == rightNum | junctions(:,2) == rightNum, :);

currents = x(min(groundResistors(:,1:2),[],2)) ./groundResistors(:,7);


% the index in the x vector corresponds to the wire number at that voltage.


totalCurrent = abs(sum(currents));

totalResistance = 1/totalCurrent

    