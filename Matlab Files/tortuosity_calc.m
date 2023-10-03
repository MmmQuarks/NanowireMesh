sys_size = 100;

[network point] = thomasMakeNanowireMesh(sys_size,sys_size,10, 3);

point = [point ones(length(point),1)];
%points structure: [wire1 wire2 Xint Yint cluster1 cluster2]

points_internalres = InternalRes(point);

%kludge to find top electrode number
topCandidates = points_internalres( (points_internalres(:,4) == sys_size), :);

topElectrodeNum = intersect(topCandidates(1,1:2),topCandidates(2,1:2));

props = [0 0 0 0 0 0 0 ];
toSpiceTxt(points_internalres, props, "spice_test2.cir", topElectrodeNum)

