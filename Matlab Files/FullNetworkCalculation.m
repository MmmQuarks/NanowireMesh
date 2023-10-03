%Generate network and points
[network, point] = thomasMakeNanowireMesh(200,200,20,3);

%add internal resistances
point = InternalRes(point);

%infile = 'test1';
%send to spice-readable file
toSpiceTxt(point, [0 0 0 0 0 0 0], 'test1', 'test1one');

%call python script

setOfVolts = csvread([infile + '_voltages.csv']);

%perform tortuosity calc
meshSize = 50;
[tortuosity, contourData, contourObject, Voltages]=tortuosityCalculation(network,point,setOfVolts,meshSize);