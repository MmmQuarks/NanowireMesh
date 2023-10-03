classdef junctionObject
	properties
		x	% junction coordinates		
		y
        
        % wires it's connected to
        node1
        node2
        
        R   % junction Resistance

		
	end

	methods
        function junction = junctionObject(node1,node2,x,y,R) %class constructor
            junction.node1 = node1;
            junction.node2 = node2;
            junction.x = x;
            junction.y = y;
            junction.R = R;
        end
	end
end
