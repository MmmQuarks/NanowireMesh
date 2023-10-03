function outVector = rodriguesRotate(inVector, rotationAxis, angle)
	%implementation of the rodrigues rotation
	if ~(0.999 < norm(rotationAxis) < 1.001) %normalizes vector if it isn't normalized
		rotationAxis = 1/norm(rotationAxis) * rotationAxis;
	end

	outVector = cos(angle) * inVector + sin(angle) * cross(rotationAxis, inVector) + (1 - cos(angle)) * dot(rotationAxis, inVector) * rotationAxis;

end
