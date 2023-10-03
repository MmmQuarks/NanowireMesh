function coeffs = splineSequence(points)
    % points will be a matrix of 3 rows and at least 2 columns (should 
    % put in error handling here)
    
    %put points in order by x values which are in first row
    [~, ind] = sort(points(1,:));
    points = points(:,ind);
    
    numPoints = size(points);
    numPoints = numPoints(2); %(num points is number of cols)
    
    coeffs = [];
    
    for nn = [1:(numPoints - 1)]
        coeffs = [coeffs singleSpline(points(:, nn:nn+1) ) ];
    end
        
    
    coeffs;
    
end

function coeffs = singleSpline(points)
    %labeling points for clarity
    
    x1 = points(1,1);
    x2 = points(1,2);
    y1 = points(2,1);
    y2 = points(2,2);
    z1 = points(3,1);
    z2 = points(3,2);
    
    %each curve is parametrized by t. t = 0 gives the first point.
    % t= 1 gives the second point. these coefficients are labeled such
    % coeffs.zt2 is the coefficient of the t^2 term in the z equation. and
    % etc.
    coeffs.xt0 = x1;
    coeffs.xt1 = x2 - x1;
    coeffs.yt0 = y1;
    coeffs.yt1 = y2 - y1;
    coeffs.zt0 = z1;
    coeffs.zt1 = 0;
    coeffs.zt2 = -3*(z1-z2);
    coeffs.zt3 = 2*(z1 - z2);
    
    
end



