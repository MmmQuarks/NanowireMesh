
%test of spline sequence

numPoints = 5;

m = (rand - 1) * 10;
b = (rand - 1) * 10;

xPoints = (rand(1, numPoints)-1)*10;
yPoints = m * xPoints + b;
zPoints = rand(1, numPoints);
points = [xPoints; yPoints; zPoints];

x1 = points(1,1);
x2 = points(1,2);
y1 = points(2,1);
y2 = points(2,2);
z1 = points(3,1);
z2 = points(3,2);

coeffs = splineSequence(points);

xcurve = [];
ycurve = [];
zcurve = [];


t = 0:0.01:1;

for nn = [1:(numPoints-1)]
   xcurve = [xcurve (coeffs(nn).xt0 + coeffs(nn).xt1*t)];
   ycurve = [ycurve (coeffs(nn).yt0 + coeffs(nn).yt1*t)];
   zcurve = [zcurve (coeffs(nn).zt0 + coeffs(nn).zt1 * t + coeffs(nn).zt2 * t.^2 + coeffs(nn).zt3 *t.^3)];
end


%below is the form of the equations
% xcurve = x1 + (x2-x1)*t;
% ycurve = y1 + (y2-y1)*t;
% zcurve = 2*(z1 - z2)*t.^3 + -3*(z1-z2) * t.^2 + z1;

hold on

box on
plot3(xcurve,ycurve,zcurve)

[~, ind] = sort(points(1,:));
points = points(:,ind);

plot3(points(1,:),points(2,:),points(3,:))

hold off

