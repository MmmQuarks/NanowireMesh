function [clen, ptsmat, clev] = contourlengths(a);

% CONTOURLENGTHS  Length of contour lines. 
%   [clen, clev] = contourlengths(a,numlevels) computes the length of each
%   contour line in the contour matrix a (returned from executing CONTOUR
%   
%   The length of the contour line is the sum of the lengths
%   of line segments between the consecutive vertices in a contour line.
%   
%   a          -   Contour matrix returned by CONTOUR or CONTOURC
%   
%   clen       -   Vector containing the contour lines lengths
%   clev       -   Vector containing the corresponding contour line levels

len = zeros(1,length(a));    % Create a matrix to hold the length of each contour line
lev = len;  % Create another matrix of the same size to hold the corresponding level

linenum = 1;
k = 0;

while linenum <= length(a)
    k = k+1;
    lev(k) = a(1,linenum);  % lev(k) is the level value for the contour line #linenum
    v = a(2,linenum);    % v is the number of vertices at the level lev(k)
    ptsmat(k) = v;
    for i = 1:v-1
        len(k) = len(k) + ptdist(a,linenum+i,linenum+i+1);  % Increment the length of the contour line by the distance between the two vertices
    end
    linenum = linenum+v+1;  % Move on to the next contour line
end

clen = len(1:k)';
clev = lev(1:k)';


% The following function ptdist calculates the distance between two points
% in cartesian coordinates

function l = ptdist(a,i1,i2)

dx = a(1,i1) - a(1,i2);
dy = a(2,i1) - a(2,i2);
l = sqrt((dx)^2 + (dy)^2);