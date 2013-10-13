function [ X1, Y1 ] = mvcenter( X, Y )
% Recenter the trajectories in the X,Y matrices upon the body corresponding to
% the first column of X and Y. Each collumn of X and Y are interpreted as
% the set of coordinates describing the trajectory of a body i
n = size(X,2);
if n ~= size(Y,2)
  return
end
X1 = X;
Y1 = Y;

for i = 1:n
  X1(:,i) = X(:,i) - X(:,1);
  Y1(:,i) = Y(:,i) - Y(:,1);
end
