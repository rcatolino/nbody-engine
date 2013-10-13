function [ X1, Y1 ] = mvmasscenter( X, Y, BM )
% Recenter the trajectories in the X,Y matrices upon the barycenter of
% all the bodies. Each column of X and Y are interpreted as the set
% of coordinates describing the trajectory of a body i
n = size(X,2);
if n ~= size(Y,2)
  return
end
X1 = X;
Y1 = Y;
Cx = zeros(size(X,1),1);
Cy = zeros(size(X,1),1);
for i = 1:n
  Cx = Cx + X(:,i)*BM(i);
  Cy = Cy + Y(:,i)*BM(i);
end
Cx = Cx/sum(BM);
Cy = Cy/sum(BM);
for i = 1:n
  X1(:,i) = X(:,i) - Cx(:);
  Y1(:,i) = Y(:,i) - Cy(:);
end
