function [ X1, Y1 ] = mvaxis( X, Y )
% Compute the coordinates of the trajectories in (X,Y) in the new grid
% having the same center and an x axis the normalized vector from the
% first body to the second
Cosg = X(:,1);
Sing = Y(:,1);
Cgn = sqrt(Cosg.^2 + Sing.^2);
Cosg = Cosg./Cgn;
Sing = Sing./Cgn;

X1 = Cosg;
Y1 = Sing;

n = size(X,2);
if n ~= size(Y,2)
  return
end
X1 = X;
Y1 = Y;

for i = 1:n
  X1(:,i) = X(:,i).*Cosg + Y(:,i).*Sing;
  Y1(:,i) = X(:,i).*Sing - Y(:,i).*Cosg;
end

end
