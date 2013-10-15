function [] = plotngrid( X, Y, BM )
n = size(X,2);
if n < 3
  return;
end
hold on
[X1, Y1] = mvmasscenter(X, Y, BM);
[Xa, Ya] = mvaxis(X1, Y1);
plot(Xa(:,1),Ya(:,1),'xr')
plot(Xa(:,2),Ya(:,2),'xg')
plot(Xa(:,3:n),Ya(:,3:n))
grid on
axis equal
hold off
end
