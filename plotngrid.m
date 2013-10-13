function [] = plotngrid( X, Y, BM )
hold on
[X1, Y1] = mvmasscenter(X, Y, BM);
[Xa, Ya] = mvaxis(X1, Y1);
plot(Xa(:,1),Ya(:,1),'xr')
plot(Xa(:,2),Ya(:,2),'xg')
plot(Xa(:,3),Ya(:,3),'b')
grid on
axis equal
hold off
end
