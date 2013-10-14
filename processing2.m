clear all;
close all;
clc;

% Compute the positions of n body interacting with each other through gravitational
% forces, given initial positions and speed, using the second order finite difference
% approximation.

% Gravitational constant :
G = 6.67300*10^-11; % Unit : m³.kg¯¹.s¯²

% intial time : (in seconds)
t0 = 0;
% length of the simulation (in days)
l = 404;
errmax = 10*10^-9;
dtmin = 5; %(seconds)
lend = l*3600*24;
lcurrent = 0;

% time increments : (in seconds)
dt = 3600; %(in seconds)

% Bodies description
% Number of bodies
n = 2;
% Masses (in kg);
%BM = [ 1.988435*10^30, 5.9721986*10^24, 1.9*10^27];
BM = [ 30*1.988435*10^30, 1.988435*10^30, 1.9721986*10^24];
[lx, ly, vlx, vly] = l4char(200*10^9, BM(1), BM(2));
% Initial positions and speed (respectively in meters and meters/second) (at t = t0 (k = 0))
BP0 = [ 0,        0, 0; % Central body (big sun)
        200*10^9, 0, 0; % Second star.
        lx, ly, 0]; % Orbiting body (planet)
        %150*10^9, 0, 0; % Orbiting body (earth)
        %0, 800*10^9, 0]; % Second orbiting body (jupiter)
BV0 = [ 0, 0, 0; % Initialy immobile sun
        0, orbitalspeed(200*10^9, BM(1), BM(2)), 0; % Second star speed.
        vlx, vly, 0]; % Stationary planet.
        %0, 3978, 0; % Earth speed.
        %-13007, 0, 0]; % Second planet speed (jupiter).
        %0, 29780, 0]; % Earth speed.

% Compute the initial positions at t = t0 + dt (k=1).
BP1 = ics(BP0, BV0, dt);

% Position vectors and initialisation :
X = zeros(2,n);
Y = zeros(2,n);
Z = zeros(2,n);
err = zeros(2,n);
k = 2;

X(1,1:n) = BP0(1:n,1)';
Y(1,1:n) = BP0(1:n,2)';
Z(1,1:n) = BP0(1:n,3)';
X(2,1:n) = BP1(1:n,1)';
Y(2,1:n) = BP1(1:n,2)';
Z(2,1:n) = BP1(1:n,3)';
X = X;
Y = Y;
Z = Z;
% Simulation :
p = 0;
while k < 3
  k = k+1;
  % display percentage
  if floor(100*lcurrent/lend) == p
    p
    p = p+1;
  end
  lcurrent = lcurrent + dt;
  X = [X; zeros(1,n)];
  Y = [Y; zeros(1,n)];
  Z = [Z; zeros(1,n)];
  S = zeros(n,3);
  for j = 1:n
    Sx = 0; Sy = 0; Sz = 0;
    for i = 1:n
      if i ~= j
        Ci = G*BM(i);
        [Fx, Fy, Fz] = f(X(k-1, i), Y(k-1, i), Z(k-1, i), X(k-1, j), Y(k-1, j), Z(k-1, j));
        Sx = Sx + Ci*Fx;
        Sy = Sy + Ci*Fy;
        Sz = Sz + Ci*Fz;
      end
    end
    Xkj = finite_diff(Sx, X(k-2,j), X(k-1,j), dt)
    Ykj = finite_diff(Sy, Y(k-2,j), Y(k-1,j), dt)
    Zkj = finite_diff(Sz, Z(k-2,j), Z(k-1,j), dt)
    X(k,j) = Xkj;
    Y(k,j) = Ykj;
    Z(k,j) = Zkj;

    % Store the result of f for each dimension.
    S(j,:) = [Sx; Sy; Sz];
  end

  % Second pass for heun's method
  for j = 1:n
    Sx = 0; Sy = 0; Sz = 0;
    for i = 1:n
      if i ~= j
        Ci = G*BM(i);
        [Fx, Fy, Fz] = f(X(k, i), Y(k, i), Z(k, i), X(k, j), Y(k, j), Z(k, j));
        Sx = Sx + Ci*Fx;
        Sy = Sy + Ci*Fy;
        Sz = Sz + Ci*Fz;
      end
    end
    Xkj = finite_diff((Sx+S(j,1))/2, X(k-2,j), X(k-1,j), dt)
    Ykj = finite_diff((Sy+S(j,2))/2, Y(k-2,j), Y(k-1,j), dt)
    Zkj = finite_diff((Sz+S(j,3))/2, Z(k-2,j), Z(k-1,j), dt)
    err(k,j) = norm([X(k,j)-Xkj, Y(k,j)-Ykj, Z(k,j)-Zkj])/norm([Xkj, Ykj, Zkj])
    X(k,j) = Xkj;
    Y(k,j) = Ykj;
    Z(k,j) = Zkj;
  end

  % Check the error for each j
  if k > 1000 && max(err(k,:)) > 10^-9
    %dt = dt*0.8
    if dt < dtmin
      plotngrid(X(1:k,:),Y(1:k,:), BM)
      figure;
      hold on
      for j = 1:n
        plot(err(1000:k,j));
      end
      hold off
      error('error too big');
    end
  elseif k > 1000 && max(err(k,:)) < 10^-10
    %dt = dt*1.1
  end
end

plotngrid(X, Y, BM)
figure
[X1,Y1] = mvmasscenter(X, Y, BM);
plot(X1, Y1);
figure

plot(err(1000:k,:));
