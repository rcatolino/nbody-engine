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
dtmin = 50; %(seconds)

% time increments : (in seconds)
dt = 600; %(in seconds)
kmax = floor((24*3600*l)/dt)

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
        0, 0.9*orbitalspeed(200*10^9, BM(1), BM(2)), 0; % Second star speed.
        vlx, vly, 0]; % Stationary planet.
        %0, 3978, 0; % Earth speed.
        %-13007, 0, 0]; % Second planet speed (jupiter).
        %0, 29780, 0]; % Earth speed.

% Compute the initial positions at t = t0 + dt (k=1).
BP1 = ics(BP0, BV0, dt);

% Position vectors and initialisation :
X = zeros(kmax,n);
Y = zeros(kmax,n);
Z = zeros(kmax,n);
err = zeros(kmax,n);

X(1,1:n) = BP0(1:n,1)';
Y(1,1:n) = BP0(1:n,2)';
Z(1,1:n) = BP0(1:n,3)';
X(2,1:n) = BP1(1:n,1)';
Y(2,1:n) = BP1(1:n,2)';
Z(2,1:n) = BP1(1:n,3)';
% Simulation :

for k = 3:kmax
  if mod(k,1000) == 0
    k
  end
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
    X(k,j) = finite_diff(Sx, X(k-2,j), X(k-1,j), dt);
    Y(k,j) = finite_diff(Sy, Y(k-2,j), Y(k-1,j), dt);
    Z(k,j) = finite_diff(Sz, Z(k-2,j), Z(k-1,j), dt);

    % Store the result of f for each dimension.
    S(j,:) = [Sx; Sy; Sz];
  end
  % Second pass for heun's method
  for j = 1:n
    Shx = 0; Shy = 0; Shz = 0;
    for i = 1:n
      if i ~= j
        Ci = G*BM(i);
        [Fx, Fy, Fz] = f(X(k, i), Y(k, i), Z(k, i), X(k, j), Y(k, j), Z(k, j));
        Shx = Shx + Ci*Fx;
        Shy = Shy + Ci*Fy;
        Shz = Shz + Ci*Fz;
      end
    end
    X(k,j) = finite_diff((Shx+S(j,1))/2, X(k-2,j), X(k-1,j), dt);
    Y(k,j) = finite_diff((Shy+S(j,2))/2, Y(k-2,j), Y(k-1,j), dt);
    Z(k,j) = finite_diff((Shz+S(j,3))/2, Z(k-2,j), Z(k-1,j), dt);
    err(k,j) = norm([X(k,j)-X(k,j), Y(k,j)-Y(k,j), Z(k,j)-Z(k,j)])/norm([X(k,j), Y(k,j), Z(k,j)]);
    if k > 1000 && err(k,j) > 4*10^-8
      plotngrid(X(1:k,:),Y(1:k,:), BM)
      figure;
      hold on
      for j = 1:n
        plot(err(1000:k,j));
      end
      hold off
      return
    end
  end
end

%plotngrid(X, Y, BM)
plot(X,Y);
figure
[X1,Y1] = mvmasscenter(X, Y, BM);
plot(X1, Y1);
figure

hold on
for j = 1:n
  plot(err(1000:kmax,j));
end
hold off
