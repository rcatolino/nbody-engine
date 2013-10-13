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
l = 360;

% time increments : (in seconds)
dt = 600; %(in seconds)
kmax = (24*3600*l)/dt

% Bodies description
% Number of bodies
n = 2;
% Initial positions and speed (respectively in meters and meters/second) (at t = t0 (k = 0))
BP0 = [ 0,        0, 0; % Central body (sun)
        %0, 75*10^9, 0; % Second orbiting body (other planet)
        150*10^9, 0, 0; % Orbiting body (earth)
        0, 800*10^9, 0]; % Second orbiting body (jupiter)
        %150*10^9, 0, 0]; % Orbiting body (earth)
BV0 = [ 0, 0, 0; % Initialy immobile sun
        %-30000, 0, 0; % Second planet speed.
        0, 3978, 0; % Earth speed.
        -13007, 0, 0]; % Second planet speed (jupiter).
        %0, 29780, 0]; % Earth speed.
% Masses (in kg);
BM = [ 1.988435*10^30, 5.9721986*10^24, 1.9*10^27];

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
Xh = X;
Yh = Y;
Zh = Z;
% Simulation :

for k = 3:kmax
  S = zeros(n,3);
  for j = 1:n
    Sx = 0; Sy = 0; Sz = 0;
    for i = 1:n
      if i ~= j
        Ci = G*BM(i);
        [Fx, Fy, Fz] = f(Xh(k-1, i), Yh(k-1, i), Zh(k-1, i), Xh(k-1, j), Yh(k-1, j), Zh(k-1, j));
        Sx = Sx + Ci*Fx;
        Sy = Sy + Ci*Fy;
        Sz = Sz + Ci*Fz;
      end
    end
    X(k,j) = finite_diff(Sx, Xh(k-2,j), Xh(k-1,j), dt);
    Y(k,j) = finite_diff(Sy, Yh(k-2,j), Yh(k-1,j), dt);
    Z(k,j) = finite_diff(Sz, Zh(k-2,j), Zh(k-1,j), dt);

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
    Xh(k,j) = finite_diff((Shx+S(j,1))/2, Xh(k-2,j), Xh(k-1,j), dt);
    Yh(k,j) = finite_diff((Shy+S(j,2))/2, Yh(k-2,j), Yh(k-1,j), dt);
    Zh(k,j) = finite_diff((Shz+S(j,3))/2, Zh(k-2,j), Zh(k-1,j), dt);
    err(k,j) = norm([Xh(k,j)-X(k,j), Yh(k,j)-Y(k,j), Zh(k,j)-Z(k,j)])/norm([Xh(k,j), Yh(k,j), Zh(k,j)]);
    if err(k,j) > 0.02
      hold on
      plot(Xh(1:k,:),Yh(1:k,:),'b')
      plot(X(1:k,:),Y(1:k,:),'r')
      hold off
      figure
      plot(err(1:k,j));
      return
    end
  end
end

hold on
plot(X,Y,'r')
plot(Xh,Yh,'b')
hold off
grid on
figure

hold on
for j = 1:n
  plot(err(1000:kmax,j));
end
hold off
