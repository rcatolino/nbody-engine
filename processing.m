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
% number of step (length of the simulation in days)
l = 1600;

% time increments : (in seconds)
dt = 3600; %(in seconds)
kmax = (24*3600*l)/dt;

% Bodies description
% Number of bodies
n = 3;
% Initial positions and speed (respectively in meters and meters/second) (at t = t0 (k = 0))
BP0 = [ 0,        0, 0; % Central body (sun)
        %0, 75*10^9, 0; % Second orbiting body (other planet)
        150*10^9, 0, 0; % Orbiting body (earth)
        0, 800*10^9, 0]; % Second orbiting body (jupiter)
        %150*10^9, 0, 0]; % Orbiting body (earth)
BV0 = [ 0, 0, 0; % Initialy immobile sun
        %-30000, 0, 0; % Second planet speed.
        0, 29780, 0; % Earth speed.
        -13007, 0, 0]; % Second planet speed (jupiter).
        %0, , 0]; % Earth speed.
% Masses (in kg);
BM = [ 1.988435*10^30, 5.9721986*10^24, 1.9*10^27];

% Compute the initial positions at t = t0 + dt (k=1).
BP1 = ics(BP0, BV0, dt);

% Position vectors and initialisation :
X = zeros(kmax,n);
Y = zeros(kmax,n);
Z = zeros(kmax,n);

X(1,1:n) = BP0(1:n,1)';
Y(1,1:n) = BP0(1:n,2)';
Z(1,1:n) = BP0(1:n,3)';
X(2,1:n) = BP1(1:n,1)';
Y(2,1:n) = BP1(1:n,2)';
Z(2,1:n) = BP1(1:n,3)';

% Simulation :

for k = 3:kmax
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
  end
end

plot(X,Y)
grid on
