function [ fx, fy, fz ] = f( xi, yi, zi, xj, yj, zj)
% Compute the f part of the equation, corresponding to the distance term
% in the gravity equation between the bodies i and j.
abs2 = (xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2;
d = abs2 * sqrt(abs2);

fx = (xi - xj) / d;
fy = (yi - yj) / d;
fz = (zi - zj) / d;
end
