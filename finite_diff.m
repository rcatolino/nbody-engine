function [ x2 ] = finite_diff( s, x0, x1, dt)
% Compute the next discrete term on the time scale
% x(k+1) = S + 2x(k) - x(k-1)
x2 = (dt^2)*s + 2*x1 - x0;
end
