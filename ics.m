function [ P1 ] = ics(P0, V0, dt)
% Compute the intial conditions for the bodies, and the time k = 1 : x(j,1), y(j,1), z(j,1)
% from the initial position and speed of the body j.

P1 = P0 + (dt*V0);
end
