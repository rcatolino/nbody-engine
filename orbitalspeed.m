function [ vo ] = orbitalspeed( r, m1, m2 )
% m1 is the mass of the orbited body, and m2 the mass of the orbiting body
G = 6.67300*10^-11; % Unit : m³.kg¯¹.s¯²
vo = sqrt((G*m1^2)/((m1+m2)*r));%*1.0340;
end
