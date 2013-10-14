function [ x, y, vx, vy ] = l4char( R, m1, m2)
% R is the distance between the two massive bodies, m1 and m2 their respective masses.
% B1 is considered to sit at 0,0 and B2 at R,0.
if m2 > m1
  return;
end
% L4 is at the 3rd vertex of the equilateral triangle of base B1-B2 :
x = R/2;
y = sqrt(3/4)*R;
r = m1/m2;
alpha = atan((r-2)/(r*sqrt(3)))
speed = orbitalspeed(R, m1, m2);
vy = speed*sin(alpha);
vx = -speed*cos(alpha);
end
