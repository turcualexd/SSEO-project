function kep = car2kep(rr,vv,mu)

%------------------------Cartesian to Keplerian----------------------------
%
% This function returns the Keplerian elements of an orbiting object given
% its position and velocity vectors and the gravitational constant of the
% attrator.
%
%--------------------------------------------------------------------------
%
% INPUT:
% rr              [3x1]: position vector in cartesian frame, in km
% vv              [3x1]: velocity vector in cartesian frame, in km/s
% mu                [1]: gravitational constant [km^3/s^2]
%
% OUTPUT:
% kep             [6x1]: Keplerian elements of the orbit
%                        (a, e, i, OM, om, theta)
%
%--------------------------------------------------------------------------
%
% AUTHOR: Alex Cristian Turcu
%
%--------------------------------------------------------------------------

if size(rr,1) == 1
    rr = rr.';
end

if size(vv,1) == 1
    vv = vv.';
end

r = norm(rr);
v = norm(vv);
k = [0 0 1]';

a = -mu / (v^2 - 2*mu/r);

hh = cross(rr,vv);
h = norm(hh);

ee = cross(vv,hh)/mu - rr/r;
e = norm(ee);
if e == 0
    ee = [1; 0; 0];
end

i = acos(hh(3)/h);
N = cross(k,hh) / norm(cross(hh,k));

if i == 0
    i = eps;
    N = [1; 0; 0];
end

Om = acos(N(1));
if N(2) < 0
    Om = 2*pi - Om;
end

om = acos(dot(N,ee)/e);
if ee(3) < 0
    om = 2*pi - om;
end

theta = acos(dot(rr,ee)/(r*e));
if dot(rr,vv) < 0
    theta = 2*pi - theta;
end

kep = [a; e; i; Om; om; theta];