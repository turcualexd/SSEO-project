function [r, v] = kep2car(kep, mu)

% kep2car.m - Calculate position and velocity vectors from Keplerian
%             elements
%
% PROTOTYPE:
%   [r, v] = kep2car(kep, mu)
%
% DESCRIPTION: 
%   Calculates the position and velocity vectors in the ECI reference frame
%   given the Keplerian elements and planetary constant of the main body
%
% INPUTS: 
%   Kep:        Vector containing the Keplerian elements in the form: 
%               [a e i Om om theta] [Km, rad]
%
%   mu:         Planetary constant of the main body [km^3/s^2]
%
% OUTPUTS:
%   r:          Position vector in the ECI reference frame [Km]
%
%   v:          Velocity vector in the ECI reference frame [Km/s]
%
% -------------------------------------------------------------------------

% Recover single Keplerian elements
a = kep(1);
e = kep(2);
i = kep(3);
OM = kep(4);
om = kep(5);
theta = kep(6);

% Calculate position and velocity vectors in perifocal reference frame
p = a * (1 - e^2);
r_pf = [p*cos(theta)/(1+e*cos(theta)); p*sin(theta)/(1+e*cos(theta)); 0];
v_pf = sqrt(mu/p) * [-sin(theta); (e + cos(theta)); 0];

% Define rotation matrix
R_i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];
R_OM = [cos(OM) sin(OM) 0; -sin(OM) cos(OM) 0; 0 0 1];
R_om = [cos(om) sin(om) 0; -sin(om) cos(om) 0; 0 0 1];
R = (R_om*R_i*R_OM)';

% Calculate position and velocity vectors in ECI reference frame
r = R * r_pf;
v = R * v_pf;