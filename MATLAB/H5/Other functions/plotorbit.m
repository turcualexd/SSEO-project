function plotorbit(kep, th0, thf, dth, mu, color)

%------------------------Plot orbit arc in 3D plot-------------------------
%
% This function plots an orbit giving its Keplerian elements, the starting
% and the ending theta, the resolution, the gravitational constant and the
% color.
%
%--------------------------------------------------------------------------
%
% INPUT:
% kep             [6x1]: Keplerian elements of the orbit
%                        (a, e, i, OM, om, theta)
% th0               [1]: initial value of true anomaly, in rad
% thf               [1]: final value of true anomaly, in rad
% dth               [1]: resolution of the plotted arc, in rad
% mu                [1]: gravitational constant [km^3/s^2]
% color           [str]: color for the plotted arc
%
%--------------------------------------------------------------------------
%
% AUTHOR: Alex Cristian Turcu
%
%--------------------------------------------------------------------------

if th0 <= thf
    theta = th0:dth:thf;
else
    theta = [th0:dth:2*pi,0:dth:thf];
end

R = zeros(length(theta), 3);

for t = 1:length(theta)
    th = theta(t);
    [rr, ~] = kep2car([kep(1:5); th], mu);
    R(t, :) = rr / astroConstants(2);
end

plot3(R(:, 1), R(:, 2), R(:, 3), 'Color', color, 'LineWidth', 1.5)