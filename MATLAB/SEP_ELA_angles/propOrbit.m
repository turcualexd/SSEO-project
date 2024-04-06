function [r, v] = propOrbit(r0, v0, mu, t)
y0 = [r0; v0];
options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[~, y] = ode113(@(t, y) ode_2bp(t, y,mu), t, y0, options);
y = y';
r = y(1:3, :);
v = y(4:6, :);