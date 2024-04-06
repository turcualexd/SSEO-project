function rp = compute_rp(mu, v_minus, v_plus, delta, R_plan)
%
% This function calculates the common radius of pericentre by solving a
% non-linear equation using fsolve.m. It takes as input the planetary
% constant and the radius of celestial body around which the fly-by is
% performed. The other 3 inputs are: the incoming and outgoing 
% excess hyperbolic velocities and the turn angle.
% ------------------------------------------------------------------------
%
% INPUT:
% mu             [1]:      in km^3/s^2;
% v_minus        [1]:      in km/s;
% v_plus         [1]:      in km/s;
% delta          [1]:      in rad;
% R_plan         [1]:      in km
% 
% OUTPUT:
% rp             [1]: in km
%
% Author: Marcello Pareschi.
%
%--------------------------------------------------------------------------

e_m = @(rp) 1 + rp .* v_minus .^ 2 ./ mu;

delta_m = @(rp) 2 * asin(1 ./ e_m(rp));

e_p = @(rp) 1 + rp .* v_plus .^ 2 ./ mu;

delta_p = @(rp) 2 * asin(1 ./ e_p(rp));

Delta = @(rp) delta_m(rp)/2 + delta_p(rp)/2-delta;

fun = @(rp) Delta(rp) ;

options.Display = 'off';

rp = fsolve(fun, R_plan, options);

return