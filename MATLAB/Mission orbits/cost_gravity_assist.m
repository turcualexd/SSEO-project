function [dv,rp] = cost_gravity_assist(time, id, V_m, V_p)
%
% This function calculates the deltaV at the pericentre necessary to perform 
% the requested fly by. Also the common pericentre radius is computed. The
% function uses the date and planet of fly-by, the two heliocentric
% velocities before and after the fly-by. The id of the planet is referred
% to uplanet.m, also astroConstants.m and kep2car.m are used.
% ------------------------------------------------------------------------
%
% INPUT:
% time [fly-by]  [1]:  in mjd2000;
% id             [1]:  -;
% V_m            [3x1]:  km/s;
% V_p            [3x1]:  km/s;
% 
% OUTPUT:
% dv             [1]: in km/s;
% rp             [1]: in km
%
% Author: Marcello Pareschi.
%--------------------------------------------------------------------------

mu_E = astroConstants(13);

RE = astroConstants(23);

[kep1,ksun1] = uplanet(time, id);

[~, VP] = kep2car(kep1, ksun1);

v_m = V_m - VP;

v_p = V_p - VP;

delta = acos(dot(v_m, v_p) / (norm(v_m) * norm(v_p)));

rp = compute_rp(mu_E, norm(v_m), norm(v_p), delta, RE);

v_mn = norm(v_m);

v_pn = norm(v_p);

vp1 = sqrt(v_mn^2 + 2*mu_E/rp);

vp2 = sqrt(v_pn^2 + 2*mu_E/rp);

dv = norm(vp1 - vp2);

    if rp < (astroConstants(23) + 500)

        dv = 100;

    end

end