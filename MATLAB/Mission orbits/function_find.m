function [dv_min, dv_1, dv_2, dv_3, dep_date, fb_date, arr_date] = function_find(T_dep, T_arr, T_fb)

%------------------------Brute-force algorithm-----------------------------
%
% This function calculates the minimum deltaV for a interplanetary
% trajectory from Mars to 1036 Ganymed with a fly-by on Earth. The function
% uses three time-based DOF: departure date range, time of flight of
% first leg and time of flight of second leg. The output is the calculated
% minimum total cost, computed as the injection cost on first Lambert
% arc, assist cost at pericentre of flyby hyperbola and exit cost from second
% Lambert arc. No injection or exit hyperbolae are considered
%
%--------------------------------------------------------------------------
% 
% INPUT:
% T_dep        [1xN]:  in mjd2000;
% TOF_1        [1xM]:  in mjd2000;
% TOF_2        [1xK]:  in mjd2000;
%
% OUTPUT:
% dv_min         [1]: in km/s;
% dep_date       [1]: in mjd2000;
% fb_date        [1]: in mjd2000;
% arr_date       [1]: in mjd2000;
%
%--------------------------------------------------------------------------
%
% AUTHORs: Marcello Pareschi; Daniele Paternoster
%
%--------------------------------------------------------------------------


dv_min = 1e4;
wait = waitbar(0, 'Minimum search... (0%)');
for i = 1 : length(T_dep)
   
    waitbar(i/length(T_dep), wait, sprintf('Minimum search... (%g%%)', i/length(T_dep)* 100));
   
   for j = 1 : length(T_fb)

        [dv1, ~, V_m] = cost_lambert(T_dep(i), T_fb(j), 3, 3, 0);
        
        for k = 1 : length(T_arr)
            
            [dv2, V_p, ~] = cost_lambert_2(T_fb(j), T_arr(k), 3, 5, 0);
            [dv3, rp]     = cost_gravity_assist(T_fb(j), 3, V_m, V_p);
            
            dvtot = dv1 + dv2 + dv3;

            if (dvtot < dv_min) && rp > (astroConstants(23) + 500)

                dv_min = dvtot;
                dv_1    = dv1;
                dv_2    = dv2;
                dv_3    = dv3;
                dep_date = T_dep(i);
                fb_date =  T_fb(j);
                arr_date = T_arr(k);

            end
        end
    end
end

close(wait);
