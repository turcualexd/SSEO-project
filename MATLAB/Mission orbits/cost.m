function dV = cost(x)


T_dep = x(1);
T_dsm = x(2);
T_fb  = x(3);
T_arr = x(4);

a = x(5); e = x(6); i = x(7); OM = x(8); om = x(9); th = x(10);

    % [dv1, V_m] = cost_lambert(T_dep, T_fb, 3, 3, 0);
    % [dv2, V_p] = cost_lambert_2(T_fb, T_arr, 3, 5, 0);
    % [dv3, ~]      = cost_gravity_assist(T_fb, 3, V_m, V_p);
    % dV = dv1 + dv2 + dv3;

    [r_dsm, ~]                  = kep2car([a, e, i, OM, om, th], astroConstants(4));
    [dV_1, dV_dsm, dV_fb, dV_2] = cost_lambert_3(T_dep, T_dsm, T_fb, T_arr, r_dsm);

    dV = dV_1 + dV_dsm + dV_fb + dV_2;


end