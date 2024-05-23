function [dV_1, dV_dsm, dV_fb, dV_2] = cost_lambert_3(T_dep, T_dsm, T_fb, T_arr, r_dsm)

    
    TOF1                         = (T_dsm - T_dep)*3600*24;
    TOF2                         = (T_fb - T_dsm)*3600*24;
    
    [kep1,ksun1]                 = uplanet(T_dep, 3);
    [RI, V1]                     = kep2car(kep1, ksun1);
    [~,~,~,~,VI,V_dsm_min,~,~]   = lambertMR(RI, r_dsm, TOF1, ksun1);
    
    [kep_fb,~]                   = uplanet(T_fb, 3);
    [R_fb, ~]                    = kep2car(kep_fb, ksun1);
    [~,~,~,~,V_dsm_plu,V_m,~,~]  = lambertMR(r_dsm, R_fb, TOF2, ksun1);
    
    [dV_2, V_p]                  = cost_lambert_2(T_fb, T_arr,3,5);
    [dV_fb, ~]                   = cost_gravity_assist(T_fb, 3, V_m', V_p);

    dV_1 = norm(VI' - V1);
    dV_dsm = norm(V_dsm_plu' - V_dsm_min');
    



end