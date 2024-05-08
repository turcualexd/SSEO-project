function [V1, V2, V3] = dsm_maneuver(a,e,i,OM,om,th,T_dep,T_dsm)

    [R,~] = kep2car([a,e,i,OM,om,th],astroConstants(4));
    [kep1,ksun1] = uplanet(T_dep, 3);
    [RI, V1] = kep2car(kep1, ksun1);
    [~,~,~,~,V2,V3,~,~] = lambertMR(RI,R,(T_dep - T_dsm)*24*3600, ksun1, 0);


end