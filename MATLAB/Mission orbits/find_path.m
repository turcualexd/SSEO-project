
%% Genetic Algorithm 2
clear; close all; clc;

dep_1 = [2011 08 05 00 00 00];
dep_2 = [2011 08 26 00 00 00];

 % dep_1 = [2011 08 05 17 20 00];
 % dep_2 = [2011 08 05 17 30 00];

dsm_1 = [2012 08 20 00 00 00];
dsm_2 = [2012 09 10 00 00 00];

fb_1  = [2013 09 20 00 00 00];
fb_2  = [2013 10 20 00 00 00];

arr_1 = [2016 06 20 00 00 00];
arr_2 = [2016 07 10 00 00 00];

lb = [date2mjd2000(dep_1) date2mjd2000(dsm_1) date2mjd2000(fb_1) ...
    date2mjd2000(arr_1) 2.4E+08 0.3740 0 deg2rad(30) deg2rad(180) deg2rad(170)];

ub = [date2mjd2000(dep_2) date2mjd2000(dsm_2) date2mjd2000(fb_2) ...
date2mjd2000(arr_2) 2.45E+08 0.4500 deg2rad(1) deg2rad(140) deg2rad(300) deg2rad(190)];

A = [];
b = [];

Aeq = [];
beq = [];

opts = optimoptions('ga','ConstraintTolerance',1e-6);
[x, fval] = ga(@cost,10,A,b,Aeq,beq,lb,ub,[],opts);




dep_gregGA = mjd20002date(x(1));
dsm_gregGA = mjd20002date(x(2));
fb_gregGA  = mjd20002date(x(3));
arr_gregGA = mjd20002date(x(4));
[r_dsm,~] = kep2car(x(5:end),astroConstants(4)); 
[dV_1, dV_dsm, dV_fb, dV_2] = cost_lambert_3(x(1), x(2), x(3), x(4),r_dsm);
%%

orbitType = 0;
Nrev = 0;
Ncase = 0;
optionsLMR = 0;

[kep_E1,ksun] = uplanet(x(1), 3);
[R_E1, ~] = kep2car(kep_E1, ksun);
[a1,p1,e1,~,VI,VF1,~,dth1] = lambertMR(R_E1,r_dsm,(x(2)-x(1))*3600*24,ksun,orbitType,Nrev,Ncase,optionsLMR);
kep_arc1 = car2kep(R_E1, VI, ksun);

[kep_E2,ksun] = uplanet(x(3), 3);
[R_E2, ~] = kep2car(kep_E2, ksun);
[a2,p2,e2,~,VI2,VF2,~,dth2] = lambertMR(r_dsm,R_E2,(x(3)-x(2))*3600*24,ksun,orbitType,Nrev,Ncase,optionsLMR);
kep_arc2 = car2kep(r_dsm, VI2, ksun);

[kep_J,ksun] = uplanet(x(4), 5);
[R_J, ~] = kep2car(kep_J, ksun);
[a3,p3,e3,~,VI3,VF3,~,dth3] = lambertMR(R_E2,R_J,(x(4)-x(3))*3600*24,ksun,orbitType,Nrev,Ncase,optionsLMR);
kep_arc3 = car2kep(R_E2, VI3, ksun);

dth = pi / 100;
Terra_3D(1e-4)
plotorbit(kep_arc1, kep_arc1(6), kep_arc1(6) + dth1, dth, ksun, 'b')
plotorbit(kep_arc2, kep_arc2(6), kep_arc2(6) + dth2, dth, ksun, 'r')
plotorbit(kep_arc3, kep_arc3(6), kep_arc3(6) + dth3, dth, ksun, 'g')

%%
dv_peri = dV_2*sqrt((1-0.9816)/2)