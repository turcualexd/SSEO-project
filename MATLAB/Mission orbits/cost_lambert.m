function [deltav_tot, VF] = cost_lambert(departure, arrival, ibody1, ibody2, orbitType )
%
% This function calculates the deltaV of the injection into the 1st
% heliocentric leg between 2 planets. The IDs of the planets are referred
% to the uplanet.m function. The kep2car.m and lambertMR.m are used.
% Also VI and VF respectively the initial and final velocities on Lambert's
% arc, are calculated.
% ------------------------------------------------------------------------
%
% INPUT:
% departure      [1]:  in mjd2000;
% arrival        [1]:  in mjd2000;
% ibody1         [1]:  -;
% ibody2         [1]:  -;
% orbitType      [1]:  -;
%
% OUTPUT:
% deltav_tot     [1]: in km/s;
% VI             [1]: in km/s;
% VF             [1]: in km/s;
%
% Author: Marcello Pareschi.
%--------------------------------------------------------------------------

TOF = (arrival - departure)*24*3600;

[kep1,ksun1] = uplanet(departure, ibody1);

[RI, V1] = kep2car(kep1, ksun1);

[kep2,ksun2] = uplanet(arrival, ibody2);

[RF, ~] = kep2car(kep2, ksun2);

[~,~,~,~,VI,VF,~,~] = lambertMR(RI,RF,TOF,ksun1,orbitType,1,0);

deltav_tot =  norm(VI'-V1);

VI = VI';

VF = VF';