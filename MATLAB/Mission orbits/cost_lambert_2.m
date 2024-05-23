function [deltav_tot, VI] = cost_lambert_2(departure, arrival, ibody1, ibody2)
%
% This function calculates the deltaV of the exit from the 2nd
% heliocentric leg between a planet and a NEO object. The IDs  are referred
% to the uplanet.m function and ephNEO.m. The kep2car.m and lambertMR.m are used.
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

[RI, ~] = kep2car(kep1, ksun1);

[kep2,~] = uplanet(arrival, ibody2);

[RF, V2] = kep2car(kep2, ksun1);

[~,~,~,~,VI,VF,~,~] = lambertMR(RI,RF,TOF,ksun1);

deltav_tot = norm(VF'-V2);

VI = VI';

VF = VF';