clear, clc, close all

el = load("Elementi Juno.txt");
el = [el(:, 2:3) deg2rad(el(:, 4:7))];
dsm1 = 390;
dsm2 = 406;
el_dsm1 = el(dsm1, :);
el_dsm2 = el(dsm2, :);
el_dsm12 = el(dsm1:dsm2, :);
mu = astroConstants(4);
r_dsm12 = zeros(3,17);
for i = 1:17
    r_dsm12(:,i)  = kep2car(el_dsm12(i,:),mu);
end
    r_dsm1 = kep2car(el_dsm1, mu);
    r_dsm2 = kep2car(el_dsm2, mu);
r_e_dsm1 = [1.388828997654905E+08 -5.936533358921350E+07 1.159179357510060E+03];
r_e_dsm2 = [1.487946423457917E+08 -2.248191275651955E+07 8.331244925418869E+02];
drawPlanet(10, [0 0 0], 15, 1)
title("DSM-1 and DSM-2")
drawPlanet(3, r_e_dsm1, 750)
plot3(r_dsm1(1), r_dsm1(2), r_dsm1(3), "r.", "MarkerSize", 14)
lw = 0.25;
plot3([r_e_dsm1(1) 0], [r_e_dsm1(2) 0], [r_e_dsm1(3) 0], 'b--', 'LineWidth', lw)
plot3([r_e_dsm1(1) r_dsm1(1)], [r_e_dsm1(2) r_dsm1(2)], [r_e_dsm1(3) r_dsm1(3)], 'r--', 'LineWidth', lw)
% drawPlanet(10, [0 0 0], 15, 1)
% title("DSM-2")
drawPlanet(3, r_e_dsm2, 750)
plot3(r_dsm2(1), r_dsm2(2), r_dsm2(3), "k.", "MarkerSize", 14)
plot3([r_e_dsm2(1) 0], [r_e_dsm2(2) 0], [r_e_dsm2(3) 0], 'b--', 'LineWidth', lw)
plot3([r_e_dsm2(1) r_dsm2(1)], [r_e_dsm2(2) r_dsm2(2)], [r_e_dsm2(3) r_dsm2(3)], 'r--', 'LineWidth', lw)
legend('','','DSM1', '', '','','DSM2','','');
a1 = -r_e_dsm1';
b1 = r_dsm1 + a1;
sep1 = rad2deg(acos(dot(a1, b1)/(norm(a1)*norm(b1))))
a2 = -r_e_dsm2';
b2 = r_dsm2 + a2;
sep2 = rad2deg(acos(dot(a2, b2)/(norm(a2)*norm(b2))))