clear, clc, close all;

T = readtable("Juno_Earth_EGA_original.txt");
T = T(3:4323,:);

av  = table2array(T(:,12));
ev  = table2array(T(:,3));
iv  = wrapTo2Pi(deg2rad(table2array(T(:,5))));
OMv = wrapTo2Pi(deg2rad(table2array(T(:,6))));
omv = wrapTo2Pi(deg2rad(table2array(T(:,7))));
thv = wrapTo2Pi(deg2rad(table2array(T(:,11))));

muS = astroConstants(4);
muE = astroConstants(13);

rv_JE_EGA = nan(length(av),3);
for i = 1:length(av)
    [rv_JE_EGA(i,:), ~] = kep2car([av(i) ev(i) iv(i) OMv(i) omv(i) thv(i)], muE);
end

clear T i muS muE

save("Juno_Earth_EGA.mat", "rv_JE_EGA")