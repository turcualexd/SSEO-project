clear, clc, close all;

name = "Science\Juno_Sun_science";

T = readtable("Original\" + name + ".txt");
%return
T = T(3:20347,:);

av  = table2array(T(:,12));
ev  = table2array(T(:,3));
iv  = wrapTo2Pi(deg2rad(table2array(T(:,5))));
OMv = wrapTo2Pi(deg2rad(table2array(T(:,6))));
omv = wrapTo2Pi(deg2rad(table2array(T(:,7))));
thv = wrapTo2Pi(deg2rad(table2array(T(:,11))));

muS = astroConstants(4);
muE = astroConstants(13);
muJ = astroConstants(15);

rv_JS_science = nan(length(av),3);
for i = 1:length(av)
    [rv_JS_science(i,:), ~] = kep2car([av(i) ev(i) iv(i) OMv(i) omv(i) thv(i)], muS);
end

%return
save("Radius\" + name + ".mat", "rv_JS_science")