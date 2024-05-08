clear; close all; clc;
T = readtable("horizons_results.txt",MissingRule="omitrow");
%
clc
ev  = table2array(T(:,3));
rpv = table2array(T(:,4));
iv  = wrapTo2Pi(deg2rad(table2array(T(:,5))));
OMv = wrapTo2Pi(deg2rad(table2array(T(:,6))));
omv = wrapTo2Pi(deg2rad(table2array(T(:,7))));
av  = table2array(T(:,12));
thv = wrapTo2Pi(deg2rad(table2array(T(:,11))));
mu  = astroConstants(4);
%%
R = nan(3,length(ev));

for k = 1 : length(ev)

    e  = ev(k);
    i  = iv(k);
    OM = OMv(k);
    om = omv(k);
    a = av(k);
    th = thv(k);
    vect = [a e i OM om th];
    [r, v] = kep2car(vect, mu);
    R(:, k) = r;

end
rn = vecnorm(R, 2);
%%
figure
hold on
axis equal
grid minor
plot3(R(1, :),R(2,:),R(3,:))
%%
figure
plot(rn/astroConstants(2))
%%
figure
hold on
plot(R(1,:))
plot(R(2,:))
plot(R(3,:))