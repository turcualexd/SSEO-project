function [R av ev iv OMv omv thv]  = loadData(filename)

% takes as an input a file name from nasa horizon with the following modes:
% OOE, ICRF, ecliptic x-y plane, gregorian, km/s, CSV ON

T = readtable(filename,MissingRule="omitrow");
clc
ev  = table2array(T(:,3));
rpv = table2array(T(:,4));
iv  = wrapTo2Pi(deg2rad(table2array(T(:,5))));
OMv = wrapTo2Pi(deg2rad(table2array(T(:,6))));
omv = wrapTo2Pi(deg2rad(table2array(T(:,7))));
av  = table2array(T(:,12));
thv = wrapTo2Pi(deg2rad(table2array(T(:,11))));
mu  = astroConstants(4);

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