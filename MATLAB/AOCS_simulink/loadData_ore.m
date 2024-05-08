function [R av ev iv OMv omv thv]  = loadData_ore(filename)

% takes as an input a file name from nasa horizon with the following modes:
% OOE, ICRF, ecliptic x-y plane, gregorian, km/s, CSV ON

T = readtable(filename,MissingRule="omitrow");
clc
ev  = table2array(T(2:end,5));
rpv = table2array(T(2:end,6));
iv  = wrapTo2Pi(deg2rad(table2array(T(2:end,5+2))));
OMv = wrapTo2Pi(deg2rad(table2array(T(2:end,6+2))));
omv = wrapTo2Pi(deg2rad(table2array(T(2:end,7+2))));
av  = table2array(T(2:end,12+2));
thv = wrapTo2Pi(deg2rad(table2array(T(2:end,11+2))));
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