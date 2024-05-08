clear, clc, close all

%% Input values

g0 = 9.81;              % m/s^2
M_dry = 1593;           % kg
Isp_me = 317;           % s
Isp_rcs = 215;          % s
T = 662;                % N

%% Import delta-V for main engine and RCS

manovre = load("Manovre.txt");
dv_vect = manovre(:, 1);
type_vect = manovre(:, 2);

dv_me = 0;
dv_rcs = 0;
for i = 1 : length(dv_vect)
    if type_vect(i)
        dv_me = dv_me + dv_vect(i);
    else
        dv_rcs = dv_rcs + dv_vect(i);
    end
end

%% Compute total mass for fuel and oxidizer

M_rcs = 1.02 * 1.2 * M_dry * ( exp((1.05*dv_rcs) / (Isp_rcs*g0)) - 1 );
M_me = 1.02 * (1.2 * M_dry + M_rcs) * ( exp((1.05*dv_me) / (Isp_me*g0)) - 1 );
M = 1.2*M_dry + M_me + M_rcs;

%% Compute burn time for me maneuvers

t_burn = [];
for i = 1 : length(dv_vect)
    if type_vect(i)
        M_p = M*(1 - exp(-dv_vect(i)/(Isp_me*g0)));
        M = M - M_p;
        t_burn = [t_burn M_p*Isp_me*g0/T];
    else
        M_p = M*(1 - exp(-dv_vect(i)/(Isp_rcs*g0)));
        M = M - M_p;
    end
end

% Duration in minutes of DSM-1, DSM-2, JOI and PRM respectively 
t_burn = t_burn([1 2 6 7])/60