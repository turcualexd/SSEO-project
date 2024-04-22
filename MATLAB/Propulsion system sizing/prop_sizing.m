clear, clc, close all

%% Input values

g0 = 9.81;              % m/s^2
M_dry = 1593;           % kg
Isp_me = 317;           % s
Isp_rcs = 215;          % s
OF = 0.85;              % -
rho_f = 1010;           % kg/m^3
rho_ox = 1442;          % kg/m^3
n_tank_f = 4;           % -
n_tank_ox = 2;          % -
n_tank_He = 2;          % -
p_tank = 2.1 + 5e-2;    % MPa
gamma_He = 7/5;         % -
R_He = 8314 / 4;        % J/(kg*K)
T_tank = 293;           % K

% Materials for tanks
sigma_Ti = 950;         % MPa
rho_Ti = 4500;          % kg/m^3
sigma_Al = 510;         % MPa
rho_Al = 2810;          % kg/m^3
sigma_Fe = 1400;        % MPa
rho_Fe = 8100;          % kg/m^3

sigma_tank = sigma_Ti;
rho_tank = rho_Ti;
sigma_tank_He = sigma_Ti;
rho_tank_He = rho_Ti;


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

M_f = 1 / (OF + 1) * M_me + M_rcs;
M_ox = OF / (OF + 1) * M_me;


%% Compute volume for tanks of propellant

V_f_tot = 1.1 * M_f / rho_f;
V_ox_tot = 1.1 * M_ox / rho_ox;

V_tank_f = V_f_tot / n_tank_f;
V_tank_ox = V_ox_tot / n_tank_ox;
V_tank = max(V_tank_f, V_tank_ox);


%% Compute geometry and mass for tanks of propellant

r_tank = (3/(4*pi) * V_tank) ^ (1/3);
t_tank = r_tank * p_tank / (2*sigma_tank);
M_tank = 4/3 * pi * rho_tank * ( (r_tank + t_tank)^3 - r_tank^3 );


%% Compute total mass and volume for helium

M_He = 1.2 * p_tank*1e6 * (n_tank_f + n_tank_ox) * V_tank * gamma_He / (0.9 * R_He * T_tank);
V_He_tot = M_He * R_He * T_tank / (10*p_tank*1e6);
V_tank_He = V_He_tot / n_tank_He;


%% Compute geometry and mass for tanks of helium

r_tank_He = (V_tank_He / (2*pi)) ^ (1/3);
h_tank_He = V_tank_He / (r_tank_He^2 * pi);

t_tank_He = r_tank_He * 10*p_tank / (2*sigma_tank_He);
M_tank_He = rho_tank_He * h_tank_He * pi * ( (r_tank_He + t_tank_He)^2 - r_tank_He^2 ) + ...
            2 * rho_tank_He * t_tank_He * r_tank_He^2 * pi;

%%

X = 2.493424054970258E+04; Y = 1.821577789789455E+04; Z =-2.132517705954675E+03;
VX= 4.154120504420773E+00; VY= 6.253043930624797E+00; VZ=-7.603324497195763E-02;

c3_centaur = 2*(0.5*(VX^2 + VY^2 + VZ^2) - astroConstants(13)/sqrt(X^2 + Y^2 + Z^2))

X = 2.495950334920351E+04; Y = 1.825353174540417E+04; Z =-2.188233066043516E+03;
VX= 4.168895408026662E+00; VY= 6.277232780781191E+00; VZ=-1.176509285910803E-01;

c3_juno = 2*(0.5*(VX^2 + VY^2 + VZ^2) - astroConstants(13)/sqrt(X^2 + Y^2 + Z^2))
