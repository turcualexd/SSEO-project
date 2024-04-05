clear, clc, close all

%% Initial data

manovre = load("Manovre.txt");
dv = manovre(:, 1);
type = manovre(:, 2);

OF = 0.85;
m_sc_in = 3625;
g0 = 9.81;
p_prop = 2.1 + 5e-2;    % MPa

Is_me = 318.6;
Is_rcs = 220;

rho_f = 1010;
rho_ox = 1442;

S_Ti = 950;     % MPa
rho_Ti = 4500;  % kg/m^3

S_Al = 503;     % MPa
rho_Al = 2810;  % kg/m^3

S_Fe = 1100;   % MPa
rho_Fe = 7900; % kg/m^3

%% Fuel and oxidizer tanks sizing

dv_me = 0;
dv_rcs = 0;

for i = 1 : length(dv)
    if type(i)
        dv_me = dv_me + dv(i);
    else
        dv_rcs = dv_rcs + dv(i);
    end
end

dv_me = 1.1 * dv_me;
dv_rcs = 1.1 * dv_rcs;

m_me = 1.055 * m_sc_in * (1 - exp( -dv_me / (Is_me*g0) ));
m_rcs = 1.055 * m_sc_in * (1 - exp( -dv_rcs / (Is_rcs*g0) ));

m_f = m_me / (1 + OF) + m_rcs;
m_ox = m_me * OF / (1 + OF);

V_tot_f = 1.1 * m_f / rho_f;
V_tot_ox = 1.1 * m_ox / rho_ox;

r_f = (3 * V_tot_f / (16*pi)) ^ (1/3);
r_ox = (3 * V_tot_ox / (8*pi)) ^ (1/3);
r_1tank = max([r_f r_ox]);

V_tot_tank = 6 * 4/3 * pi * r_1tank^3;

t_Ti = 1000 * p_prop * r_1tank / (2*S_Ti);  % mm
t_Al = 1000 * p_prop * r_1tank / (2*S_Al);  % mm
t_Fe = 1000 * p_prop * r_1tank / (2*S_Fe);  % mm

m_tank_Ti = 4/3 * pi * rho_Ti * ((r_1tank + t_Ti*1e-3)^3 - r_1tank^3); % kg
m_tank_Al = 4/3 * pi * rho_Al * ((r_1tank + t_Al*1e-3)^3 - r_1tank^3); % kg
m_tank_Fe = 4/3 * pi * rho_Fe * ((r_1tank + t_Fe*1e-3)^3 - r_1tank^3); % kg


%% Helium tank sizing

R_He = 8314 / 4;    % J/(kg*K)
T_He = 100;         % K
n = 3;              % number of tanks

V_tot_tank_He_vec = linspace(0.001, 0.5, 10000).';                  % m^3
p_tank_He_vec = p_prop .* (1 + V_tot_tank ./ V_tot_tank_He_vec);    % MPa

r_He = (3 ./ (4*pi*n) .* V_tot_tank_He_vec) .^ (1/3);               % m
t_He_Ti = 1000 .* p_tank_He_vec .* r_He ./ (2*S_Ti);                % mm

m_He_Ti = 4/3 * pi * rho_Ti * ((r_He + t_He_Ti*1e-3).^3 - r_He.^3); % kg
m_He = p_tank_He_vec * 1e6 .* V_tot_tank_He_vec/n / (R_He * T_He);  % kg
m_He_tot = (m_He_Ti + m_He);  % mass of a tank full of He at T_He

[min_mass, index] = min(m_He_tot);
index = index + 2000;
min_radius = r_He(index);
min_volume = V_tot_tank_He_vec(index);
min_t = t_He_Ti(index)
min_pressure = p_tank_He_vec(index)
min_He_mass = m_He(index);

tot_height = (min_radius + r_1tank) * 2