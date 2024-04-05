clear, clc, close all

% Load DeltaVs
des = load("Manovre design.txt");
dv_des = des(:, 1);
tipo_des = des(:, 2);
vero = load("Manovre vere.txt");
dv_vero = vero(:, 1);
tipo_vero = vero(:, 2);

% Dati
g0 = 9.81;              % m/s^2
m = 1593;               % kg
Is_me = 318.6;          % s
Is_rcs = 220;           % s
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

sigma_tank = sigma_Fe;
rho_tank = rho_Fe;
sigma_tank_He = sigma_Fe;
rho_tank_He = rho_Fe;
%%
% Backwards design masses with margins (from dry mass)
m_p_rcs = 0;
m_p_me = 0;
for i = length(dv_des) : -1 : 1
    if tipo_des(i)
        m_p = m*(exp(1.05*dv_des(i)/(Is_me*g0)) - 1);
        m = m + m_p;
        m_p_me = m_p_me + m_p;
    else
        m_p = m*(exp(2*dv_des(i)/(Is_rcs*g0)) - 1);
        m = m + m_p;
        m_p_rcs = m_p_rcs + m_p;
    end
end

% Actual used masses up to now
m = 3625;
m_f_vero = 0;
m_ox_vero = 0;
for i = 1 : length(dv_vero)
    if tipo_vero(i) % main engine
        m_p = m*(1 - exp(-dv_vero(i)/(Is_me*g0)));
        m = m - m_p;
        m_f_vero = m_f_vero + m_p/(1 + OF);
        m_ox_vero = m_ox_vero + m_p*OF/(1 + OF);
    else % rcs
        m_p = m*(1 - exp(-dv_vero(i)/(Is_rcs*g0)));
        m = m - m_p;
        m_f_vero = m_f_vero + m_p;
    end
end

M_f = m_p_me/(1 + OF) + m_p_rcs
M_ox = m_p_me*OF/(1 + OF)
% m_tot_back = m_f_des + m_ox_des;
% m_tot_vero = m_f_vero + m_ox_vero;


V_f_tot = 1.1 * M_f / rho_f;
V_ox_tot = 1.1 * M_ox / rho_ox;

V_tank_f = V_f_tot / n_tank_f;
V_tank_ox = V_ox_tot / n_tank_ox;
V_tank = max(V_tank_f, V_tank_ox);


%% Compute geometry and mass for tanks of propellant

r_tank = (3/(4*pi) * V_tank) ^ (1/3);
t_tank = r_tank * p_tank / (2*sigma_tank) * 1000
M_tank = 4/3 * pi * rho_tank * ( (r_tank + t_tank/1000)^3 - r_tank^3 )


%% Compute total mass and volume for helium

M_He = 1.2 * p_tank*1e6 * (n_tank_f + n_tank_ox) * V_tank * gamma_He / (0.9 * R_He * T_tank);
V_He_tot = M_He * R_He * T_tank / (10*p_tank*1e6);
V_tank_He = V_He_tot / n_tank_He;


%% Compute geometry and mass for tanks of helium

r_tank_He = (V_tank_He / (2*pi)) ^ (1/3);
h_tank_He = V_tank_He / (r_tank_He^2 * pi);

t_tank_He = r_tank_He * 10*p_tank / (2*sigma_tank_He) * 1000
t_tank_He = t_tank_He/1000;
M_tank_He = rho_tank_He * h_tank_He * pi * ( (r_tank_He + t_tank_He)^2 - r_tank_He^2 ) + ...
            2 * rho_tank_He * t_tank_He * r_tank_He^2 * pi
