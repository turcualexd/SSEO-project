clear; close all; clc;

T = readtable("horizons_results.txt",MissingRule="omitrow");
% ------ calcolo distanze dal sole nella inner cruise

giorni = [1 60];
clc
ev  = table2array(T(giorni(1):giorni(end),3));
rpv = table2array(T(giorni(1):giorni(end),4));
iv  = wrapTo2Pi(deg2rad(table2array(T(giorni(1):giorni(end),5))));
OMv = wrapTo2Pi(deg2rad(table2array(T(giorni(1):giorni(end),6))));
omv = wrapTo2Pi(deg2rad(table2array(T(giorni(1):giorni(end),7))));
av  = table2array(T(giorni(1):giorni(end),12));
thv = wrapTo2Pi(deg2rad(table2array(T(giorni(1):giorni(end),11))));
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

rn = vecnorm(R/astroConstants(2), 2); % valore già in AU
Sr_vect = 1367./(rn).^2;
Sr_medio = trapz(Sr_vect)/length(rn);

% figure
% plot(Sr_vect)
% hold on
% plot(0*rn + Sr_medio)
% 
% title('Sr approximation')
% grid minor
% legend('Sr as function of time', 'Sr medium' )

delta_c_12 = sqrt(5.8^2 - 1.5^2);
A_cs = 70; % area frontale 
q = 0.55; % riflettività da slides
% siamo in sun poiting quasi tutto il tmepo, quindi cos(I) = 1;
% distanze calcolate sul modello solidworks semplificato 
T_srp_inner = Sr_medio/astroConstants(5) * A_cs/3 * ( 1 + q ) * ( 2 * delta_c_12 + 4.5);
I = [11245.08 0 0; 0 10044.71 0; 0 0 17593.63];             % inertia matrix, solo assi principali considerati, dire che usiamo i pannelli per ottenerli corretti
I_inv = inv(I);
A0 = [0 0 1; 0 1 0; -1 0 0];                                % assetto iniziale 
n = 2 * pi/60;                                              % spin rate 1 rpm
w0 = [0 0 n];    
w0_nominal = [0 0 n];  % velocitè angolare iniziale
braccio_x = 2.7;
braccio_y = 2.7;
braccio_z = 3.2;

kx = (I(3,3) - I(2,2))/I(1,1);
ky = (I(3,3) - I(1,1))/I(2,2); 
kz = (I(2,2) - I(1,1))/I(3,3);

% *------------- parametri simulazione senza controllo SISTEMALA ORA ----------------
giorni_no_controllo = 14;
t_s = 1;
toll = deg2rad(15);
t0 = 0;
step_t = 0.5*t_s;
t_f =  giorni_no_controllo*24*60*60; 
gain_omega = [0 0 -95];
gain_alphas = 0;
gain_integ = 0;

model = "disturbi_inner_sim";
load_system(model);
set_param("disturbi_inner_sim", 'SimulationMode', 'accelerator');
sim_options.SolverType = 'Fixed-step';
sim_options.Solver = 'ode4';
sim_options.FixedStep = 'step_t';
sim_options.StartTime = 't0';
sim_options.StopTime = 't_f';

tic
results = sim("disturbi_inner_sim", sim_options);
toc

w = results.w';
alphas = results.alphas;
t = results.tout/60;
A_b_n = results.A_bn;
M_c = results.M_c;
target = results.target;
rot_mat = results.rot_mat;
w_err = results.w_err;

h_b = 0 * w;
h_n = h_b;

for i = 1 : length(t)
 
    h_b(:,i) = I * w(:,i); 
    h_n(:,i) = A_b_n(:,:,i)' * h_b(:,i);
      
end

figure
hold on
grid minor
plot(t, w)
legend('w_x', 'w_y', 'w_z')
title('Angular velocity')

figure
hold on
grid minor
plot(t, alphas*180/pi)
legend('err_x', 'err_y', 'err_z')
title('Errors in attitude')
xlabel('time [min]')

figure
hold on
grid minor
plot(t, w_err*180/pi)
legend('w_{err_x}', 'w_{err_y}', 'w_{err_z}')
title('Errors angular speed')
xlabel('time [min]')

figure
hold on
grid minor
plot(t, h_n)
legend('h_{xn}', 'h_{yn}', 'h_{zn}')
xlabel('time [min]')
title('Angular momentum (inertial)')

figure
hold on
grid minor 
plot(t, h_b)
legend('h_{xb}', 'h_{yb}', 'h_{zb}')
xlabel('time [min]')
title('Angular momentum (body)')

figure
hold on
grid minor 
plot(t,M_c)
legend('M_{xb}', 'M_{yb}', 'M_{zb}')
xlabel('time [min]')
title('Control moment')

th = 4.5;
I_sp = 250;

consumo_z = sum(abs(M_c(:,3)))/(th*braccio_z);
tempo_z = consumo_z * step_t;
m_prop_z_continuo = 2*2*th * tempo_z/(I_sp * 9.81)

A0 = A_b_n(:,:,end);
w0 = w(:,end);

gain_omega = -950;
gain_alphas = [-100 -100 -1]';
gain_integ = -5;

model = "disturbi_inner_sim";
load_system(model);
set_param("disturbi_inner_sim", 'SimulationMode', 'accelerator');
sim_options.SolverType = 'Fixed-step';
sim_options.Solver = 'ode4';
sim_options.FixedStep = 'step_t';
sim_options.StartTime = 't0';
sim_options.StopTime = 't_f';

t_f =  20*60; 

tic
results = sim("disturbi_inner_sim", sim_options);
toc

w = results.w';
alphas = results.alphas;
t = results.tout/60;
A_b_n = results.A_bn;
M_c = results.M_c;
target = results.target;
rot_mat = results.rot_mat;
w_err = results.w_err;

h_b = 0 * w;
h_n = h_b;

for i = 1 : length(t)
 
    h_b(:,i) = I * w(:,i); 
    h_n(:,i) = A_b_n(:,:,i)' * h_b(:,i);
      
end

figure
hold on
grid minor
plot(t, w)
legend('w_x', 'w_y', 'w_z')
title('Angular velocity')

figure
hold on
grid minor
plot(t, alphas*180/pi)
legend('err_x', 'err_y', 'err_z')
title('Errors in attitude')
xlabel('time [min]')

figure
hold on
grid minor
plot(t, w_err*180/pi)
legend('w_{err_x}', 'w_{err_y}', 'w_{err_z}')
title('Errors angular speed')
xlabel('time [min]')

figure
hold on
grid minor
plot(t, h_n)
legend('h_{xn}', 'h_{yn}', 'h_{zn}')
xlabel('time [min]')
title('Angular momentum (inertial)')

figure
hold on
grid minor 
plot(t, h_b)
legend('h_{xb}', 'h_{yb}', 'h_{zb}')
xlabel('time [min]')
title('Angular momentum (body)')

figure
hold on
grid minor 
plot(t,M_c)
legend('M_{xb}', 'M_{yb}', 'M_{zb}')
xlabel('time [min]')
title('Control moment')

th = 4.5;
I_sp = 250;

consumo_x = sum(abs(M_c(:,1)))/(th*braccio_x)
consumo_y = sum(abs(M_c(:,2)))/(th*braccio_y)
consumo_z = sum(abs(M_c(:,3)))/(th*braccio_z)

tempo_x = consumo_x * step_t;
tempo_y = consumo_y * step_t;
tempo_z = consumo_z * step_t;

m_prop_x = 2*2*th * tempo_x/(I_sp * 9.81)
m_prop_y = 2*2*th * tempo_y/(I_sp * 9.81)
m_prop_z = 2*2*th * tempo_z/(I_sp * 9.81)
m_prop_singola_ora_inner = (m_prop_x + m_prop_y + m_prop_z)
%%
massa_totale_inner1 = (giorni(2)-giorni(1))/giorni_no_controllo * ( m_prop_singola_ora_inner + m_prop_z_continuo) 

