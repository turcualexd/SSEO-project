%% inner cruise 1
clear
clc
close all


% ------ calcolo distanze dal sole nella inner cruise 

T = readtable("horizons_results.txt",MissingRule="omitrow");
giorni = [1 61];
clc
ev  = table2array(T(1:60,3));
rpv = table2array(T(1:60,4));
iv  = wrapTo2Pi(deg2rad(table2array(T(1:60,5))));
OMv = wrapTo2Pi(deg2rad(table2array(T(1:60,6))));
omv = wrapTo2Pi(deg2rad(table2array(T(1:60,7))));
av  = table2array(T(1:60,12));
thv = wrapTo2Pi(deg2rad(table2array(T(1:60,11))));
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

r_medio = trapz(rn)/(length(rn));

figure
plot(rn)
hold on
plot(0*rn + r_medio)
title('DIstance for SRP approximation')
grid minor

Sr = 1367/5.2^2; % W/m2, calcolato con la media della distanza

Sr_vect = 1367./(rn).^2;
Sr_medio = trapz(Sr_vect)/length(rn); % calcolato come media delle radiazioni
% ----- Sr_medio è quello da usare per la inner----

figure
plot(Sr_vect)
hold on
plot(0*rn + Sr_medio)
plot(0*rn + Sr)
title('Sr approximation')
grid minor
legend('Sr as function of time', 'Sr medium', 'Sr with medium distance' )

%close all
delta_c_12 = sqrt(5.8^2 - 1.5^2);

A_cs = 70; % area frontale 
q = 0.55; % riflettività da slides
% siamo in sun poiting quasi tutto il tmepo, quindi cos(I) = 1;
% distanze calcolate sul modello solidworks semplificato 
T_srp_inner = Sr/astroConstants(5) * A_cs/3 * ( 1 + q ) * ( 2 * delta_c_12 + 4.5);

% ------------ simulazione -------------
I = [11245.08 0 0; 0 10044.71 0; 0 0 17593.63];             % inertia matrix, solo assi principali considerati, dire che usiamo i pannelli per ottenerli corretti
I_inv = inv(I);
A0 = [0 0 1; 0 1 0; -1 0 0];                                % assetto iniziale 
n = 2 * pi/60;                                              % spin rate 1 rpm
w0 = [0 0 n];    
w0_nominal = [0 0 n];  % velocitè angolare iniziale

braccio_x = 2.7;
braccio_y = 2.7;
braccio_z = 3.2;

giorni_no_controllo = 20;

kx = (I(3,3) - I(2,2))/I(1,1);
ky = (I(3,3) - I(1,1))/I(2,2); 
kz = (I(2,2) - I(1,1))/I(3,3);
% *------------- parametri simulazione senza controllo  ----------------
t_s = 1;
toll = deg2rad(15);
t0 = 0;
step_t = 0.5*t_s;
t_f =  (19 * 24 + 23.5) * 3600;  
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

% ---------- simulaizone per allineare assi

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

t_f =  30*60; 

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

% fuel consumption computation in inner 1

% max_x = max(abs(M_c(:,1)));
% max_y = max(abs(M_c(:,2)));
% max_z = max(abs(M_c(:,3)));

% max_x = trapz(abs(M_c(:,1)))/length(M_c(:,1));
% max_y = trapz(abs(M_c(:,2)))/length(M_c(:,2));
% max_z = trapz(abs(M_c(:,3)))/length(M_c(:,3));

% consumo_x = 0;
% consumo_y = 0;

% for i = 1:length(M_c)
% 
%     if abs(M_c(i, 1)) > max_x %&& max_x > 1e-1
%         consumo_x = consumo_x + 1;
%     end
% 
%      if abs(M_c(i, 2)) > max_y %&& max_y > 1e-1
%         consumo_y = consumo_y + 1;
%      end
% 
%      if abs(M_c(i, 3)) > max_z %&& max_z > 1e-1
%        consumo_z = consumo_z + 1;
%      end
% end

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
massa_totale_inner1 = (giorni(2)-giorni(1))/giorni_no_controllo * ( m_prop_singola_ora_inner + m_prop_z_continuo)

%% inner cruise 2 
% clear
clc
close all

% ------ calcolo distanze dal sole nella inner cruise 

T = readtable("horizons_results.txt",MissingRule="omitrow");
giorni = [61 661];
clc
ev  = table2array(T(61:659,3));
rpv = table2array(T(61:659,4));
iv  = wrapTo2Pi(deg2rad(table2array(T(61:659,5))));
OMv = wrapTo2Pi(deg2rad(table2array(T(61:659,6))));
omv = wrapTo2Pi(deg2rad(table2array(T(61:659,7))));
av  = table2array(T(61:659,12));
thv = wrapTo2Pi(deg2rad(table2array(T(61:659,11))));
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

r_medio = trapz(rn)/(length(rn));

figure
plot(rn)
hold on
plot(0*rn + r_medio)
title('DIstance for SRP approximation')
grid minor

Sr = 1367/r_medio^2; % W/m2, calcolato con la media della distanza

Sr_vect = 1367./(rn).^2;
Sr_medio = trapz(Sr_vect)/length(rn); % calcolato come media delle radiazioni
% ----- Sr_medio è quello da usare per la inner----

figure
plot(Sr_vect)
hold on
plot(0*rn + Sr_medio)
plot(0*rn + Sr)
title('Sr approximation')
grid minor
legend('Sr as function of time', 'Sr medium', 'Sr with medium distance' )

%close all
delta_c_12 = sqrt(5.8^2 - 1.5^2);

A_cs = 70; % area frontale 
q = 0.55; % riflettività da slides
% siamo in sun poiting quasi tutto il tmepo, quindi cos(I) = 1;
% distanze calcolate sul modello solidworks semplificato 
T_srp_inner = Sr_medio/astroConstants(5) * A_cs/3 * ( 1 + q ) * ( 2 * delta_c_12 + 4.5);

% ------------ simulazione -------------
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
giorni_no_controllo = 20;
t_s = 1;
toll = deg2rad(15);
t0 = 0;
step_t = 0.5*t_s;
t_f =  (19*24 + 23.5 )*60*60; 
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

% ---------- simulaizone di un'ora per allineare assi

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

t_f =  30*60; 

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

% fuel consumption computation in inner 1

% max_x = max(abs(M_c(:,1)));
% max_y = max(abs(M_c(:,2)));
% max_z = max(abs(M_c(:,3)));

% max_x = trapz(abs(M_c(:,1)))/length(M_c(:,1));
% max_y = trapz(abs(M_c(:,2)))/length(M_c(:,2));
% max_z = trapz(abs(M_c(:,3)))/length(M_c(:,3));

% consumo_x = 0;
% consumo_y = 0;

% for i = 1:length(M_c)
% 
%     if abs(M_c(i, 1)) > max_x %&& max_x > 1e-1
%         consumo_x = consumo_x + 1;
%     end
% 
%      if abs(M_c(i, 2)) > max_y %&& max_y > 1e-1
%         consumo_y = consumo_y + 1;
%      end
% 
%      if abs(M_c(i, 3)) > max_z %&& max_z > 1e-1
%        consumo_z = consumo_z + 1;
%      end
% end

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
massa_totale_inner2 = (giorni(2)-giorni(1))/giorni_no_controllo * (m_prop_singola_ora_inner + m_prop_z_continuo )

%% inner cruise 3

% clear
clc
close all

% ------ calcolo distanze dal sole nella inner cruise 

T = readtable("horizons_results.txt",MissingRule="omitrow");
giorni = [660 820];
clc
ev  = table2array(T(660:821,3));
rpv = table2array(T(660:821,4));
iv  = wrapTo2Pi(deg2rad(table2array(T(660:821,5))));
OMv = wrapTo2Pi(deg2rad(table2array(T(660:821,6))));
omv = wrapTo2Pi(deg2rad(table2array(T(660:821,7))));
av  = table2array(T(660:821,12));
thv = wrapTo2Pi(deg2rad(table2array(T(660:821,11))));
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

r_medio = trapz(rn)/(length(rn));

figure
plot(rn)
hold on
plot(0*rn + r_medio)
title('DIstance for SRP approximation')
grid minor

Sr = 1367/r_medio^2; % W/m2, calcolato con la media della distanza

Sr_vect = 1367./(rn).^2;
Sr_medio = trapz(Sr_vect)/length(rn); % calcolato come media delle radiazioni
% ----- Sr_medio è quello da usare per la inner----

figure
plot(Sr_vect)
hold on
plot(0*rn + Sr_medio)
plot(0*rn + Sr)
title('Sr approximation')
grid minor
legend('Sr as function of time', 'Sr medium', 'Sr with medium distance' )

%close all
delta_c_12 = sqrt(5.8^2 - 1.5^2);

A_cs = 70; % area frontale 
q = 0.55; % riflettività da slides
% siamo in sun poiting quasi tutto il tmepo, quindi cos(I) = 1;
% distanze calcolate sul modello solidworks semplificato 
T_srp_inner = Sr_medio/astroConstants(5) * A_cs/3 * ( 1 + q ) * ( 2 * delta_c_12 + 4.5);

% ------------ simulazione -------------
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
giorni_no_controllo = 20;
t_s = 1;
toll = deg2rad(15);
t0 = 0;
step_t = 0.5*t_s;
t_f =  (19 * 24 + 23.5) *60*60; 
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

% ---------- simulaizone di un'ora per allineare assi

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

t_f =  30*60; 

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

% fuel consumption computation in inner 1

% max_x = max(abs(M_c(:,1)));
% max_y = max(abs(M_c(:,2)));
% max_z = max(abs(M_c(:,3)));

% max_x = trapz(abs(M_c(:,1)))/length(M_c(:,1));
% max_y = trapz(abs(M_c(:,2)))/length(M_c(:,2));
% max_z = trapz(abs(M_c(:,3)))/length(M_c(:,3));

% consumo_x = 0;
% consumo_y = 0;

% for i = 1:length(M_c)
% 
%     if abs(M_c(i, 1)) > max_x %&& max_x > 1e-1
%         consumo_x = consumo_x + 1;
%     end
% 
%      if abs(M_c(i, 2)) > max_y %&& max_y > 1e-1
%         consumo_y = consumo_y + 1;
%      end
% 
%      if abs(M_c(i, 3)) > max_z %&& max_z > 1e-1
%        consumo_z = consumo_z + 1;
%      end
% end

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
massa_totale_inner3 = (giorni(2)-giorni(1))/giorni_no_controllo * ( m_prop_singola_ora_inner + m_prop_z_continuo)

massa_totale_inner1
massa_totale_inner2
massa_totale_inner3

massa_totale_carbuarnte = massa_totale_inner1 + massa_totale_inner2 + massa_totale_inner3

