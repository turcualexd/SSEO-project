clear
close all
clc

% ----------- geometria del satellite -------------------
I = [11245.08 0 0; 0 10044.71 0; 0 0 17593.63];             % inertia matrix, solo assi principali considerati, dire che usiamo i pannelli per ottenerli corretti
I_inv = inv(I);
A0 = [0 0 1; 0 1 0; -1 0 0];                                % assetto iniziale 
n = 2 * pi/60;                                              % spin rate 1 rpm
w0 = [0 0 n];                                               % velocitè angolare iniziale

braccio_x = 2.7;
braccio_y = 2.7;
braccio_z = 3.2;

kx = (I(3,3) - I(2,2))/I(1,1);
ky = (I(3,3) - I(1,1))/I(2,2); 
kz = (I(2,2) - I(1,1))/I(3,3);
% *------------- parametri simulazione ----------------
t_s = 1;
toll = deg2rad(15);
t0 = 0;
step_t = 0.5*t_s;
t_f =  40*60; 
gain_omega = -9500;
gain_alphas = [-100 -100 -1]';


% --------------------------------- simulink inizialization
% -------------------------

model = "simulink_aocs";
load_system(model);
set_param("simulink_aocs", 'SimulationMode', 'accelerator');
sim_options.SolverType = 'Fixed-step';
sim_options.Solver = 'ode4';
sim_options.FixedStep = 'step_t';
sim_options.StartTime = 't0';
sim_options.StopTime = 't_f';

tic
results = sim("simulink_aocs", sim_options);
toc

w = results.w';
alphas = results.alphas;
t = results.tout/60;
A_b_n = results.A_bn;
M_c = results.M_c;
target = results.target;

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

%% fuel consumption computation

max_x = max(abs(M_c(:,1)));
max_y = max(abs(M_c(:,2)));
max_z = max(abs(M_c(:,3)));

max_x = trapz(abs(M_c(:,1)))/length(M_c(:,1));
max_y = trapz(abs(M_c(:,2)))/length(M_c(:,2));
max_z = trapz(abs(M_c(:,3)))/length(M_c(:,3));

consumo_x = 0;
consumo_y = 0;
consumo_z = 0;

for i = 1:length(M_c)

    if abs(M_c(i, 1)) > max_x
        consumo_x = consumo_x + 1;
    end

     if abs(M_c(i, 2)) > max_y
        consumo_y = consumo_y + 1;
     end

     if abs(M_c(i, 3)) > max_z
       consumo_z = consumo_z + 1;
     end
end


tempo_x = consumo_x * step_t;
tempo_y = consumo_y * step_t;
tempo_z = consumo_z * step_t;


th = 4.5;
I_sp = 250;

m_prop_x = 2*2*th * tempo_x/(I_sp * 9.81)
m_prop_y = 2*2*th * tempo_y/(I_sp * 9.81)
m_prop_z = 2*2*th * tempo_z/(I_sp * 9.81)
m_prop_singola_dsm = 2*(m_prop_x + m_prop_y + m_prop_z)

%% SPIN UP-down mode

clc
close all
clear

spin = [1 5; 5 1; 1 5; 5 1; 1 2; 2 1; 1 5; 5 1; 1 5; 5 2]; % RPM inizio e fine per riga
% spin = [1.1 1]
massa_singolo_spin = zeros(1, size(spin, 1));
massa_totale = 0;
I = [11245.08 0 0; 0 10044.71 0; 0 0 17593.63];            
I_inv = inv(I);
braccio_x = 2.7;
braccio_y = 2.7;
braccio_z = 3.2;
    
kx = (I(3,3) - I(2,2))/I(1,1);
ky = (I(3,3) - I(1,1))/I(2,2); 
kz = (I(2,2) - I(1,1))/I(3,3);

t_s = 1;
toll = deg2rad(15);
t0 = 0;
step_t = 0.5*t_s;
t_f =  20*60; 
gain_omega = -9500;
gain_alphas = -90;  
model = "simulink_aocs_spin";
load_system(model);
set_param("simulink_aocs_spin", 'SimulationMode', 'accelerator');
sim_options.SolverType = 'Fixed-step';
sim_options.Solver = 'ode4';
sim_options.FixedStep = 'step_t';
sim_options.StartTime = 't0';
sim_options.StopTime = 't_f';

for j = 1:size(spin, 1)

    spin_0 = spin(j,1) * 2*pi/60; 
    spin_fin = spin(j,2) * 2*pi/60; 
             
    w0 = [0 0 spin_0]; 
    
    %------------- parametri simulazione ----------------
        
    tic
    results = sim("simulink_aocs_spin", sim_options);
    toc
    
    w = results.w';
    t = results.tout/60;
    % A_b_n = results.A_bn;
    M_c = results.M_c;
    % 
    % h_b = 0 * w;
    % h_n = h_b;
    % 
    % for i = 1 : length(t)
    % 
    %     h_b(:,i) = I * w(:,i); 
    %     h_n(:,i) = A_b_n(:,:,i)' * h_b(:,i);
    % 
    % end
     

    % figure
    % hold on
    % grid minor
    % plot(t, w)
    % legend('w_x', 'w_y', 'w_z')
    % title('Angular velocity')
    % hold off
    % 
    % figure
    % hold on
    % grid minor 
    % plot(t,M_c)
    % legend('M_{xb}', 'M_{yb}', 'M_{zb}')
    % xlabel('time [min]')
    % title('Control moment')
    % hold off

    max_x = max(abs(M_c(:,1)));
    max_y = max(abs(M_c(:,2)));
    max_z = max(abs(M_c(:,3)));
    
    consumo_x = 0;
    consumo_y = 0;
    consumo_z = 0;
    
    for i = 1:length(M_c)
    
        if abs(M_c(i, 1)) == max_x && max_x ~= 0
            consumo_x = consumo_x + 1;
        end
    
         if abs(M_c(i, 2)) == max_y && max_y ~= 0
            consumo_y = consumo_y + 1;
         end
    
         if abs(M_c(i, 3)) == max_z && max_z ~= 0
           consumo_z = consumo_z + 1;
         end
    end
    
    tempo_x = consumo_x * step_t;
    tempo_y = consumo_y * step_t;
    tempo_z = consumo_z * step_t;
    
    th = 4.5;
    I_sp = 250;
    
    m_prop_x = 2 * 2 * th * tempo_x/(I_sp * 9.81);
    m_prop_y = 2 * 2 * th * tempo_y/(I_sp * 9.81);
    m_prop_z = 2 * 2 * th * tempo_z/(I_sp * 9.81);
    m_prop_singola_spin = (m_prop_x + m_prop_y + m_prop_z);

    massa_singolo_spin(j) = m_prop_singola_spin;
    massa_totale = massa_totale + m_prop_singola_spin;

end

massa_singolo_spin
massa_totale

%% disturbance model, vecchia versione

% ------ calcolo distanze dal sole nella inner cruise 

T = readtable("horizons_results.txt",MissingRule="omitrow");
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
rn = vecnorm(R/astroConstants(2), 2); % valore già in AU

r_medio = trapz(rn)/(length(rn));

% figure
% plot(rn)
% hold on
% plot(0*rn + r_medio)
% title('DIstance for SRP approximation')
% grid minor

Sr = 1367/r_medio^2; % W/m2, calcolato con la media della distanza

Sr_vect = 1367./(rn).^2;
Sr_medio = trapz(Sr_vect)/length(rn); % calcolato come media delle radiazioni
% ----- Sr_medio è quello da usare per la inner----

% figure
% plot(Sr_vect)
% hold on
% plot(0*rn + Sr_medio)
% plot(0*rn + Sr)
% title('Sr approximation')
% grid minor
% legend('Sr as function of time', 'Sr medium', 'Sr with medium distance' )

 close all
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
w0 = [0 0 n];                                               % velocitè angolare iniziale

braccio_x = 2.7;
braccio_y = 2.7;
braccio_z = 3.2;

kx = (I(3,3) - I(2,2))/I(1,1);
ky = (I(3,3) - I(1,1))/I(2,2); 
kz = (I(2,2) - I(1,1))/I(3,3);
% *------------- parametri simulazione senza controllo SISTEMALA ORA ----------------
t_s = 1;
toll = deg2rad(15);
t0 = 0;
step_t = 0.5*t_s;
t_f =  24*60*60; 
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

max_x = max(abs(M_c(:,1)));
max_y = max(abs(M_c(:,2)));
max_z = max(abs(M_c(:,3)));

max_x = trapz(abs(M_c(:,1)))/length(M_c(:,1));
max_y = trapz(abs(M_c(:,2)))/length(M_c(:,2));
max_z = trapz(abs(M_c(:,3)))/length(M_c(:,3));

consumo_x = 0;
consumo_y = 0;
consumo_z = 0;

for i = 1:length(M_c)

    if abs(M_c(i, 1)) > max_x %&& max_x > 1e-1
        consumo_x = consumo_x + 1;
    end

     if abs(M_c(i, 2)) > max_y %&& max_y > 1e-1
        consumo_y = consumo_y + 1;
     end

     if abs(M_c(i, 3)) > max_z %&& max_z > 1e-1
       consumo_z = consumo_z + 1;
     end
end

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
massa_totale = (822-3) * m_prop_singola_ora_inner



