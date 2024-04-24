clear
clc
close all

% juno = load("juno_table_jupiter.txt");
% 
% 
% for i = 1 : size(R_juno_jupiter, 2)
%     B(i) = 2*4.3e-4 * 71372^3 / norm(R_juno_jupiter(:,i))^3;
% end
% 
% Mm = B * 20;
% 
% R_juno_sun = loadData("horizons_results_sun.txt");
% 
% 
% for i = 1 : size(R_juno_sun, 2)
%     R_js =  norm(R_juno_sun(:,i));
% end

I = diag([11245.08 10044.71 17593.63]);             % inertia matrix, solo assi principali considerati, dire che usiamo i pannelli per ottenerli corretti
I_inv = inv(I);
A0 = [1 0 0; 0 0 1; 0 -1 0]*[0 0 1; 0 1 0; -1 0 0];                                % assetto iniziale 
n = 2 * 2 * pi/60;                                              % spin rate 1 rpm
w0 = [0 0 n];    
w0_nominal = [0 0 n];
mu = astroConstants(15);
vet = load("Vet.txt");
rvet = vet(1, 1:3)';
vvet = vet(1, 4:6)';
el = car2kep(rvet, vvet, mu);
T = 11*24*3600;
a = el(1);
e = el(2);
rp = a*(1 - e);
a0 = (mu*T^2/(4*pi^2))^(1/3);
e0 = 1 - rp/a0;
i0 = el(3);
p = a0*(1 - e0^2);
vp = (1 + e)*sqrt(mu/p);

th0 = deg2rad(180);
braccio_x = 2.7;
braccio_y = 2.7;
braccio_z = 3.2;
j_b = [0.05; 0.05; 0.05];
toll = deg2rad(15);
t0 = 0;
step_t = 0.5;
gain_omega = -1;
gain_alphas = [-1000; -1000; -10];
gain_integ = 0;
&
dist = 1;

model = "planetary_phase";
%setsim
load_system(model);
set_param(model, "StopTime", "T", "SolverName", "ode4", "FixedStep", "step_t", "SimulationMode", "accelerator");
s = sim(model);
w = s.w';
M_c = s.M_c;
alphas = s.alphas;
t = s.tout;
A_b_n = s.A_bn;

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
plot(t, rad2deg(alphas))
legend('err_x', 'err_y', 'err_z')
title('Errors in attitude')
xlabel('time [min]')

figure
hold on
grid minor 
plot(t,M_c)
legend('M_{xb}', 'M_{yb}', 'M_{zb}')
xlabel('time [min]')
title('Control moment')

figure
hold on
grid minor
plot(t, h_n)
legend('h_x', 'h_y', 'h_z')
title('h inerziale')
xlabel('time [min]')

th = 4.5;
I_sp = 250;

consumo_x = sum(abs(M_c(:,1)))/(th*braccio_x);
consumo_y = sum(abs(M_c(:,2)))/(th*braccio_y);
consumo_z = sum(abs(M_c(:,3)))/(th*braccio_z);

tempo_x = consumo_x * step_t;
tempo_y = consumo_y * step_t;
tempo_z = consumo_z * step_t;

n_orb = 33;

m_prop_x = 2*2*th * tempo_x/(I_sp * 9.81);
m_prop_y = 2*2*th * tempo_y/(I_sp * 9.81);
m_prop_z = 2*2*th * tempo_z/(I_sp * 9.81);
m_prop_orb = (m_prop_x + m_prop_y + m_prop_z);
m_p_tot = n_orb*m_prop_orb