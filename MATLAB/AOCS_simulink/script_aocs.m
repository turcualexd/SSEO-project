clear
close all
clc

% I = [11245.08 0 0; -61.66 10044.71 0; -61.66 0 17593.63];
I = [11245.08 0 0; 0 10044.71 0; 0 0 17593.63];
I_inv = inv(I);
A0 = [-1 0 0; 0 -1 0; 0 0 1];
n = 2* pi/60;                                             % spin rate 2 rpm
w0 = [0 0 n];

t_s = 1;
toll = deg2rad(15);
t0 = 0;
step_t = 0.5*t_s;
t_f =  20*60; 

braccio = 2;
p_slew = [-.05 + .002i, -.05 - .002i, -.03 + .003i, -.03 - .003i, -.08 + 0.001i, -.08 - 0.001i];

k_close = 1e3 * [ -1.3597    0.0285    0.0609   -0.0422    0.0015    0.0023;
    0.0002   -0.0698    0.0065    0.0000   -0.0016    0.0003;
    0.1195    0.0177   -0.8508    0.0056   -0.0001   -0.0158];

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
t = results.tout;
A_b_n = results.A_bn;

h_b = 0 * w;
h_n = h_b;

for i = 1 : length(t)
 
    h_b(:,i) = I * w(:,i); 
    h_n(:,i) = A_b_n(:,:,i)' * h_b(:,i);
    % h_t(:,i) = ( AY * A_l_n(time(i)) * AI)' * h_t_b;
    
end

%% 

figure
hold on
grid minor
plot(t, w)
legend('w_x', 'w_y', 'w_z')
title('Angular velocity')

figure
hold on
grid minor
plot(t, alphas)
legend('err_x', 'err_y', 'err_z')

figure
hold on
grid minor
plot(t, h_n)
plot(t, h_b)
legend('h_{xb}', 'h_{yb}', 'h_{zb}', 'h_{xn}', 'h_{yn}', 'h_{zn}')

