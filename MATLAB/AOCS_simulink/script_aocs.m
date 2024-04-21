clear
close all
clc

I = [11245.08 0 0; -61.66 10044.71 0; -61.66 0 17593.63];
I_inv = inv(I);
A0 = [-1 0 0; 0 1 0; 0 0 1];
n = 2* 2*pi/60;    % spin rate 2 rpm
w0 = [0 0 n];

t_s = 1;
toll = deg2rad(15);
t0 = 0;
step_t = 0.5*t_s;
t_f =  20*60; 

braccio = 2;
p_slew = [-.05 + .002i, -.05 - .002i, -.03 + .003i, -.03 - .003i, -.08 + 0.001i, -.08 - 0.001i];

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

%% 
w = results.w';
alphas = results.alphas;
t = results.tout;

figure
hold on
grid minor
plot(t, w)

figure
hold on
grid minor
plot(t, alphas)
legend('err_x', 'err_y', 'err_z')

