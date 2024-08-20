clear, clc, close all;

load("Ephemeris\Radius\Science\Juno_Jupiter_science_extended.mat")
load("Ephemeris\Radius\Science\Juno_Sun_science_extended.mat")

r_mod_JJ = vecnorm(rv_JJ_science_ext, 2, 2);    % km
r_mod_JS = vecnorm(rv_JS_science_ext, 2, 2);    % km

AU = astroConstants(2);                         % km
q0 = astroConstants(31);                        % W/m^2
R_J = astroConstants(25);                       % km
sigma = 5.67e-8;                                % W/(m^2*K^4)

a = 0.45;                                       % -
T_J = 110;                                      % K


%% q calculation

q_sun = q0 * (AU ./ r_mod_JS).^2;
q_alb = a .* q_sun .* (R_J ./ r_mod_JJ).^2;
q_ir  = sigma * T_J^4 * (R_J ./ r_mod_JJ).^2;

q = q_sun + q_alb + q_ir;

[q_min, i_min] = min(q);


%% Plot

time = (0:length(r_mod_JJ)-1) / 8;              % days

linewdth = 1;
fontsz = 12;

figure
hold on
plot(time, q, 'LineWidth', 1)
plot(time(i_min), q_min, 'bx', 'MarkerSize', 8, 'LineWidth', linewdth)
xline(time(2930), 'r--', 'LineWidth', linewdth)
box on
grid minor
xlabel("Time [days]")
ylabel("q_{tot} [W/m^2]")
legend('Total flux', 'Minimum flux for TP-7', 'JOI')
set(gca, 'FontSize', fontsz)