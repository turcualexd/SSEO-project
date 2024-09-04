clear, clc, close all;

load("Ephemeris\Radius\EGA\Juno_Earth_EGA.mat")
load("Ephemeris\Radius\EGA\Juno_Sun_EGA.mat")

r_mod_JE = vecnorm(rv_JE_EGA, 2, 2);        % km
r_mod_JS = vecnorm(rv_JS_EGA, 2, 2);        % km

time = 0:length(rv_JS_EGA)-1;               % min

AU = astroConstants(2);                     % km
q0 = astroConstants(31);                    % W/m^2
R_E = astroConstants(23);                   % km
sigma = 5.67e-8;                            % W/(m^2*K^4)

a = 0.35;                                   % -
T_E = 218;                                  % K


%% q calculation

q_sun = q0 * (AU ./ r_mod_JS).^2;
q_alb = a * q0 * (R_E ./ r_mod_JE).^2;
q_ir  = sigma * T_E^4 * (R_E ./ r_mod_JE).^2;

shadow_entry_index = 1*24*60 + 19*60 + 19 + 1;
shadow_exit_index = shadow_entry_index + 19;

q_sun = [q_sun(1:shadow_entry_index-1); zeros(19,1);
    q_sun(shadow_exit_index:end)];
q_alb = [q_alb(1:shadow_entry_index-1); zeros(19,1);
    q_alb(shadow_exit_index:end)];

q = q_sun + q_alb + q_ir;
min(q)


%% Plot

SOI_entry = 1144 - 1;
SOI_exit = 4061 - 1;

time = time - SOI_entry;

linewdth = 1;
fontsz = 10;

figure
hold on
plot(time, q, 'LineWidth', linewdth)
plot(1455, max(q), 'rx', 'MarkerSize', 8, 'LineWidth', linewdth)
plot(1474, min(q), 'bx', 'MarkerSize', 8, 'LineWidth', linewdth)
box on
grid minor
xlabel('Time [min]')
ylabel('q_{tot} [W/m^2]')
%legend('', 'Maximum flux for EGA', 'Minimum flux for EGA')
xlim([0 SOI_exit - SOI_entry])
set(gca, 'FontSize', fontsz)

axes('position',[.25 .18 .2 .5])
box on
hold on
grid minor
indexofinterest= (time>1430) & (time<1500);
plot(time(indexofinterest), q(indexofinterest), 'LineWidth', linewdth)
plot(1455, max(q), 'rx', 'MarkerSize', 8, 'LineWidth', linewdth)
plot(1474, min(q), 'bx', 'MarkerSize', 8, 'LineWidth', linewdth)
yline(1759.2288, 'r--', 'LineWidth', linewdth)
yline(45.6173, 'b--', 'LineWidth', linewdth)
lgnd = legend('', 'Maximum flux for EGA', 'Minimum flux for EGA', ...
    'Selected hot case (TP-3)', 'Selected cold case (TP-7)');
lgnd.FontSize = fontsz;