%clear, clc, close all;

AU = astroConstants(2);                     % km
q0 = astroConstants(31);                    % W/m^2

load("Ephemeris\Radius\Juno_Sun.mat")

r_mod_JS = vecnorm(rv_JS, 2, 2);            % km

q_sun = q0 * (AU ./ r_mod_JS).^2;

% 0.88 AU case - at perihelion
q_ph_tot = max(q_sun);


%% Plot

time = 0:length(q_sun)-1;                   % days
time = time - 661;

linewdth = 1;
fontsz = 12;

figure
hold on
plot(time, q_sun, 'LineWidth', 1)
xline(757-661,'r--', 'LineWidth', linewdth)
plot(757-661, max(q_sun), 'rx', 'MarkerSize', 8, 'LineWidth', linewdth)
%xline(767,'g--')    % EGA
%xline(3,'--')       % Begin IC-1
%xline(63,'--')      % IC-1 -> IC-2
%xline(661,'--')     % IC-2 -> IC-3
%xline(822,'--')     % IC-3 -> OC
%xline(1620,'--')    % End OC
xlim([661 767]-661)
box on
grid minor
xlabel("Time [days]")
ylabel("q_{sun} [W/m^2]")
legend('Solar flux', 'Perihelion', 'Maximum flux for TP-3')
set(gca, 'FontSize', fontsz)