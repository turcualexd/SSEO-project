clear, clc, close all;

% Ephemeris of Juno wrt Sun (whole mission, 1 day resolution)
load("Ephemeris\Juno_Sun.mat");

% Constants
R_earth = astroConstants(23);
q0 = astroConstants(31);
muS = astroConstants(4);
a = 0.35;


%% q_sun for whole mission

r_vec = zeros(3, length(av));
r_mod = zeros(1, length(av));
q_sun  = zeros(1, length(av));

for i = 1:length(av)
    [r_vec(:,i), ~] = kep2car([av(i) ev(i) iv(i) OMv(i) omv(i) thv(i)], muS);
    r_mod(i) = norm(r_vec(:,i)) / astroConstants(2);    % norm distance [AU]
    q_sun(i) = q0 / r_mod(i)^2;                         % flux of power from Sun
end


%% Perihelion

% 0.8 AU case - at perihelion
q_ph_tot = max(q_sun);


%% Plot

days = 0:length(av)-1;
days = days - 661;

linewdth = 1;
fontsz = 12;

figure
hold on
plot(days, q_sun, 'LineWidth', 1)
xline(757-661,'r--', 'LineWidth', linewdth)
plot(757-661, max(q_sun), 'rx', 'MarkerSize', 8, 'LineWidth', linewdth)
%xline(767,'g--')    % EGA
%xline(3,'--')       % Begin IC-1
%xline(63,'--')      % IC-1 -> IC-2
%xline(661,'--')     % IC-2 -> IC-3
%xline(822,'--')     % IC-3 -> OC
%xline(1620,'--')    % End OC
xlim([0 767-661])
box on
grid minor
xlabel("Time [days]")
ylabel("q_{sun} [W/m^2]")
legend('Solar flux', 'Perihelion', 'Maximum flux for TP-3')
set(gca, 'FontSize', fontsz)