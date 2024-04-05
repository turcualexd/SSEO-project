%% Estrai vettore posizione
% clear, clc, close all
% 
% F = fopen("Effemeridi Giove.txt", "r");
% G = fopen("Vettore Giove.txt", "w");
% if F ~= -1 && G ~= -1
%     while ~feof(F)
%         line = fgetl(F);
%         line = fgetl(F);
%         r = [str2double(line(5:26)), str2double(line(31:52)), str2double(line(57:78))];
%         line = fgetl(F);
%         v = [str2double(line(5:26)), str2double(line(31:52)), str2double(line(57:78))];
%         line = fgetl(F);
%         fprintf(G, "%f %f %f %f %f %f\n", r(1), r(2), r(3), v(1), v(2), v(3));
%     end
% end
% fclose(F);
% fclose(G);


%% Earth-Juno distance
clear, clc, close all

r_Juno = load("Vettore Juno.txt");
r_Juno = r_Juno(:, 1:3)';
r_Terra = load("Vettore Terra.txt");
r_Terra = r_Terra(:, 1:3)';
r_Giove = load("Vettore Giove.txt");
r_Giove = r_Giove(:, 1:3)';
dist_TJ = zeros(1, size(r_Juno, 2));
dist_SJ = dist_TJ;
dist_GJ = dist_TJ;
t = (1 : size(r_Juno, 2))/12 + 11/12;
for i = 1 : size(r_Juno, 2)
    dist_TJ(i) = norm(r_Juno(:, i) - r_Terra(:, i))/astroConstants(2);
    dist_SJ(i) = norm(r_Juno(:, i))/astroConstants(2);
    if t(i) < 1796
        dist_GJ(i) = norm(r_Juno(:, i) - r_Giove(:, i))./astroConstants(2);
    else
        dist_GJ(i) = norm(r_Juno(:, i) - r_Giove(:, i))./(100*astroConstants(25));
    end
end
plot(t, dist_TJ, "LineWidth", 1)
grid minor
hold on
plot(t, dist_GJ, "LineWidth", 1)
plot(t, dist_SJ, "LineWidth", 1)
xline(1796, "k-.", "LineWidth", 1)
xline(2264, "k--", "LineWidth", 1)
legend("Terra-Juno", "Giove-Juno", "Sole-Juno", "JOI", "EOM Nominale", "FontSize",3)

%% SPE and SEP angles
clear, clc, close all

r_Juno = load("Vettore Juno.txt");
r_Juno = r_Juno(:, 1:3)';
r_Terra = load("Vettore Terra.txt");
r_Terra = r_Terra(:, 1:3)';
spe = zeros(1, size(r_Juno, 2));
sep = zeros(1, size(r_Juno, 2));
for i = 1 : size(r_Juno, 2)
    r_TJ = r_Juno(:, i) - r_Terra(:, i);
    spe(i) = rad2deg(acos(dot(r_Juno(:,  i), r_TJ)/(norm(r_Juno(:,  i))*norm(r_TJ))));
    sep(i) = rad2deg(acos(dot(-r_Terra(:,  i), r_TJ)/(norm(r_Terra(:,  i))*norm(r_TJ))));
end
t = (1 : size(r_Juno, 2))/12 + 11/12;
plot(t, spe)
hold on
grid minor
plot(t, sep)
xline(796, "k--")
legend("SPE", "SPE", "Fly-by")

%% Sizing
clearvars -except dist_TJ 
clc, close all

% Dati
L_max =9.660609252036554e+08*1e3;   % m
R = 2e5;                    % bps
d_ant = 2.5;                % m
d_ground = 34;              % m
f_x = 8.4*1e9;              % Hz (X band)
f_ka = 34.4*1e9;            % Hz (Ka band)
alpha_enc = 6;              % -
BER_min = 1e-6;             % ??
mu_ant = 0.55;              % -
c = 299792458;              % m/s
eta = 0.1;                  % deg
L_hga_x = 0.25;             % dB
L_hga_ka = 0.5;             % dB
P_hga = 25;                 % W


% Calcoli
R_real = alpha_enc*R;
lambda_x = c/f_x;
lambda_ka = c/f_ka;
G_ant_x = 10*log10(pi^2*d_ant^2*mu_ant/lambda_x^2);
G_ant_ka = 10*log10(pi^2*d_ant^2*mu_ant/lambda_ka^2);
G_ground_x = 10*log10(pi^2*d_ground^2*mu_ant/lambda_x^2);
G_ground_ka = 10*log10(pi^2*d_ground^2*mu_ant/lambda_ka^2);
theta_x = 65.3*lambda_x/d_ground;
theta_ka = 65.3*lambda_ka/d_ground;
L_space_x = 20*log10(lambda_x/(4*pi*L_max));
L_space_ka = 20*log10(lambda_x/(4*pi*L_max));
L_point_x = -12*log10((eta/theta_x)^2);
L_point_ka = -12*log10((eta/theta_ka)^2);

