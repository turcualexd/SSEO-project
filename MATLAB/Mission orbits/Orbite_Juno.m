%% Creo file Elementi Juno
%{
Creo il file con solo gli elementi e mjd2000 dalle effemeridi tirate giù da
nasa horizons, in teoria non serve runnarlo
%}
clear, clc, close all

F = fopen("Effemeridi Juno.txt", "r");
G = fopen("Elementi Juno.txt", "w");
if F ~= -1 && G ~= -1
    while ~feof(F)
        line = fgetl(F);
        year = str2double(line(26:29));
        switch line(31:33)
            case "Jan"
                month = 1;
            case "Feb"
                month = 2;
            case "Mar"
                month = 3;
            case "Apr"
                month = 4;
            case "May"
                month = 5;
            case "Jun"
                month = 6;
            case "Jul"
                month = 7;
            case "Aug"
                month = 8;
            case "Sep"
                month = 9;
            case "Oct"
                month = 10;
            case "Nov"
                month = 11;
            case "Dec"
                month = 12;
        end
        day = str2double(line(35:36));
        hrs = str2double(line(38:39));
        min = str2double(line(41:42));
        sec = str2double(line(44:50));
        date = date2mjd2000([year month day hrs min sec]);
        line = fgetl(F);
        e = str2double(line(6:26));
        i = str2double(line(58:78));
        line = fgetl(F);
        OM = str2double(line(6:26));
        om = str2double(line(32:52));
        line = fgetl(F);
        theta = str2double(line(58:78));
        line = fgetl(F);
        a = str2double(line(6:26));
        fprintf(G, "%f %f %f %f %f %f %f\n", date, a, e, i, OM, om, theta);
    end
end

fclose(G);
fclose(F);

%% Eventi principali
%{

Launch        : Aug 5,  2011 16:25 UTC 
  Separation    : Aug 5,  2011 17:24:56.12  UTC
  DSM-1         : Aug 30, 2012 22:29        UTC (344.95 m/s)
  DSM-2         : Sep 3,  2012 22:29        UTC (387.77 m/s)
  Shadow entry  : Oct 9,  2013 19:19:36     UTC (alt.= 7069 km, umbral)
  Earth Flyby   : Oct 9,  2013 19:21:25     UTC (alt.= 558.848) 
                   (lat=-34.170 deg,E.long= 34.008 deg
  Shadow exit   : Oct 9,  2013 19:39:01     UTC (alt.= 14896 km, umbral)
  Arrive Jupiter: Jul 5,  2016 02:30 UTC    (JOI to 53.5-day capture orbit, 
                                             delta-v= 541.7 m/s)
  Period Reduction Maneuver:
                  Oct 19, 2016 (CANCELLED)  (Target 14-day orbit period,
                                             delta-v= 395.2 m/s)
  Ganymede flyby: Jun 07, 2021 (PJ34)       (1038 km from surface, per. 53->43 d) 
  Europa flyby  : Sep 29, 2022 (PJ45)       Period from 43->38 days
  Io flybys     : Dec 30, 2023 (PJ57)
                  Feb  3, 2024 (PJ58)       Period reduced to 33 days
  End of mission: Sep   , 2025 (PJ76)       Jupiter impact, 700 km below
                                             1-bar atmospheric pressure level)

%}
%% Plot singoli elementi
%{
Plot di tutti gli elementi dal giorno dopo il lancio al giorno di arrivo a
giove con passo di un giorno. Le linee tratteggiate indicano i momenti in
cui sono avvenute le due deep space maneuvres a il flyby secondo le date
riportate sopra. Sembrano tornare a parte DSM2
%}

clear, clc, close all

el = load("Elementi Juno.txt");
t = 1 : size(el, 1);
label = {"a [km]", "e [-]", "i [deg]", "OM [deg]", "om [deg]", "theta [deg]"};
tit = {"Semiasse maggiore", "Eccentricità", "Inclinazione", "RAAN", "Anomalia del preicentro", "Anomalia vera"};
dsm1 = 390;
dsm2 = 394; 
fb = 795; 
for i = 2 : 7
    figure
    plot(t, el(:, i), "LineWidth", 1);
    grid minor
    hold on
    xline(dsm1, 'k--', 'LineWidth', 1);
    xline(dsm2, 'k--', 'LineWidth', 1);
    xline(fb, 'k--', 'LineWidth', 1);
    xlabel("time [days]")
    ylabel(label{i-1})
    title(tit{i-1})
end

%% Plot traiettoria
clear, clc, close all

ch = 16;
lw = 2;
sun = 30;
jup = 200;
earth = 1500;
assi = 16;

el = load("Elementi Juno.txt");
el = [el(:, 2:3) deg2rad(el(:, 4:7))];
mu = astroConstants(4);
rvet = zeros(3, size(el, 1));
for i = 1 : length(rvet)
    rvet(:, i) = kep2car(el(i, :), mu);
end
plotOrbitVet(rvet, lw, "#41bf13")
hold on
grid minor
axis equal
drawPlanet(10, [0 0 0], sun, 0)
dep = date2mjd2000([2011 7 5 16 25 0]);
fb = date2mjd2000([2013 10 9 19 19 36]);
arr = date2mjd([2016 10 5 2 30 00]);
[dep_pos, v_e] = kep2car(uplanet(dep, 3), mu);
fb_pos = kep2car(uplanet(fb, 3), mu);
[arr_pos, v_j] = kep2car(uplanet(arr, 5), mu);
drawPlanet(3, rvet(:, 1), earth, 0);
drawPlanet(3, fb_pos, earth, 0);
drawPlanet(5, rvet(:, end), jup, 0);
t_e = linspace(0, 3600*24*365, 1e5);
r_e = propOrbit(dep_pos, v_e, mu, t_e);
t_j = linspace(0, 1.6*3600*24*365, 1e5);
r_j = propOrbit(arr_pos, v_j, mu, t_j);
plotOrbitVet(r_e, lw, 'b')
t_m = linspace(0, 2*3600*24*365);
[r_m, v_m] = kep2car(uplanet(dep, 4), mu);
r_m = propOrbit(r_m, v_m, mu, t_m);
plotOrbitVet(r_m, lw, 'r')
plotOrbitVet(r_j, lw, [0.9290 0.6940 0.1250])
legend("Juno trajectory", "", "", "", "", "Earth orbit", "Mars orbit", "Jupiter orbit", 'FontSize', ch)
xlabel("[km]", 'Interpreter', 'latex', 'FontSize', assi)
ylabel("[km]", 'Interpreter', 'latex', 'FontSize', assi)