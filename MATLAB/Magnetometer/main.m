%% ORBIT PARAMETERS
RAAN    =  0;                % Right Ascension of Ascendent Node [deg]
w       =  350;                % Argument of perigee               [deg]
v0      =  0;                % True anomaly at the departure     [deg]
i       =  50;             % inclination                       [deg]
a       =  6828;              % Major semi-axis           (>6378) [km]
e       =  0.001;             % Eccentricity
start_time = datetime(2019,05,15,00,00,00); % UTC time of sattelite starting point
norb = 100;                     % number of orbits
time_step = 60;               % Calculate point every time_step   [s],
                              %   decrease for faster calculation

%% CALCULATION
[lla, time] = orbit_calc(RAAN, w, v0, i, a, e, start_time, norb, ...
    time_step);

%  B   Column vectors of Magnetic field vector     [nT]
%  H   Horizontal component of Magnetic field  B   [nT]
%  D   Magnetic field Declination                  [deg]
%  I   Magnetic field Inclanation                  [deg]
%  F   Magnetic field Intensity                    [nT] 
[B, H, D, I, F] = b_calc(lla, time);

%% PLOTTING
subplot(3,1,1)
plot(F);
title('F - Magnetic field Intensity ')
xlabel('point num') 
ylabel('[nT]') 

subplot(3,1,2)
plot(D);
title('D - Magnetic field Declination')
xlabel('point num') 
ylabel('[deg]')  

subplot(3,1,3)
plot(I);
title('I - Magnetic field Inclanation')
xlabel('point num') 
ylabel('[deg]') 