close all
clear all

dtusat = readtable('Sun.csv');
data = dtusat;
% Number of time steps to display
steps = 800;
% Get data from file
theta = data(:,1);
theta = table2array(theta);
phi = data(:,2);
phi = table2array(phi);
% Conversion of angles from degrees to radians
theta_rad = theta*pi/180;
phi_rad = phi*pi/180;
for n = 1:length(theta)
% Determination of which angles to calculate
    if theta(n) <= 90 && theta(n) >= -90
     x = 1;
    else
     x = -1;
    end
    if theta(n) >= 0 && theta(n) < 180
     y = 1;
    else
     y = -1;
    end
    if phi(n) >= 0
     z = 1;
    else
     z = -1;
    end
    % Calculation of sun sensor angles
    if x == 1
     xy(n) = theta_rad(n);
     xz(n) = atan(tan(phi_rad(n))/cos(theta_rad(n)));
     mxy(n) = -pi/2;
     mxz(n) = -pi/2;
    else
     xy(n) = -pi/2;
     xz(n) = -pi/2;
     mxy(n) = pi-theta_rad(n);
     if mxy(n) > pi
     mxy(n) = -(2*pi-mxy(n));
     end
     mxz(n) = -atan(tan(phi_rad(n))/cos(theta_rad(n)));
    end
    if y == 1
     yx(n) = pi/2-theta_rad(n);
     yz(n) = atan(tan(phi_rad(n))/sin(theta_rad(n)));
     myx(n) = -pi/2;
     myz(n) = -pi/2;
    else
     yx(n) = -pi/2;
     yz(n) = -pi/2;
     myx(n) = pi/2-theta_rad(n);
     if myx(n) > pi/2
     myx(n) = pi-myx(n);
     end
     myz(n) = -atan(tan(phi_rad(n))/sin(theta_rad(n)));
    end
    if z == 1
     zx(n) = pi/2-atan(tan(phi_rad(n))/cos(theta_rad(n)));
     if zx(n) > pi/2
     zx(n) = zx(n)-pi;
     end
     zy(n) = pi/2-atan(tan(phi_rad(n))/sin(theta_rad(n)));
     if zy(n) > pi/2
     zy(n) = zy(n)-pi;
     end
     mzx(n) = -pi/2;
     mzy(n) = -pi/2;
     else
     zx(n) = -pi/2;
     zy(n) = -pi/2;
     mzx(n) = pi/2-atan(tan(-phi_rad(n))/cos(theta_rad(n)));
     if mzx(n) > pi/2
     mzx(n) = mzx(n)-pi;
     end
     mzy(n) = pi/2-atan(tan(-phi_rad(n))/sin(theta_rad(n)));
     if mzy(n) > pi/2
     mzy(n) = mzy(n)-pi;
     end
    end

end


% Conversion of angles from radians to degrees
xy = xy*180/pi;
xz = xz*180/pi;
yx = yx*180/pi;
yz = yz*180/pi;
zx = zx*180/pi;
zy = zy*180/pi;
mxy = mxy*180/pi;
mxz = mxz*180/pi;
myx = myx*180/pi;
myz = myz*180/pi;
mzx = mzx*180/pi;
mzy = mzy*180/pi;


% Constraints on sensors due to exceeded field of view in sensor direction
for n = 1:length(xy)
 % Plane x

 if xy(n) < -45 || xy(n) > 45
 xy(n) = -90;
 end
 if xz(n) < -45 || xz(n) > 45
 xz(n) = -90;
 end
 % Plane y

 if yx(n) < -45 || yx(n) > 45
 yx(n) = -90;
 end
 if yz(n) < -45 || yz(n) > 45
 yz(n) = -90;
 end
 % Plane z

 if zx(n) < -45 || zx(n) > 45
 zx(n) = -90;
 end
 if zy(n) < -45 || zy(n) > 45
 zy(n) = -90;
 end
 % Plane -x
 if mxy(n) < -45 || mxy(n) > 45
 mxy(n) = -90;
 end
 if mxz(n) < -45 || mxz(n) > 45
 mxz(n) = -90;
 end
 % Plane -y
 if myx(n) < -45 || myx(n) > 45
 myx(n) = -90;
 end
 if myz(n) < -45 || myz(n) > 45
 myz(n) = -90;
 end
 % Plane -z
 if mzx(n) < -45 || mzx(n) > 45
 mzx(n) = -90;
 end
 if mzy(n) < -45 || mzy(n) > 45
 mzy(n) = -90;
 end

end

% Constraints on sensors due to exceeded field of view perpendicular to sensor direction
for n = 1:length(xy)
 % Plane x
 tmp = xy(n);
 if xz(n) < -55 || xz(n) > 55
 xy(n) = -90;
 end
 if tmp < -55 || tmp > 55
 xz(n) = -90;
 end
 % Plane y
 tmp = yx(n);
 if yz(n) < -55 || yz(n) > 55
 yx(n) = -90;
 end
 if tmp < -55 || tmp > 55
 yz(n) = -90;
 end
 % Plane z
 tmp = zx(n);
 if zy(n) < -55 || zy(n) > 55
 zx(n) = -90;
 end
 if tmp < -55 || tmp > 55
 zy(n) = -90;
 end
 % Plane -x
 tmp = mxy(n);
 if mxz(n) < -55 || mxz(n) > 55
 mxy(n) = -90;
 end
 if tmp < -55 || tmp > 55
 mxz(n) = -90;
 end
 % Plane -y
 tmp = myx(n);
 if myz(n) < -55 || myz(n) > 55
 myx(n) = -90;
 end
 if tmp < -55 || tmp > 55
 myz(n) = -90;
 end
 % Plane -z
 tmp = mzx(n);
 if mzy(n) <-55 || mzy(n) > 55 
    mzx(n) =-90;
 end
 if tmp <-55 || tmp > 55
    mzy(n) =-90;
 end
end

% Create plots
t = tiledlayout(2,6);

nexttile
plot(xy,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('XY');
axis([-Inf steps -45 45])
nexttile
plot(xz,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('XZ');
axis([-Inf steps -45 45])
nexttile
plot(yx,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('YX');
axis([-Inf steps -45 45])
nexttile
plot(yz,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('YZ');
axis([-Inf steps -45 45])
nexttile
plot(zx,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('ZX');
axis([-Inf steps -45 45])
nexttile
plot(zy,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('ZY');
axis([-Inf steps -45 45])
nexttile


plot(mxy,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('-XY');
axis([-Inf steps -45 45])
nexttile
plot(mxz,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('-XZ');
axis([-Inf steps -45 45])
nexttile
plot(myx,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('-YX');
axis([-Inf steps -45 45])
nexttile
plot(myz,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('-YZ');
axis([-Inf steps -45 45])
nexttile
plot(mzx,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('-ZX');
axis([-Inf steps -45 45])
nexttile
plot(mzy,'r+');
grid on
xlabel('time / [step size]');
ylabel('degrees');
title('-ZY');
axis([-Inf steps -45 45])