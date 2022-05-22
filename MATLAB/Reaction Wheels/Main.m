% Erich Bender
% AERO 463/464 -- Senior Project
close all
clear all
clc
%% Mass Moments of Inertia and Stability Chart
% Reference frame: The moments of inertia will be calculated along a body
% frame with:
% - x pointing in the RAM direction
% - y pointing in the cross-track direction to the right of the RAM
% - z pointing in the Nadir direction
Tmass = .0; % tip mass (kg)
Theight = .00; % height of cylindrical tip mass (m)
Tradius = .0; % radius of cylindrical tip mass (m)
Rmass = .00; % mass of boom (kg)
Rlength = .0; % length of boom (m)
Bheight = .300; % height of rectangular body (m)
Bwidth = .20; % width of rectangular body (m)
Bdepth = .10; % depth of rectangular body (m)
Rtip1 = [0; 0; -.5]; % vector describing position (x,y,z) of tip mass 1 from C.M. to component in m
Rtip2 = [0; 0; .5]; % vector describing position (x,y,z) of tip mass 2 from C.M. to component in m
Rrod1 = [0; 0; -.335125]; % vector describing position (x,y,z) of deployment device 1 from C.M. to component in m
Rrod2 = [0; 0; .335125]; % vector describing position (x,y,z) of deployment device 2 from C.M. to component in m
global Ix;
global Iy;
[Ix, Iy, Iz] = massmoment3U(Tmass, Tradius, Theight, Rmass, Rlength,Bheight, Bwidth, Bdepth, Rtip1, Rtip2, Rrod1, Rrod2)
global K1
K1 = ((Ix) - (Iy))/(Iz)
% K1 = -.6;
global K2
K2 = ((Iy) - (Iz))/(Ix)
% K2 = .8824;
global K3
K3 = -((K1 + K2)/(1 +K1*K2))
% Stability Chart
x1=linspace(-1,1,1000);
for i=1:length(x1)
 y1(i)=-x1(i);
end
figure(1)
plot(x1,y1, K1, K2, 'xr','MarkerSize',10)
axis([-1 1 -1 1])
grid off
hold on
line([-1,1],[0,0])
line([0,0],[-1,1])
hold off
title('Stability Chart')

xlabel('K1')
ylabel('K2')
%% Orbital Parameters
R = 6378 + 450; %km
mu = 398600; %km^3/s^2
T = ((2*pi)/sqrt(mu))*(R^(3/2)); % period in seconds
omega = sqrt(mu/(R^3)); % spin of the earth m/s
%% Initial Position (Euler Angles)
phi=0;
theta=0;
psi=0;
[e1i, e2i, e3i, e4i, a, E] = dcm2quat(phi, theta, psi);
%% Time Interval
ti = 0;
tf = 16*T; %propagate over about 24 hours
%% Inertia Properties of Gyro/Reaction Wheel
m = .180; % kg
r = .044; % m
global J
J = (m*r^2)/2;
%% Gravity Gradient Simulation
x = [.01*omega, .01*omega, omega, e1i, e2i, e3i, e4i];
options=odeset('RelTol',1e-10);
[t,y]=ode45('simulation',[ti,tf],x,options);
for i = 1:length(t)
 e1n = y(i,4);
 e2n = y(i,5);
 e3n = y(i,6);
 e4n = y(i,7);
 a(i) = sqrt(((y(i,4))^2)+((y(i,5))^2)+((y(i,6))^2)+((y(i,7)^2)));
 E = [(1 - 2*(e2n^2+e3n^2)),2*(e1n*e2n + e3n*e4n),2*(e1n*e3n - e2n*e4n);2*(e1n*e2n - e3n*e4n),(1 - 2*(e1n^2+e3n^2)),2*(e2n*e3n + e1n*e4n);2*(e1n*e3n + e2n*e4n),2*(e2n*e3n - e1n*e4n),(1 - 2*(e1n^2+e2n^2))];
 yaw(i,1)=acos(E(3,3))*180/pi;

 pitch(i,1)=acos(E(2,2))*180/pi;
 roll(i,1)=acos(E(1,1))*180/pi;
end

% Create plots
tt = tiledlayout(1,3);

nexttile

plot(t, yaw)
title('Yaw (RAM)')
xlabel('Time (s)')
ylabel('(deg)')
nexttile
plot(t, pitch)
title('Pitch (Crosstrack)')
xlabel('Time (s)')
ylabel('(deg)')
nexttile
plot(t, roll)
title('Roll (Nadir)')
xlabel('Time (s)')
ylabel('(deg)')


figure(5)
plot(t,a)
axis([0 tf .99 1.01])
title('Validity Check For Quaternions')
xlabel('Time (s)')
ylabel('Sum of the Squares')
R