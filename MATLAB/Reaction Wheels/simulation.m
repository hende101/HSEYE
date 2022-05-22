% Erich Bender
% AERO 463/464 -- Senior Project
%% Non Linear Simulation
% This function utilizes ode45, quaternions and the non-linear equations of motion to
% predict angular rates of a body with supplied mass moments of inertia.
function wdot = simulation(t,x)
% Global parameters from main file:
global Ix
global Iy
global K1
global K2
global K3
global J
% Orbital Parameters
R = 6378 + 450; %km
mu = 398600; %km^3/s^2
omega = sqrt(mu/(R^3));
% Reaction Wheel rotation rate
 %sigma = 0;
sigma = (2200*2*pi)/(60); % RPM to rad/sec
E = [(1 - 2*(x(5)^2+x(6)^2)), 2*(x(4)*x(5) + x(6)*x(7)), 2*(x(4)*x(6) -(x(5)*x(7))); 2*(x(4)*x(5) - x(6)*x(7)),(1 - 2*(x(4)^2+x(6)^2)),2*(x(5)*x(6) + x(4)*x(7)); 2*(x(4)*x(6) + x(5)*x(7)), 2*(x(5)*x(6) -x(4)*x(7)), (1 - 2*(x(4)^2 + x(5)^2))];
% Equations of motion
wdot(1) = K1*(x(2)*x(3) - 3*omega^2*E(2,1)*E(3,1)) - sigma*x(2)*(J/Ix);
% w1
wdot(2) = K2*(x(1)*x(3) - 3*omega^2*E(3,1)*E(1,1)) + sigma*x(1)*(J/Iy);
% w2
wdot(3) = K3*(x(1)*x(2) - 3*omega^2*E(1,1)*E(2,1)); % w3
wdot(4) = -.5*(-(x(3) + omega)*x(5) + x(2)*x(6) - x(1)*x(7)); %e1
wdot(5) = -.5*(-x(1)*x(6) - x(2)*x(7) + (x(3) + omega)*x(4)); %e2
wdot(6) = -.5*(x(1)*x(5) - x(2)*x(4) - (x(3) - omega)*x(7)); %e3
wdot(7) = -.5*(x(1)*x(4) + x(2)*x(5) + (x(3) - omega)*x(6)); %e4
wdot = wdot';
end