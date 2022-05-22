% Erich Bender
% AERO 463/464 -- Senior Project
%% Direction Cosine Matrix --> Quaternions
% This is a function used to develop quaternions based on given Euler
% angles phi, theta and psi, and a full 3-1-3 rotation.
function [e1, e2, e3, e4, a, E] = dcm2quat(phi, theta, psi)
% Construct direction cosine matrix based on Euler angles:
DCM1 = [cosd(phi)*cosd(psi)-sind(phi)*cosd(theta)*sind(psi),cosd(phi)*cosd(theta)*sind(psi)+sind(phi)*cosd(psi),sind(theta)*sind(psi)];
DCM2 = [-cosd(phi)*sind(psi)-sind(phi)*cosd(theta)*cosd(psi),cosd(phi)*cosd(theta)*cosd(psi)-sind(phi)*sind(psi),sind(theta)*cosd(psi)];
DCM3 = [sind(phi)*sind(theta), -cosd(phi)*sind(theta), cosd(theta)];
DCM = [DCM1;DCM2;DCM3];
% Formulae for extracting quaternions out of direction cosine matrix:
e4 = .5*sqrt(DCM(1,1)+DCM(2,2)+DCM(3,3)+1);
e3 = (DCM(1,2) - DCM(2,1))/(4*e4);
e2 = (DCM(3,1) - DCM(1,3))/(4*e4);
e1 = (DCM(2,3) - DCM(3,2))/(4*e4);
% Fundamental check on the validity of the quaternions. “a” should beequal
% to 1, or very close to it.
a = (e1)^2 + (e2)^2 + (e3)^2 + (e4)^2;
% Using the quaternion values e1, e2, e3 and e4, construct the
quaternion
% matrix, E:
E1 = [(1 - 2*(e2^2+e3^2)), 2*(e1*e2 + e3*e4), 2*(e1*e3 - e2*e4)];
E2 = [2*(e1*e2 - e3*e4),(1 - 2*(e1^2+e3^2)), 2*(e2*e3 + e1*e4)];
E3 = [2*(e1*e3 + e2*e4), 2*(e2*e3 - e1*e4), (1 - 2*(e1^2+e2^2))];
E = [E1; E2; E3];
