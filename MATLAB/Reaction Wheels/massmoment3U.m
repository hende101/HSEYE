% Erich Bender
% AERO 463/464 -- Senior Project
%% Mass Moment of Inertia
% This function roughly approximates the mass moment of inertia for a 3U
% CubeSat.
function [Ix, Iy, Iz] = massmoment3U(Tmass, Tradius, Theight, Rmass, Rlength, Bheight, Bwidth, Bdepth, Rtip1, Rtip2, Rrod1, Rrod2)
% All masses in kilograms. Heights, widths, depths and radii are in meters.
%Inertia tensor for tip mass:
Itip = [1/12*Tmass*(3*(Tradius)^2+(Theight)^2), 0, 0; 0, 1/2*Tmass*(Tradius)^2, 0; 0, 0, 1/12*Tmass*(3*(Tradius)^2+(Theight)^2)];
%Inertia tensor for deployment rods:
Irod = [1/12*Rmass*(Rlength)^2, 0, 0; 0, 1/12*Rmass*(Rlength)^2, 0; 0, 0, 0];
%Inertia tensor for main body (3U in our case):
Bmass = 12-2*Rmass-2*Tmass; % kilograms
Ibody = [1/12*Bmass*((Bheight)^2+(Bdepth)^2), 0, 0; 0, 1/12*Bmass*((Bwidth)^2+(Bheight)^2), 0; 0, 0, 1/12*Bmass*((Bwidth)^2+(Bdepth)^2)];
% Placing gravity gradient booms and tip masses properly with respect to the
% center of mass. Current gravity gradient configuration is tip masses 1
% meter apart:
% Parallel axis theorem:
Itip1_offset = Itip + Tmass*(dot(Rtip1,Rtip1)*eye(3) - Rtip1*Rtip1');
Itip2_offset = Itip + Tmass*(dot(Rtip2,Rtip2)*eye(3) - Rtip2*Rtip2');
Irod1_offset = Irod + Rmass*(dot(Rrod1,Rrod1)*eye(3) - Rrod1*Rrod1');
Irod2_offset = Irod + Rmass*(dot(Rrod2,Rrod2)*eye(3) - Rrod2*Rrod2');
% Total inertia tensor:
Itot = Itip1_offset + Itip2_offset + Irod1_offset + Irod2_offset + Ibody;
Ix = Itot(1,1);
Iy = Itot(2,2);
Iz = Itot(3,3);
end
