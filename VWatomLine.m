%Nov 10th 2016

%Find the minimum VW Energy for an atom distance 'a' above an
%   infinite line using 'fminsearch'
%   Limiting case for atom and ring (R -> inf)

clear all
close all
clc
format long;
eps = 1;        %Depth of Potential Well
sig = 1;        %Distance for which inter-particle potential is zero
rho = 1;        %Density
aGuess = 1;

%LJ 12-6
LJ = @(d) eps*((sig./d).^12 - 2*(sig./d).^6);
%LJ evaluated at distance from (0,h) to line of atoms along x-axis
V = @(x,a) rho*LJ(sqrt(x.^2 + a.^2));
        
%Integral that evalutes total Van der Waal interaction between
%   point at (a,0) and circle with radius R
E = @(a) integral(@(x) V(x,a),-100,100);
options = optimset('TolFun',1e-8);
amin = fminsearch(E,aGuess,options)
