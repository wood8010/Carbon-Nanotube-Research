%October 25th 2016

%Plot surface of total enegry of a point at (0,h) and the line of atoms
%   along the x-axis

clear all
close all
clc

eps = 1;        %Depth of Potential Well
sig = 1;        %Distance for which inter-particle potential is zero
rho = 1;        %Density
m = 1000;

%Varrying h
a = linspace(.8,1.2,m);
E = zeros(1,m);

for i = 1:m
    %LJ 12-6
    LJ = @(d) eps*((sig./d).^12 - 2*(sig./d).^6);
    %LJ evaluated at distance from (0,h) to line of atoms along x-axis
    V = @(x) rho*LJ(sqrt(x.^2 + a(i).^2));
        
    %Integral that evalutes total Van der Waal interaction between
    %   point at (a,0) and circle with radius R
    E(i) = integral(V,-100,100);
    
end

[Em,k] = min(E);
amin = a(k)

figure(1)
hold on
plot(a(:),E(:))
axis([.75 1.5 -1.9 2.5]);
title('Total V-W Energy E as a Function of a');
xlabel('h');
ylabel('E');

figure(2)
hold on
plot(a(:),E(:))
axis([.9 1 -1.725 -1.5]);
title('Total V-W Energy E as a Function of a');
xlabel('h');
ylabel('E');