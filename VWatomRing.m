%Nov 3rd 2016

%Plots R-amin to show it approaches .95 in the limiting case (R Large)
%   Defining energy as an integral function w.r.t. 'a' 
%   using 'fminsearch'

clear all
close all

eps = 1;        %Depth of Potential Well
sig = 1;        %Distance for which inter-particle potential is zero
rho = 1;        %Density
m = 100;        %number of R-values
Rleft = 100;
Rright = 200;
R = linspace(Rleft,Rright,m);   %R from 8 to 10 with 100 values
aGuess = Rleft-0.95;        %scale for number of a-values wrt largest R
deltaR = (Rright-Rleft)/(m-1);
Emin = zeros(1,m);
amin = zeros(1,m);

for i = 1:m
    LJ = @(d) eps*((sig./d).^12 - 2*(sig./d).^6);
        
    %LJ evaluated at distance from (a,0) to R(cos(theta),sin(theta))
    V = @(s,a) rho*LJ(sqrt((a-R(i).*cos(s/R(i))).^2 ...
            + R(i).^2.*sin(s/R(i)).^2));
        
    %Integral that evalutes total Van der Waal interaction between
    %   point at (a,0) and circle with radius R
    E = @(a) integral(@(s) V(s,a),0,2*pi*R(i));
    options = optimset('TolFun',1e-8);
    [amin(i),Emin(i)] = fminsearch(E,aGuess,options);
    aGuess = amin(i)+deltaR;
end

% figure(1)
% plot(R,Emin);
% title('Minimum Energy');
% xlabel('R');
% ylabel('E_{min}');

figure(2)
plot(R,R-amin);
title('Value of R-a_{min} corresponding to the energy minimum');
xlabel('R');
ylabel('R-a_{min}');

% figure(3)
% plot(R,amin);
% title('Value of a_{min} corresponding to the energy minimum');
% xlabel('R');
% ylabel('a_{min}');

%save('resutlsLimitingTest.mat')