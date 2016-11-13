%October 25th 2016

%Plots of total enegry of a point at (a,0) and circle radius R
%   as a function of a and R

clear all
close all

eps = 1;        %Depth of Potential Well
sig = 1;        %Distance for which inter-particle potential is zero
rho = 1;        %Density
m = 100;

%Varrying a and R
R = ones(m,1)*linspace(2,10,m);
a = zeros(m,m);
for i=1:m
    a(:,i)=linspace(0,R(1,i)-0.9,m)';
end
E = zeros(m,m);
Emin = zeros(1,m);
amin = zeros(1,m);


for i = 1:m
    for j=1:m
        LJ = @(d) eps*rho*((sig./d).^12 - 2*(sig./d).^6);
        
        %LJ evaluated at distance from (a,0) to R(cos(theta),sin(theta))
        V = @(s) LJ(sqrt((a(j,i)-R(j,i).*cos(s/R(j,i))).^2 ...
            + R(j,i).^2.*sin(s/R(j,i)).^2));
        
        %Integral that evalutes total Van der Waal interaction between
        %   point at (a,0) and circle with radius R
        E(j,i) = integral(V,0,2*pi*R(j,i));
    end
%     figure(1)
%     hold on
%     plot(a(:,i),E(:,i))
%     axis([0 1.1 -6.5 3]);
%     tl=sprintf('Radius=%5.3f',R(1,i));
%     title(tl);
%     xlabel('a');
%     ylabel('E');
%     %pause(0.25)

    [Em,k]=min(E(:,i));
    Emin(i)=E(k,i);
    amin(i)=a(k,i);
end


% figure(2)
% surf(R,a,E)
% colormap winter
% xlabel('R');
% ylabel('a');
% zlabel('E');


figure(3)
plot(R(1,:),Emin);
title('Minimum Energy');
xlabel('R');
ylabel('E_{min}');

figure(4)
hold on
plot(R(1,:),R(1,:)-amin);
title('Value of R-a_{min} corresponding to the energy minimum');
xlabel('R');
ylabel('R-a_{min}');

%save('resultsLimitingTest.mat')
