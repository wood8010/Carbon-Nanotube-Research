%September 8th 2016

%Compare difference 
%Van der Waals interaction between
%   c1 and c2, circles of radii
%   row and R respectively

close all

%Paramenters
eps = 1;
sig = 1;
rho = 1;
r = 10^(-6);
Qmin = 1000;

%Computing VW interation for concentric circles
%   with one of the radii changing from .7 to 3.7 in steps of .01
%   and fixing the other radius at 1
for R=r+.7:.01:r+2
    
    %I don't know what we should be using for theta
    %I got this from evaluating bending energy from early on
    %theta = @(s) 2*pi/L*s;
    
    %This 'length' just comes from the distance from one cirlce to the
    %other
    length = @(x,y) sqrt(r^2 + R^2 - ...
        2*r*R.*cos(x-y));
      
    LJ = @(d) eps*rho*((sig./d).^12-2*(sig./d).^6);
    V = @(theta1,theta2) LJ(length(theta1,theta2));
    
    Q = r*R*integral2(V,0,2*pi,0,2*pi);
    
    if Q < Qmin
        Qmin = Q;
        Rmin = R;
    end
    
%     figure(1)
%     grid on
%     title('Van der Waal Interaction for Two Concentric Circles')
%     hold on
%     plot(R-r,Q,'o')
end
Rmin
Rmin-r


%Computing VW interaction for concentric ellipses
%   where ellipse1 = <A*cos(x),B*sin(x)> and
%   ellipse2 = <C*cos(x),D*sin(x)> arbitrarily

% %Parameters for ellipse1
% A = 2;
% B = 1;
% %Parameters for ellipse2
% alpha = 2;
% beta = 1;


% for i=1:.01:3
%     R = .7+i;
%     C = alpha*R;
%     D = beta*R;
% 
%     %This 'length' just comes from the distance from one ellipse to the
%     %other
%     length = @(x,y) sqrt((A*cos(x)-C*cos(y)).^2 + (B*sin(x)-D*sin(y)).^2);
%       
%     LJ = @(d) eps*((sig./d).^12-2*(sig./d).^6);
%     V = @(s1,s2) LJ(length(s1,s2));
%     
%   %Find J for change of variables for this case
%     Q = integral2(V,0,2*pi,0,2*pi);
%   
%     figure(2)
%     hold on
%     title('Van der Waal Interaction for Two Concentric Ellipses')
%     plot(R,Q,'o')
% end


