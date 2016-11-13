%Computing VW interaction for concentric ellipses
%   where ellipse1 = <A*cos(x),B*sin(x)> and
%   ellipse2 = <C*cos(x),D*sin(x)> arbitrarily

close all 

sig = 1;
eps = 1;
row = 1;
L = 2*pi;
%Parameters for ellipse1
A = row;
B = row;
%Parameters for ellipse2
alpha = 1;
beta = 1;


for i=1:.01:3
    R = .7+i;
    C = alpha*R;
    D = beta*R;

    %This 'length' just comes from the distance from one ellipse to the
    %other
    length = @(x,y) sqrt((A*cos(x)-C*cos(y)).^2 + (B*sin(x)-D*sin(y)).^2);
      
    LJ = @(d) eps*((sig./d).^12-2*(sig./d).^6);
    V = @(s1,s2) LJ(length(s1,s2));
    
    Q = integral2(V,0,L,0,L);
  
    figure(2)
    hold on
    title('Van der Waal Interaction for Two Concentric Ellipses')
    plot(R,Q,'o')
end