%MCEN90018 - Advanced Fluid Dynamics - Assignment 2 - Q3
% This script answers Question 3 in the Assignment 2 sheet.

% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

%% Clear workspace
% clear all;
% clc;
% %close all;

%% Create the boundary & problem defintion
Uinf = 1;
radius = 1;
numPanels = 64;

[x,y] = generateCircleGeometry(radius, numPanels);

%% Create panels
N = numPanels;
panel = zeros(N,4); % (x_i, y_i, x_{i+1}, y_{i+1})
for i = 1:N
    panel(i,1) = x(i);
    panel(i,2) = y(i);

    if i == N % panels reconnect at the end  
        panel(i,3) = x(1);
        panel(i,4) = y(1);
    else
        panel(i,3) = x(i+1);
        panel(i,4) = y(i+1);
    end
end


%% Calculate influence coefficients
I = zeros(N);
V = zeros(N,1);

for i = 1:N
    vinf = calculateVInf(panel(i,1), panel(i,2), panel(i,3), panel(i,4), Uinf);
    V(i) = -vinf;
    for j = 1:N
        if i == j
            IElem = 0.5;
        else
            IElem = calculateIElem(panel(i,1), panel(i,2), panel(i,3), panel(i,4), panel(j,1), panel(j,2), panel(j,3), panel(j,4));
        end
        I(i,j) = IElem;
    end
end

q = I\V;

%% Calculate velocity field
[xM, yM] = meshgrid(-3:0.005:3, -3:0.005:3);
uM = ones(size(xM,1),size(xM,2)) * Uinf;
vM = zeros(size(xM,1),size(xM,2));
u = zeros(size(xM,1),size(xM,2));
v = zeros(size(xM,1),size(xM,2));

for k = 1:N
    [u,v] = calculateUV(panel(k,1),panel(k,2),panel(k,3),panel(k,4),q(k),xM,yM);
    uM = u + uM;
    vM = v + vM;
end

%% Calculate streamlines
step = 0.005;
maxvert = 300000;
yS = -2:0.25:2;
xS = ones(1,length(yS)) .* (-3);

H = stream2(xM,yM,uM,vM,xS,yS,[step, maxvert]);


%% Plot velocity contour
figure
hold on

contourf(xM, yM, uM, 25, 'LineColor', 'None')

%% Plot circle geometry
plotCircleGeometry(x,y);

%% Plot streamlines
for i = 1:size(H,2)
    streamline = cell2mat(H(i));
    plot(streamline(:,1),streamline(:,2),'-k');
end

%% Plot formatting
caxis([0 2.5])
colormap(parula)
c = colorbar;
daspect([1 1 1]);
xlabel('x [m]')
ylabel('y [m]')
ylabel(c, 'U [m/s]');
gca.FontSize = 14;
gca.FontName = 'Arial';
axis([-3 3 -2 2]);

hold off
%% Functions
function [x,y] = generateCircleGeometry(radius, numPanels)
    a = radius;
    N = numPanels;

    phi = 2*pi/N;
    phi1 = [-pi+phi/2:-phi:-3*pi+phi*3/2]';
    phi2 = [-pi-phi/2:-phi:-3*pi+phi/2]';

    x = [a*cos(phi1) a*cos(phi2)];
    y = [a*sin(phi1) a*sin(phi2)];    
end

function plotCircleGeometry(x,y,colourString)

    if nargin < 3
        colourString = 'k';
    end

    fill(x,y,colourString);
    daspect([1 1 1]);
    title(['Flow around circle geometry with ' num2str(size(x,1)) ' panels']);
    xlabel('x')
    ylabel('y')
end

function vinf = calculateVInf(Xi, Yi, Xiplus1, Yiplus1, Uinf)
    phi=atan2((Yiplus1 - Yi),(Xiplus1 - Xi)); %phi_i (eqn 24)
    vinf = Uinf*sin(2*pi - phi);
end

function IElem = calculateIElem(Xi, Yi, Xiplus1, Yiplus1, Xj, Yj, Xjplus1, Yjplus1)
    Xmi=0.5*(Xiplus1+Xi); % midpoint of panel i
    Ymi=0.5*(Yiplus1+Yi);
    Phi_i=atan2((Yiplus1 - Yi),(Xiplus1 - Xi)); %phi_i (eqn 24)

    q_j = 1;

    Xmj=0.5*(Xjplus1+Xj); % midpoint of panel j

    Ymj=0.5*(Yjplus1+Yj);
    Phi_j=atan2((Yjplus1-Yj),(Xjplus1-Xj)); %phi_j (eqn 23)

    rij = sqrt((Xmj - Xmi).^2 + (Ymj - Ymi).^2); % (eqn 22)

    beta = atan2((Ymi - Ymj),(Xmi - Xmj)); % (eqn 25)
    omega = beta - Phi_j; % (eqn 26)

    x0p = rij.*cos(omega); % (eqn 27)
    y0p = rij.*sin(omega); % (eqn 28)

    S = sqrt((Xjplus1 - Xj).^2 + (Yjplus1 - Yj).^2);
    a =-S/2;
    b =S/2; %(eqn 11)

    vj = (q_j./(2*pi)).*(atan(((S./2)-x0p)./y0p)-atan((-(S./2) - x0p)./y0p)); % eqn(30)

    uj = (q_j./(2*pi)).*((-log((y0p.^2+((S.^2)./4)- (S.*x0p)+x0p.^2))./2) + (log((y0p.^2 + ((S.^2)./4) + (S.*x0p) + x0p.^2))./2));

    vi = uj.*sin(Phi_j-Phi_i)+vj.*cos(Phi_j-Phi_i); % eqn(31)

    IElem = vi;
end

function [u,v] = calculateUV(Xi,Yi,Xiplus1,Yiplus1,q,xM,yM)
    Xmj = 0.5*(Xi + Xiplus1); % midpoints
    Ymj = 0.5*(Yi + Yiplus1);
    % equation (13)
    rij = sqrt((Xmj - xM).^2 + (Ymj - yM).^2);
    % equation (14) and (15)
    Phi = atan2((Yiplus1 - Yi),(Xiplus1 - Xi));
    beta = atan2((yM - Ymj),(xM - Xmj));
    omega = beta - Phi;
    % equations (16) and (17)
    x0p = rij.*cos(omega);
    y0p = rij.*sin(omega);
    % equation (11 & 12)
    S = sqrt((Xiplus1 - Xi).^2 + (Yiplus1 - Yi).^2);
    a = -S/2;
    b = S/2;
    % equations (18) and (19)
    vprime = (q./(2*pi)).*(atan(((S/2) - x0p)./y0p) - atan((-(S/2) - x0p)./y0p));
    uprime = (q./(2*pi)).*((-log((y0p.^2 + ((S.^2)/4) - (S.*x0p) + x0p.^2))./2) - (-log((y0p.^2 + ((S.^2)/4) + (S.*x0p) + x0p.^2))./2));
    % equations (21) and (22)
    v = vprime.*cos(Phi) + uprime.*sin(Phi);
    u = uprime.*cos(Phi) - vprime.*sin(Phi);
end