%MCEN90018 - Advanced Fluid Dynamics - Assignment 2 - Q3
% This script answers Question 3 in the Assignment 2 sheet.

% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

%% Clear workspace
clear all;
clc;
close all;

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

%% Calculate numerical Cp
U = zeros(N);
phi_i = zeros(N,1);

% Get ui
for i = 1:N
    phi_i(i) = calculatePhiI(panel(i,1),panel(i,2),panel(i,3),panel(i,4));
    for j = 1:N
        ui = calculateUi(panel(i,1),panel(i,2),panel(i,3),panel(i,4),panel(j,1),panel(j,2),panel(j,3),panel(j,4),q(j));
        U(i,j) = ui;
    end
end

% Get Cp
U_hat = sum(U,2);
Uc = U_hat + Uinf.*cos(2*pi - phi_i);
Cp = 1 - (Uc./Uinf).^2;

% Get theta values
Xm = 0.5*(panel(:,3) + panel(:,1));
Ym = 0.5*(panel(:,4) + panel(:,2));
theta = atan2(Ym,Xm);

% convert atan2 [-pi,pi] range to [0,2pi] for theta
for i = 1:size(theta)
    if theta(i) < 0
        theta(i) = theta(i) + 2*pi;
    end
end

%% Calculate theorectical solution
theta_th = 0:0.1:2*pi;
Cp_th = 1 - 4.*sin(theta_th).^2;

%% Plotting
figure
hold on
plot(theta_th, Cp_th);
scatter(theta, Cp, 'filled');
hold off

%plot formatting
xlabel("\theta [rad]")
ylabel("C_p")
gca.FontSize = 14;
gca.FontName = 'Arial';
legend("Theoretical", "Numerical")
title(['Comparison of theoretical vs numerical C_p for ' num2str(N) '-panel geometry'])
grid on
axis([0 2*pi -3 1]);

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
    title(['Panelised circle geometry with ' num2str(size(x,1)) ' panels']);
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

function ui = calculateUi(Xi, Yi, Xiplus1, Yiplus1, Xj, Yj, Xjplus1, Yjplus1, q_j)
    Xmi=0.5*(Xiplus1+Xi); % midpoint of panel i
    Ymi=0.5*(Yiplus1+Yi);
    Phi_i=atan2((Yiplus1 - Yi),(Xiplus1 - Xi)); %phi_i (eqn 24)

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

    vj = (q_j./(2*pi)).*(atan(((S./2)-x0p)./y0p)...
        - atan((-(S./2) - x0p)./y0p)); % eqn(30)

    uj = (q_j./(2*pi)).*((-log((y0p.^2+((S.^2)./4)- (S.*x0p)+x0p.^2))./2)...
         + (log((y0p.^2 + ((S.^2)./4) + (S.*x0p) + x0p.^2))./2)); % eqn(29)

    vi = uj.*sin(Phi_j-Phi_i)+vj.*cos(Phi_j-Phi_i); % eqn(31)

    ui = uj.*cos(Phi_j-Phi_i) - vj.*sin(Phi_j - Phi_i); % eqn(32)
end

function phi_i = calculatePhiI(Xi,Yi,Xiplus1,Yiplus1)
    phi_i = atan2((Yiplus1 - Yi),(Xiplus1 - Xi));
end