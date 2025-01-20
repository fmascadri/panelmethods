%MCEN90018 - Advanced Fluid Dynamics - Assignment 2 - Q6
% This script answers Question 6 in the Assignment 2 sheet.

% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

%% Clear workspace
clear all;
clc;
close all;

%% Create the boundary & problem defintion
Uinf = 1;
Vinf = 0;
radius = 1;
numPanels = 64;
Gamma = 2*Uinf*pi*radius;

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
PsiInf = zeros(N,1);

for i = 1:N
    psiinf = calculatePsiInf(panel(i,1), panel(i,2), panel(i,3), panel(i,4), Uinf, Vinf, Gamma);
    PsiInf(i) = -psiinf;
    for j = 1:N
        IElem = calculateIElem(panel(i,1), panel(i,2), panel(i,3), panel(i,4), panel(j,1), panel(j,2), panel(j,3), panel(j,4));
        I(i,j) = IElem;
    end
end

% Modify matrices to correctly constrain the vortex panel method
n = size(I,1);
I_i_Nplus1 = ones(n,1); %last column of I matrix with rows equal to N. Final N+1 row taken care of by next modification.
I_Nplus1_i = zeros(1,n+1);
I_Nplus1_i(1) = 1; %overwrite zero (as per notes matrix eq 11, sets gamma_1 == gamma_N to satisfy Kutta condition)
I_Nplus1_i(n) = 1; %overwrite zero (as per notes matrix eq 11)

% Assemble final I matrix
I = [I I_i_Nplus1];
I = [I; I_Nplus1_i];

PsiInf = [PsiInf; 0]; %add 0 to final row of PsiInf vector

% Result
gamma = I\PsiInf;

%% Calculate velocity field
[xM, yM] = meshgrid(-3:0.005:3, -3:0.005:3);
uHat = ones(size(xM,1),size(xM,2)) * Uinf; % free stream background flow
vHat = ones(size(xM,1),size(xM,2)) * Vinf; % free stream background flow

%Flow due to rotation
r = sqrt(xM.^2 + yM.^2);
u_theta = Gamma./(2.*pi.*r);

theta = atan2(yM,xM);
u_rot = u_theta.*sin(theta);
v_rot = -u_theta.*cos(theta);

uHat = uHat + u_rot;
vHat = vHat + v_rot;

%Flow due to influence of panels
for k = 1:N
    [u,v] = calculateVelocityFieldContribution(panel(k,1),panel(k,2),panel(k,3),panel(k,4),gamma(k),xM,yM);

    uHat = u + uHat;
    vHat = v + vHat;
end

%% Calculate streamlines
step = 0.005;
maxvert = 300000;
yS = -2:0.25:2;
xS = ones(1,length(yS)) .* (-3);

H = stream2(xM,yM,uHat,vHat,xS,yS,[step, maxvert]);

%% Plot velocity contour
% Filter singularity at center of vortex flow
idx = find(r < 0.8);
uHat(idx) = NaN;

figure
hold on

contourf(xM, yM, uHat, 100, 'LineColor', 'None')

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

%% Calculate numerical Cp
U = zeros(N);
phi_i = zeros(N,1);

% Get ui
for i = 1:N
    phi_i(i) = calculatePhiI(panel(i,1),panel(i,2),panel(i,3),panel(i,4));
    for j = 1:N
        ui = calculateUi(panel(i,1),panel(i,2),panel(i,3),panel(i,4),panel(j,1),panel(j,2),panel(j,3),panel(j,4),gamma(j));
        U(i,j) = ui;
    end
end
U_hat = sum(U,2);

% U contribution from uniform flow
U_uniform = Uinf.*cos(2*pi - phi_i);

% U contribution from rotation
U_rotation = (Gamma/(2*pi)); %(radius = 1 since on the surface)

% Get Cp
U_partial = (U_hat + U_uniform)*0.5; %still don't know why 0.5x factor is required to make it work...
Uc = U_partial + U_rotation;
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

%% Calculate theorectical solution - rotating cylinder
theta_th = 0:0.1:2*pi;
Cp_th = 1 - 4.*sin(theta_th).^2 - (2*Gamma/(Uinf*pi*radius))*sin(theta_th) - (Gamma/(2*Uinf*pi*radius))^2;

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
legend("Theoretical", "Numerical", 'Location','southeast');
title("Comparison of theoretical vs numerical C_p for rotating cylinder")
grid on
axis([0 2*pi -8 1]);

%% Functions
function [x,y] = generateCircleGeometry(radius, numPanels)
    a = radius;
    N = numPanels;

    phi = 2*pi/N;
    %Changed panel creation to start with the first node at 0 degrees so
    %that panel 1 and panel N satisfy the Kutta condition at the node 
    %where they connect (at theta=0deg):
    phi1 = [0:-phi:-2*pi+phi]';
    phi2 = [-phi:-phi:-2*pi]';

    x = [a*cos(phi1) a*cos(phi2)];
    y = [a*sin(phi1) a*sin(phi2)];    
end

function plotCircleGeometry(x,y,colourString)

    if nargin < 3
        colourString = 'k';
    end

    fill(x,y,colourString);
    daspect([1 1 1]);
    title(['Flow around rotating cylinder with 64 vortex panels and \Gamma = 2U_{\infty}\pi a']);
    xlabel('x')
    ylabel('y')
end

function psiinf = calculatePsiInf(Xi, Yi, Xiplus1, Yiplus1, Uinf, Vinf, Gamma)
    Xmi=0.5*(Xiplus1+Xi); % midpoint of panel i
    Ymi=0.5*(Yiplus1+Yi);

    % Addition of background vortex for rotating cylinder
    r = sqrt(Xmi^2 + Ymi^2);
    psi_vortex = -(Gamma/2*pi)*log(r);

    psiinf = Uinf*Ymi - Vinf*Xmi + psi_vortex;
end

function IElem = calculateIElem(Xi, Yi, Xiplus1, Yiplus1, Xj, Yj, Xjplus1, Yjplus1)
    Xmi=0.5*(Xiplus1+Xi); % midpoint of panel i
    Ymi=0.5*(Yiplus1+Yi);
    Phi_i=atan2((Yiplus1 - Yi),(Xiplus1 - Xi)); %phi_i (eqn 24)

    Xmj=0.5*(Xjplus1+Xj); % midpoint of panel j
    Ymj=0.5*(Yjplus1+Yj);
    Phi_j=atan2((Yjplus1-Yj),(Xjplus1-Xj)); %phi_j (eqn 23)

    rij = sqrt((Xmi - Xmj).^2 + (Ymi - Ymj).^2); % (eqn 22)

    beta = atan2((Ymi - Ymj),(Xmi - Xmj)); % (eqn 25)
    omega = beta - Phi_j; % (eqn 26)

    x0p = rij.*cos(omega); % (eqn 27)
    y0p = rij.*sin(omega); % (eqn 28)

    S = sqrt((Xjplus1 - Xj).^2 + (Yjplus1 - Yj).^2);
    a = -S/2;
    b = S/2; %(eqn 11)

    % Vortex panel modification %(eqn 9)
    psi = (1/(2*pi))*( ((x0p-b)/2) * log((x0p - b).^2 + y0p.^2) + y0p*atan((x0p - b)./y0p) + b) ...
        - (1/(2*pi))*( ((x0p-a)/2) * log((x0p - a).^2 + y0p.^2) + y0p*atan((x0p - a)./y0p) + a);

    IElem = psi;
end

function [u,v] = calculateVelocityFieldContribution(Xj,Yj,Xjplus1,Yjplus1,gamma,xM,yM)
    Xmj = 0.5*(Xj + Xjplus1); % midpoints
    Ymj = 0.5*(Yj + Yjplus1);
    % equation (13)
    r = sqrt((xM - Xmj).^2 + (yM - Ymj).^2);
    % equation (14) and (15)
    Phi = atan2((Yjplus1 - Yj),(Xjplus1 - Xj));
    beta = atan2((yM - Ymj),(xM - Xmj));
    omega = beta - Phi;
    % equations (16) and (17)
    x0p = r.*cos(omega);
    y0p = r.*sin(omega);
    % equation (11 & 12)
    S = sqrt((Xjplus1 - Xj).^2 + (Yjplus1 - Yj).^2);
    a = -S/2;
    b = S/2;

    %Panel frame
    up = calculateUVelocity(gamma, x0p, y0p, a, b);
    vp = calculateVVelocity(gamma, x0p, y0p, a, b);

    %Convert to global frame
    u = up.*cos(Phi) - vp.*sin(Phi);
    v = up.*sin(Phi) + vp.*cos(Phi);
end

function u = calculateUVelocity(gamma, x0p, y0p, a, b)
% (eqn 6 notes)
    u = (gamma./(2*pi)) .* (atan((x0p - b)./y0p) - atan((x0p - a)./y0p));
end

function v = calculateVVelocity(gamma, x0p, y0p, a, b)
% (eqn 7 notes)
    v = (-gamma./(2*pi)) .* ( 0.5* ( log( (x0p-b).^2 + y0p.^2) ) - 0.5*( log( (x0p-a).^2 + y0p.^2) ) );
end

function ui = calculateUi(Xi, Yi, Xiplus1, Yiplus1, Xj, Yj, Xjplus1, Yjplus1, gamma_j)
    Xmi=0.5*(Xiplus1+Xi); % midpoint of panel i
    Ymi=0.5*(Yiplus1+Yi);
    Phi_i=atan2((Yiplus1 - Yi),(Xiplus1 - Xi)); %phi_i (eqn 24)

    Xmj=0.5*(Xjplus1+Xj); % midpoint of panel j
    Ymj=0.5*(Yjplus1+Yj);
    Phi_j=atan2((Yjplus1-Yj),(Xjplus1-Xj)); %phi_j (eqn 23)

    rij = sqrt((Xmi - Xmj).^2 + (Ymi - Ymj).^2);

    beta = atan2((Ymi - Ymj),(Xmi - Xmj)); % (eqn 25)
    omega = beta - Phi_j; % (eqn 26)

    x0p = rij.*cos(omega); % (eqn 27)   
    y0p = rij.*sin(omega); % (eqn 28)

    S = sqrt((Xjplus1 - Xj).^2 + (Yjplus1 - Yj).^2);
    a =-S/2;
    b =S/2; %(eqn 11)

    vj = calculateVVelocity(gamma_j, x0p, y0p, a, b);
    uj = (gamma_j./(2*pi)) .* (atan2((x0p - b),y0p) - atan2((x0p - a),y0p)); %have to use atan2 to get the right results in this section but atan in the velocity field contribution calculation. Something must be conceptually wrong.

    vi = uj.*sin(Phi_j-Phi_i) + vj.*cos(Phi_j-Phi_i); % eqn(31)
    ui = uj.*cos(Phi_j-Phi_i) - vj.*sin(Phi_j-Phi_i); % eqn(32)
end

function phi_i = calculatePhiI(Xi,Yi,Xiplus1,Yiplus1)
    phi_i = atan2((Yiplus1 - Yi),(Xiplus1 - Xi));
end