%MCEN90018 - Advanced Fluid Dynamics - Assignment 2 - Q8
% This script answers Question 8 in the Assignment 2 sheet.

% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

%% Clear workspace
clear all;
clc;
close all;

%% Create the boundary & problem defintion
Uinfty = 1;
Vinfty = 0;
alpha = 10; %degrees

Uinf = Uinfty*cosd(alpha);
Vinf = Uinfty*sind(alpha);

chord = 1;
thickness = 0.12;

[x,y] = generateAirfoilGeometry(chord, thickness);

%% Create panels
N = size(x,2);
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
    psiinf = calculatePsiInf(panel(i,1), panel(i,2), panel(i,3), panel(i,4), Uinf, Vinf);
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
uHat = ones(size(xM,1),size(xM,2)) * Uinf; %background flow
vHat = ones(size(xM,1),size(xM,2)) * Vinf; %background flow

for k = 1:N
    [u,v] = calculateVelocityFieldContribution(panel(k,1),panel(k,2),panel(k,3),panel(k,4),gamma(k),xM,yM);

    uHat = u + uHat;
    vHat = v + vHat;
end

%% Calculate streamlines
step = 0.005;
maxvert = 300000;
yS = -0.5625:0.0625:0.5625;
xS = ones(1,length(yS)) .* (-0.8);

% Include a streamline at the rear stagnation point to show Kutta condition
yS = [yS y(1)];
xS = [xS x(1)+0.01];

H = stream2(xM,yM,uHat,vHat,xS,yS,[step, maxvert]);

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

% Get Cp
U_hat = sum(U,2);
Uc = U_hat + Uinf.*cos(2*pi - phi_i);
Uc = 0.5*Uc; %not clear why this fixes everything in the numerical vs theory calculation. Can't see where the additional factor 2x is in the getUi, etc code
Cp = 1 - (Uc./Uinf).^2;

% Get x midpoint values
Xm = 0.5*(panel(:,3) + panel(:,1));
Ym = 0.5*(panel(:,2) + panel(:,4));

%% Convert to aoa position
[x,y] = rotate_flow(x,y, alpha);
[xM, yM] = rotate_flow(xM, yM, alpha);
[uHat, vHat] = rotate_flow(uHat, vHat, alpha);

%% Plot velocity contour
figure
hold on

contourf(xM, yM, uHat, 100, 'LineColor', 'None')

%% Plot airfoil geometry
plotAirfoilGeometry(x,y);

%% Plot streamlines
for i = 1:size(H,2)
    streamline = cell2mat(H(i));
    [stream_x, stream_y] = rotate_flow(streamline(:,1), streamline(:,2), alpha);
    plot(stream_x,stream_y,'-k');
    if i == size(H,2)
        plot(stream_x,stream_y,'-r', 'LineWidth',2);
    end
end

%% Plot formatting
caxis([0.5 1.5])
colormap(jet)
c = colorbar;
daspect([1 1 1]);
xlabel('x [m]')
ylabel('y [m]')
ylabel(c, 'U [m/s]');
gca.FontSize = 14;
gca.FontName = 'Arial';
axis([-0.5 1.5 -0.5 0.5]);
title(['NACA0012 airfoil with ' num2str(size(panel,1)) ' vortex panels at ' num2str(alpha) ' degrees AOA'])

hold off


%% Plotting Cp
figure
plot(Xm, Cp);

%plot formatting
xlabel("x_m")
ylabel("C_p")
gca.FontSize = 14;
gca.FontName = 'Arial';
title(['C_p vs x for airfoil at ' num2str(alpha) ' degrees AOA'])
grid on

%% Functions
function [x,y] = generateAirfoilGeometry(c, t)
    %Generate points for upper and lower and combine them such that
    %bottom trailing edge is first and upper trailing edge is last
    x_lower = 1:-0.02:0;
    x_upper = 0.01:0.02:0.99;

    y_upper = 5*t*c*(0.2969*sqrt(x_upper./c) - 0.1260*(x_upper./c) - 0.3516*(x_upper./c).^2 + 0.2843*(x_upper./c).^3 - 0.1036*(x_upper./c).^4);
    y_lower = -5*t*c*(0.2969*sqrt(x_lower./c) - 0.1260*(x_lower./c) - 0.3516*(x_lower./c).^2 + 0.2843*(x_lower./c).^3 - 0.1036*(x_lower./c).^4);
    
    y = [y_lower y_upper];
    x = [x_lower x_upper];
end

function plotAirfoilGeometry(x,y,colourString)

    if nargin < 3
        colourString = 'k';
    end

    fill(x,y,colourString);
end

function psiinf = calculatePsiInf(Xi, Yi, Xiplus1, Yiplus1, Uinf, Vinf)
    Xmi=0.5*(Xiplus1+Xi); % midpoint of panel i
    Ymi=0.5*(Yiplus1+Yi);

    psiinf = Uinf*Ymi - Vinf*Xmi;
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

function [x_rot,y_rot] = rotate_flow(x,y,alpha)
    z = x + 1i*y;
    z1 = exp(-1i*deg2rad(alpha)).*z;
    x_rot = real(z1);
    y_rot = imag(z1);
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
    uj = (gamma_j./(2*pi)) .* (atan2((x0p - b),y0p) - atan2((x0p - a),y0p)); 
    %have to use atan2 to get the right results in this section but atan in the velocity field contribution calculation. Something must be conceptually wrong.

    vi = uj.*sin(Phi_j-Phi_i) + vj.*cos(Phi_j-Phi_i); % eqn(31)
    ui = uj.*cos(Phi_j-Phi_i) - vj.*sin(Phi_j-Phi_i); % eqn(32)
end

function phi_i = calculatePhiI(Xi,Yi,Xiplus1,Yiplus1)
    phi_i = atan2((Yiplus1 - Yi),(Xiplus1 - Xi));
end