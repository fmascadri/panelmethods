%MCEN90018 - Advanced Fluid Dynamics - Assignment 2 - Q4
% This script answers Question 4 in the Assignment 2 sheet.

% Author: Francesco Mascadri
% Contact: fmascadri@student.unimelb.edu.au
% April 2022

%% Clear workspace
clear all;
clc;
close all;

%% Create the boundary - just do this once then load the saved results
geometryImage = imread("cow.jpg");
geometryImage = mean(geometryImage, 3);
geometryImage = flipud(geometryImage);
% [x,y] = getGeometryPoints(geometryImage);

% best to save the x and y values as .mat files to the directory and read
% from there to save constant retracing.

%% Once cow has been traced and x and y saved from workspace, just read the .mat files
load("cow_x.mat");
load("cow_y.mat");

%% Optional plot of traced cow image
figure
hold on
ax = gca;
imagesc(ax, geometryImage);
colormap("gray")
set(ax, 'ydir','normal');

line(x,y,'Color','red','LineWidth',1);
scatter(x,y,'r','filled');
daspect([1 1 1])
axis([0 617 0 400])

title('Traced cow body with 87 panels')
xlabel('x [pixels]')
ylabel('y [pixels]')
ax.FontSize = 14;
ax.FontName = 'Arial';

hold off

%% Conversion to SI units

pixels_per_m = 600/2.5; %600px / 2.5m length of average cow

%Convert to SI units from pixel units
x = x./pixels_per_m;
y = y./pixels_per_m;

numPanels = size(x,1);
Uinf = 1; %m/s

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
[xM, yM] = meshgrid(-3:0.05:3.5, -1:0.05:2.5);
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
yS = -0.5:0.1:2.3;
xS = ones(1,length(yS)) .* (-2);

H = stream2(xM,yM,uM,vM,xS,yS,[step, maxvert]);


%% Plot velocity contour
figure
hold on

contourf(xM, yM, uM, 25, 'LineColor', 'None')

%% Plot geometry
line(x,y)
fill(x,y,'k')

%% Plot streamlines
for i = 1:size(H,2)
    streamline = cell2mat(H(i));
    plot(streamline(:,1),streamline(:,2),'-k');
end

%% Plot formatting
colormap(parula)
c = colorbar;
daspect([1 1 1]);
xlabel('x [m]')
ylabel('y [m]')
ylabel(c, 'U [m/s]');
axis([-1 3.5 -0.4 2]);
title('Top-view of adult cow in U_{\infty} = 1 m/s flow', 'FontName', 'Arial', 'FontSize', 14)
caxis([0 1.5]);
gca.FontSize = 14;
gca.FontName = 'Arial';

hold off

%% Plot actual body being modelled
figure
z = zeros(size(x));
hold on
surf([x';x'],[y';y'],[z'; z'+1])
fill3(x,y,z,'blue');
fill3(x,y,z+1,'blue');
contourf(xM, yM, uM, 25, 'LineColor', 'None')
for i = 1:size(H,2)
    streamline = cell2mat(H(i));
    plot(streamline(:,1),streamline(:,2),'-k');
end
c = colorbar;
daspect([1 1 1]);
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
ylabel(c, 'U [m/s]');
caxis([0 1.5]);
gca.FontSize = 14;
gca.FontName = 'Arial';
title('Isometric view of cow model in flow')
daspect([1 1 1])
axis([-1 3 -1 2 0 1])
view(3)

%% Functions
function [x,y] = getGeometryPoints(geometryImage)
    ax1 = axes;
    imagesc(ax1,geometryImage);
    set(ax1, 'ydir', 'normal');
    [x,y] = ginput;
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