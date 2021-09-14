%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code to determine the Integral Feasible Prestress of Tensegrity 
% Structures using the Doble Singular Value Decomposition (DSVD) method 
% developed by: [Yuan, X., Chen, L., & Dong, S. (2007). Prestress design 
% of cable domes with new forms. International Journal of Solids and 
% Structures, 44(9), 2773-2782.]
% (https://doi.org/10.1016/j.ijsolstr.2006.08.026)

	
% By: 
%   Jaswant Cobos
%   jaswant.cobos@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close
clear 
clc

%% DATA

% This code read the information from an Excel document (DOC.xlsx)
% The user should order the sheets from the document this way:

% Sheet 1: Connectivity matrix (CON) (don't insert headers)
%      |Element|Initial_Node|Final_Node|Symmetry_Group|Bar_Cable|
%   Symmetry_groups must be numbers
%   The Bar_cable column should be filled with 0 if the member is a bar
%   and with 1 when the member is a cable

% Sheet 2: Coordinate matrix (COOR) (don't insert headers)
%       |Node|x_Coordinates|y_Coordinates|y_Coordinates|
%   Nodes must be numbers

% Sheet 3: Free nodes vector (FN) (don't insert headers)
%   Free nodes must be numbers
%   Insert the number of the free nodes in the first column of this sheet

CON = xlsread('GEIGER_DOME', 1); % Connectivity matrix
COOR = xlsread('GEIGER_DOME', 2); % Coordinate matrix
FN = xlsread('GEIGER_DOME', 3); % Fre nodes vector
SN = 'GEIGER_DOME'; % Name of the structure

% Tolerances for the matrix rank during the two 
% Singular Value Decomposition 
tol_1 = 0.00001; % Set the tolerance for the first SVD
tol_2 = 0.0001; % Set the tolerance for the second SVD

% Graphical parameters
addnodenumber = false; % true = add the node numbers to the final plot
addfixednodes = true; % true = add the node fixed symbol in the final plot

%% PREVIOUS CALCULATIONS:

b = size(CON, 1); % Number of members
CON = sortrows(CON, (4)); % Ordering by symmetry group
sg = unique(CON(:,4),'stable');% Symmetry groups
n = size(sg, 1); % Number of symmetry groups
fn = size(FN, 1); % Number of free nodes
SN = strrep(SN, '_', ' '); % Changing the name of the structure

%% GRAPH OF THE ELEMENTS

close

% The next graph helps to visualize if all the memebers are connected
% properly or if the node coordinates are well defined

% Settings for the graph:

%   dx, dy, and dz represent the 10% of the distance between the maximum
%   and minimum coordinate for each axis
dx = (max(COOR(:, 2)) - min(COOR(:, 2))) / 10;
dy = (max(COOR(:, 3)) - min(COOR(:, 3))) / 10;
dz = (max(COOR(:, 4)) - min(COOR(:, 4))) / 10;

%   x, y, and z max and min, are the maximun and minimum values for the
%   graph window
xmin = min(COOR(:, 2)) - dx;
xmax = max(COOR(:, 2)) + dx;

ymin = min(COOR(:, 3)) - dy;
ymax = max(COOR(:, 3)) + dy;

zmin = min(COOR(:, 4)) - dz;
zmax = max(COOR(:, 4)) + dz;

%   Headers, labels, and legend
title(SN,'fontsize', 16, 'linewidth', 0.7)

text(xmax, ymax + dy, (zmax / 2), ...
    'CABLES', 'color', 'w', 'backgroundcolor', [190 190 190]/255, ...
    'linewidth', 1);
text(xmax, ymax, ((zmax / 2) + dz), ...
    'BARS', 'color', 'w', 'backgroundcolor', 'k', 'linewidth', 1);

xlabel('X')
ylabel('Y')
zlabel('Z')

% Plot of the members
for i = 1:b
    % Initial and final node coordinates for each element
    x = [COOR(CON(i, 2), 2), COOR(CON(i, 3), 2)];
    y = [COOR(CON(i, 2), 3), COOR(CON(i, 3), 3)];
    z = [COOR(CON(i, 2), 4), COOR(CON(i, 3), 4)];
    
    hold on;
    
    if CON(i, 5) == 1
        plot3(x, y, z, 'color', [190 190 190]/255, 'linewidth', 1);
    else
        plot3(x, y, z, 'k', 'linewidth', 1);
    end
end

% Number of the nodes
text(COOR(:, 2), COOR(:, 3), COOR(:, 4), string(COOR(:, 1)));

axis equal;

if sum(abs(COOR(:, 4))) == 0
    axis([xmin xmax ymin ymax]) % Axis limits configuration
else
    axis([xmin xmax ymin ymax zmin zmax]) % Axis limits configuration
    view([xmax ymin zmax]) % Point of view configuration
end

set(gcf,'WindowState', 'maximized')

%% CONNECTIVITY MATRIX ASSEMBLY

% This matrix was assembled using the method of the chapter 2 of the
% book "Tensegrity Structures Form, Stability, and Symmetry" written by
% Jin Yao Zhang and Makoto Ohsaki

% Connectivity matrix used to find the equilibrium matrix:
Cs = zeros(b, size(COOR, 1));
for i = 1:b
    Cs(i, min(CON(i, 2:3))) = 1;
    Cs(i, max(CON(i, 2:3))) = -1;
end

%% LENGTH OF THE ELEMENTS

% These matrices were assembled using the method of the chapter 2 of the
% book "Tensegrity Structures Form, Stability, and Symmetry" written by
% Jin Yao Zhang and Makoto Ohsaki

u = Cs * COOR(:, 2); % X coordinate difference vector
v = Cs * COOR(:, 3); % Y coordinate difference vector
w = Cs * COOR(:, 4); % Z coordinate difference vector

UU = diag(u); % Matrix with u as diagonal
VV = diag(v); % Matrix with v as diagonal
WW = diag(w); % Matrix with w as diagonal

l = sqrt((UU * u) + (VV * v) + (WW * w)); % Length of the elements
L = sqrt(UU ^ 2 + VV ^ 2 + WW ^ 2); % Matrix with l as diagonal

%% EQUILIBRIUM MATRIX ASSEMBLY

% This matrix was assembled using the method of the chapter 3 of the
% book "Computational Modeling of Tensegrity Structures Art, Nature,
% Mechanical and Biological Systems" written by Buntara Sthenly Gan

% Equilibrium matrix calculation:
Ab=[(Cs' * diag(Cs * COOR(:, 2) ./ l)); ... % Complete equilibrium matrix
    (Cs' * diag(Cs * COOR(:, 3) ./ l)); ... % it includes free nodes and
    (Cs' * diag(Cs * COOR(:, 4) ./ l))];    % fixed nodes              

% Extraction of the equilibrium matrix for the free nodes
Vefn = ismember(COOR(:,1),FN); % Vector for the extraction 1 DOF
Mfixn = COOR(~Vefn, :); % Matrix of fixed nodes
Vefn = [Vefn; Vefn; Vefn]; % Vector for the extraction 3 DOF

A = Ab(Vefn, :); % Equilibrium matrix for the free nodes

%% FIRST (SVD) SINGULAR VALUE DECOMPOSITION

m = 3 * size(FN, 1); % Total number of infinitesimal mechanism
r = rank(A, tol_1); % Rank of the equilibrium matrix
[U, S, V] = svd (A); % SVD of the equilibrium matrix

s = b - rank(A); % Number of independent self-stress modes
k = m - rank(A); % Number of independent zero-energy deformation modes or
                 % rigid body mechanisms

if s == 0 && k == 0
    disp('Statically determinate and kinematically determinate')
    disp('Check the diag(S) and diag(S2) vectors to change the')
    disp('tolerances of the rank in matrices A and T2 in %% DATA')
elseif s ==0 && k > 0
    disp('Statically determinate and kinematically indetermined')
    disp('Check the diag(S) and diag(S2) vectors to change the')
    disp('tolerances of the rank in matrices A and T2 in %% DATA')
elseif s > 0 && k == 0
    disp('Statically indeterminate and kinematically determined')
elseif s > 0 && k > 0
    disp('Statically indeterminate and kinematically indetermined')
end 

% Plot of the singular values for the first SVD
close

hold on

% Importance of each singular value
subplot(1,2,1)
semilogy(diag(S), '-ok', 'linewidth', 0.1)

title('Singular values - first SVD', 'fontsize', 16, 'linewidth', 0.7)

xlabel('Singular Value (SV)')
ylabel('log(SV)')

% Cumulative energy for each singular value
subplot(1,2,2)
plot(cumsum(diag(S))/sum(diag(S)),'-ok', 'linewidth', 1)

title('Cumulative Energy - first SVD', 'fontsize', 16, 'linewidth', 0.7)

xlabel('Singular Value (SV)')
ylabel('Cumulative Energy for the SV')

set(gcf,'WindowState', 'maximized')

% Independent self-stress mode
T = V(:, ((r+1) : b)); 

%% (DSVD) DOBLE SINGULAR VALUE DECOMPOSITION

% The matrix "e" contains the base vector (ones) of member forces composed
% of a unit stress in the ith group and zero in the other (n - 1) groups
E = zeros(b, n);

for i = 1:n
    boc = CON(CON(:, 4) == sg(i), 5);
    boc = boc(1);
    
    if boc == 0
        bct = -1;
    else
        bct = 1;
    end
    E(CON(:, 4) == sg(i), i) = bct; % Matrix of the base vectors
end

% With the following equations:
% (1) X = T1 * a1 + T2 * a2 + ... + Ts * as
% (2) X = (x1 x1 x1 ... xi xi xi ... xn ... xn)'
% If we do:
% (1)-(2) => 
% X - X = T1 * a1 + ... + Tn * an + (-e1 * x1) + ... + (-en * xn) = 0
% T2 * a2 = 0

% Then T2 can be expressed as follows:
T2 = [T, -E];
[U2, S2, V2] = svd(T2); % Second (SVD) singular value decomposition

r2  = rank(T2, tol_2); % Rank of the matrix T2
s2 = s + n - r2; % Number of integral prestress modes
b2 = size(T2, 2); % Number of columns of T2

% Plot of the singular values for the second SVD
close

hold on

% Importance of each singular value
subplot(1,2,1)
semilogy(diag(S2), '-ok', 'linewidth', 0.1)

title('Singular values - second SVD', 'fontsize', 16, 'linewidth', 0.7)

xlabel('Singular Value (SV)')
ylabel('log(SV)')

% Cumulative energy for each singular value
subplot(1,2,2)
plot(cumsum(diag(S))/sum(diag(S)),'-ok', 'linewidth', 1)

title('Cumulative Energy - second SVD', 'fontsize', 16, 'linewidth', 0.7)

xlabel('Singular Value (SV)')
ylabel('Cumulative Energy for the SV')

set(gcf,'WindowState', 'maximized')

% Prestress of members
xs = V2((size(T, 2) + 1):b2, (r2 + 1):b2); 

% Integral prestress mode
X = E * xs;

%% GRAPH OF THE STRUCTURE AND THE PRESTRESS OF ELEMENTS

close

%   Headers, labels, and legend
CON(:, 6) = X/max(abs(X));

title({SN, 'Integral Prestress Mode'},'fontsize', 16, 'linewidth', 0.7)

text(xmax, ymax + dy, (zmax / 2), ...
    'Traction', 'color', 'w', 'backgroundcolor', [190 190 190]/255, ...
    'linewidth', 1);
text(xmax, ymax, ((zmax / 2) + dz), ...
    'Compression', 'color', 'w', 'backgroundcolor', 'k', 'linewidth', 1);

xlabel('X')
ylabel('Y')
zlabel('Z')

% Plot of the members
t = abs(CON(:, 6)*5); % Line thicknesses
for i = 1:b
    % Initial and final node coordinates for each element
    x = [COOR(CON(i, 2), 2), COOR(CON(i, 3), 2)];
    y = [COOR(CON(i, 2), 3), COOR(CON(i, 3), 3)];
    z = [COOR(CON(i, 2), 4), COOR(CON(i, 3), 4)];
    
    hold on;
          
    if CON(i, 6) < 0.00001 && CON(i, 6) > -0.00001 % Non prestressed member
        plot3(x, y, z, 'r', 'linewidth', t(i));
    elseif CON(i, 6) > 0 % Tractioned member
        plot3(x, y, z, 'color', [190 190 190]/255, 'linewidth', t(i));
    elseif CON(i, 6) < 0 % Compressioned member
        plot3(x, y, z, 'k', 'linewidth', t(i));
    end
end

set(gcf,'WindowState', 'maximized')

% Fixed nodes
if addfixednodes == true
    plot3(Mfixn(:, 2), Mfixn(:, 3), Mfixn(:, 4), 'X', 'MarkerSize', dx,...
      'MarkerEdgeColor','k', 'linewidth', 1);
end

% Number of the nodes
if addnodenumber == true
    text(COOR(:, 2), COOR(:, 3), COOR(:, 4), string(COOR(:, 1)));
end

axis equal;

if sum(abs(COOR(:, 4))) == 0
    axis([xmin xmax ymin ymax]) % Axis limits configuration
else
    axis([xmin xmax ymin ymax zmin zmax]) % Axis limits configuration
    view([xmax ymin zmax]) % Point of view configuration
end

%% RESULTS FROM THE DSVD

% Integral prestress mode with symmetry group
Symmetry_group = sg;

% Unique values of prestress for each symmetry groups
Feasible_integral_prestress_mode = zeros(n, size(X,2));
for i = 1:size(X, 2)
    Feasible_integral_prestress_mode(:, i) = unique(X(:, i), 'stable');
    
    % Prestress normalized vector
    Feasible_integral_prestress_mode(:, i) = ...
    round(Feasible_integral_prestress_mode(:, i)/max(abs(...
    Feasible_integral_prestress_mode(:, i))),3);
end

Result = table(Symmetry_group, Feasible_integral_prestress_mode);
Result = sortrows(Result,1);

disp(Result)