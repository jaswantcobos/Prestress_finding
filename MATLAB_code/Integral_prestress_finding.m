%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code to determine the Integral Feasible Prestress of Tensegrity Structures
% using the Doble Singular Value Decomposition (DSVD) method developed by 
% X.F.Yuan and S.L.Dong (https://doi.org/10.1016/S0045-7949(03)00254-2)

%	Undergraduate Thesis:
%  Tensegrity as an Alternative Structural System Applied to Pedestrian Bridges
%  University of Cuenca

%	By: 
%  Jaswant Cobos

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Data:

% This program read the information from an Excel document (data.xlsx)
% The user should order the sheets from the document and the data like this:

% Sheet 1: Connectivity matrix (CON) (don't insert headers)
%  Insert data in this order: |Element|Initial_Node|Final_Node|Symmetry_Group| 
%   Symmetry groups must be ordered from least to greatest

% Sheet 2: Coordinate matrix (COOR) (don't insert headers)
%  Insert data in this order: |Node|x_Coordinates|y_Coordinates|y_Coordinates|

% Sheet 3: Free nodes vector (NL) (don't insert headers)
%  Insert data in the first column of this sheet

CON=xlsread('data.xlsx',1); % Connectivity matrix
COOR=xlsread('data.xlsx',2); % Coordinates matrix
NL=xlsread('data.xlsx',3); % Free Nodes Vector

b=size(CON,1); % Number of elements
n=size(NL,1); % Number of free nodes
Cs=zeros(b,n); % Conectivity Matrix

for i=1:size(Cs,1)
    Cs(i,min(CON(i,2:3)))=1;
    Cs(i,max(CON(i,2:3)))=-1;
end

%% Equilibrium matrix assembly:

% Lenght of elements:

for i=1:size(CON,1)
    NI=CON(i,2); % Initial node of the element
    NF=CON(i,3); % Final node of the element
    fNI = find( COOR(:,1)==NI); % Row position of the initial node
    fNF = find( COOR(:,1)==NF); % Row position of the final node
    CI=COOR(fNI,2:4); % X, Y, and Z initial node coordinates
    CF=COOR(fNF,2:4); % X, Y, and Z final node coordinates
    % Element Lenght:
    L(i,1)=((CF(1)-CI(1))^2+(CF(2)-CI(2))^2+(CF(3)-CI(3))^2)^(1/2);
end
L=diag(L);

% Equilibrium matrix:

A=[(Cs'*diag(Cs*COOR(:,2))*inv(L))',(Cs'*diag(Cs*COOR(:,3))*inv(L))',(Cs'*diag(Cs*COOR(:,4))*inv(L))']';

%% Determination of independent prestressing modes:

r1=rank(A); % A matrix rank.
s1=b-r1;
[~,~,V]=svd(A); % Singular value decomposition of the equibrium matrix
T=V(:,r1+1:b); % Independent prestressing modes:

%% Determination of integral prestressing modes:

[a,e]=size(T); % T matrix size
n=max(CON(:,4)); % Number of symmetry groups 
E=zeros(a,n); % This matrix will have the -e vectors from the Yuan method 
i=1; 
j=1; 

% Assembly of the -e vectors from the Yuan method:

for k=1:n;      
    while (CON(j,4)-k)==0; 
     E(j,i)=-1; 
     j=j+1; 
     if j==a+1; 
         j=1; 
     end 
    end 
    i=i+1; 
end

X=[T E]; % X matrix from the Yuan method 
[~,S,V2] = svd(X); % Second singular value decomposition 
r2=rank(X); % X matrix rank.
b2=size(X,2); % Number of column vectors of X matrix
s2=s1+n-r2; % Number of integral prestressing modes
Tr=V2(:,r2+1:b2); 
W1=E*Tr(size(T,2)+1:b2,:); % Integral prestressing mode

% The next "for" orders the integral presstressing modes according to the 
% symmetry group

for i=1:size(W1,1);
    if i==1
        W(i,1)=W1(i);
    else
        if W1(i)==W1(i-1) %<1.01&W1(i)/W1(i-1)>0.99
            W(i,1)=0;
        else
            W(i,1)=W1(i);
        end
    end
end
h=find(W==0); % Empty spaces in W vector
W(h)=[]; 

% The next "if" calculates the prestressing relative distribution according to
% the greater value

if abs(min(W))>max(W)
    W=W/abs(min(W));
else
    W=W/max(W);
end


%% Plot:

b=size(CON,1); 

% Configuration of the plot boundaries:

dx=(max(COOR(:,2))-min(COOR(:,2)))/10; 
dy=(max(COOR(:,3))-min(COOR(:,3)))/10; 
dz=(max(COOR(:,4))-min(COOR(:,4)))/10; 
xmin=min(COOR(:,2))-dx; 
xmax=max(COOR(:,2))+dx; 
ymin=min(COOR(:,3))-dy; 
ymax=max(COOR(:,3))+dy; 
zmin=min(COOR(:,4))-dz; 
zmax=max(COOR(:,4))+dz; 

% Structure Plot

title('Prestressing Mode','fontsize',20,'linewidth',0.7) 
text((xmax-0.5*dx),(ymax-3*dy),(zmax-3*dz),'Traction','color','b','backgroundcolor','w','linewidth',0.7);
text((xmax-0.5*dx),(ymax-5*dy),(zmax-5*dz),'Compression','color','r','backgroundcolor','w','linewidth',0.7);
text((xmax-0.5*dx),(ymax-7*dy),(zmax-7*dz),'No stress','color','k','backgroundcolor','w','linewidth',0.7);

for i=1:b; 
    x=[COOR(CON(i,2),2),COOR(CON(i,3),2)]; 
    y=[COOR(CON(i,2),3),COOR(CON(i,3),3)]; 
    z=[COOR(CON(i,2),4),COOR(CON(i,3),4)]; 
    hold on; 
    grid on; 
    axis([xmin xmax ymin ymax zmin zmax]); 
    if W(CON(i,4))<0.00001&&W(CON(i,4))>-0.00001; 
        plot3(x,y,z,'k','linewidth',2); 
    elseif W(CON(i,4))>0; 
        g=round(abs(W(CON(i,4))*10)); % Line thickness proportional to 
        if g==0;                      % the distribution value
            g=0.5; 
        end; 
        plot3(x,y,z,'b','linewidth',g); 
    elseif W(CON(i,4))<0 ; 
        g=round(abs(W(CON(i,4))*10));  % Line thickness proportional to 
        if g==0;                       % the distribution value
            g=0.5; 
        end; 
        plot3(x,y,z,'r','linewidth',g); 
    end 
end 
axis equal;