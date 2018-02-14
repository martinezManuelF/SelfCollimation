%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILENAME:         BandDiagrams.m
% COURSE:           EE5322--21st Century Electromagnetics
% INSTRUCTOR:       Raymond C. Rumpf
% NAME:             Manuel F. Martinez
% SEMESTER:         Spring 2018
% DUE DATE:         02/13/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTIALIZE MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% RESTORE STATE
clear all;
close all;
clc;

% UNITS
meters      = 1;
millimeters = 1e-3 * meters;
centimeters = 1e-2 * meters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;
degrees     = pi/180;
F           = 1;
H           = 1;

% CONSTANTS
e0 = 8.85418782e-12 * F/meters;
u0 = 1.25663706e-6 * H/meters;
N0 = sqrt(u0/e0);
c0 = 299792458 * meters/seconds;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DASHBOARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SUPERCELL PARAMETERS
a   = 1;
w   = 0.8 * a;
er1 = 9.0;
er2 = 1.0;

% GRID
Nx = 1024;
Ny = 1024;

% PWEM PARAMETERS
PQ.P    = 11;
PQ.Q    = 11;
MODE.EM = 'H';

% NUMBER OF POINTS IN BAND DIAGRAM
SYM.NP = 40;

% MISC. SIMULATION PARAMETERS
MODE.BC = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BUILD UNIT CELL ONTO HIGH RESOLUTION REAL-SPACE GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE GRID PARAMETERS
dx = a/Nx;
dy = a/Ny;
xa = [0:Nx-1]*dx; xa = xa - mean(xa);
ya = [0:Ny-1]*dy; ya = ya - mean(ya);

% BUILD UNIT CELL
DEV.UR  = ones(Nx,Ny);
DEV.ER  = er1 * ones(Nx,Ny);
h   = sqrt(w^2 - (w/2)^2);    % Height of triangle
ny  = round(h/dy);
ny1 = 1 + floor((Ny - ny)/2);
ny2 = ny1 + ny - 1;
for ny = ny1 : ny2
    f = (ny - ny1 + 1)/(ny2 - ny1 + 1);
    nx = round(f*w/dx);
    nx1 = 1 + floor((Nx - nx)/2);
    nx2 = nx1 + nx - 1;
    DEV.ER(nx1:nx2,ny) = er2;
end 

% VISUALIZE SUPER CELL
figure('Color','w');
c = imagesc(xa,ya,DEV.ER');
axis equal tight;
title('$\textrm{Unit Cell}$','Interpreter','LaTex','FontSize',15);
colormap(gray);
colorbar;
c = get(c,'Parent');
set(c,'FontSize',12);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN PWEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% COMPUTE BLOCH WAVE VECTORS

% Reciprocal Lattice Vectors
DEV.T1 = (2*pi/a) * [1 ; 0];
DEV.T2 = (2*pi/a) * [0 ; 1];

% Key points of symmetry
NP1 = SYM.NP;
D = [0 ; 0];
X = 0.5*DEV.T2;
Z = 0.5*DEV.T1 + 0.5*DEV.T2;
M = 0.5*DEV.T1;
E = 0.5*DEV.T1 - 0.5.*DEV.T2;
G = -0.5*DEV.T2;

if (MODE.BC == 1)
    % Generate list of Bloch Wave Vectors
    bx = [ linspace(D(1),X(1),NP1), linspace(X(1),Z(1),NP1) , ...
            linspace(Z(1),M(1),NP1) , linspace(M(1),E(1),NP1) , ...
            linspace(E(1),G(1),NP1) , linspace(G(1),D(1),NP1) ];
    by = [ linspace(D(2),X(2),NP1), linspace(X(2),Z(2),NP1) , ...
            linspace(Z(2),M(2),NP1) , linspace(M(2),E(2),NP1) , ...
            linspace(E(2),G(2),NP1) , linspace(G(2),D(2),NP1)];
    BETA = [ bx ; by ];
elseif (MODE.BC == 2)
    bx      = linspace(-pi,+pi,SYM.NP);
    by      = linspace(+pi,-pi,SYM.NP);
    BETA    = zeros(2,SYM.NP^2);
    C       = 1;
    for x = 1 : SYM.NP
        for y = 1 : SYM.NP
            BETA(:,C) = [bx(x) by(y)];
            C = C + 1;
        end
    end
end

% GENERATE MISSING STRUCTURE PARAMETERS
DEV.LATTICE = a;

% PERFORM PWEM SIMULATION
WN = pwem2d(DEV,SYM,BETA,PQ,MODE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% POST-PROCESS DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (MODE.BC == 1)
    figure('Color','w');
    h = plot([1:length(bx)], WN, '-b','LineWidth',1.5);
    hold on;
    ylim([0 0.6]);
    xlim([1 length(bx)]);
    h2 = get(h(1),'Parent');
    set(h2,'FontSize',13,'LineWidth',1.5);
    title([MODE.EM '$\textrm{ Mode}$'],'FontSize',16,'Interpreter','LaTex')
    
    % Format x-Axis
    XT = [1 NP1 2*NP1 3*NP1 4*NP1 5*NP1 6*NP1];
    XL = {'\Delta','X','Z','M','\Sigma','\Gamma','\Delta'};
    set(h2,'XTick',XT,'XTickLabel',XL);
    xlabel('Bloch Wave Vector ($\beta$)','Interpreter','LaTex','FontSize',14);
    
    % Format y-Axis
    ylabel('Normalized Frequency $\frac{\omega a}{2\pi c_0}$','Interpreter',...
        'LaTex','FontSize',14,'Rotation',90);
    
%     % Add Vertical Lines
    for n = 2 : length(XT)-1
        x = [ XT(n) XT(n) ];
        y = [0 1];
        line(x,y,'Color','k','LineWidth',0.75,'LineStyle','--');
    end
    
    % Draw BZ and IBZ
    SQx = [50 60 70];
    SQy = [0.05 0.1 0.15];
    vx  = [SQx(2) SQx(2) SQx(3) SQx(3) SQx(2)];
    vy  = [SQy(1) SQy(3) SQy(3) SQy(1) SQy(1)];
    fill(vx,vy,[198/255 226/255 255/255],'EdgeColor',...
        [198/255 226/255 255/255]);
    text(SQx(2)-0.5,SQy(2),'$\Delta$','Interpreter','LaTex',...
        'HorizontalAlignment','Right');
    text(SQx(2)-0.5,SQy(3)+0.005,'X','Interpreter','LaTex',...
        'HorizontalAlignment','Right');
    text(vx(3)+0.5,vy(3),'Z','Interpreter','LaTex',...
        'HorizontalAlignment','Left');
    text(SQx(3)+0.5,SQy(2),'M','Interpreter','LaTex',...
        'HorizontalAlignment','Left');
    text(SQx(3)+0.5,SQy(1),'$\Sigma$','Interpreter','LaTex',...
        'HorizontalAlignment','Left');
    text(SQx(2)-0.5,SQy(1)+0.005,'$\Gamma$','Interpreter','LaTex',...
        'HorizontalAlignment','Right');

    for x = 1 : length(SQx)
        for y = 1 : length(SQy)
            plot(SQx(x),SQy(y),'ok','MarkerSize',3,'MarkerFaceColor','k');
        end
    end
    a = line([SQx(1) SQx(3)],[SQy(1) SQy(1)]);
    set(a,'LineWidth',1,'Color','k');
    a = line([SQx(1) SQx(1)],[SQy(1) SQy(3)]);
    set(a,'LineWidth',1,'Color','k');
    a = line([SQx(3) SQx(3)],[SQy(1) SQy(3)]);
    set(a,'LineWidth',1,'Color','k');
    a = line([SQx(1) SQx(3)],[SQy(3) SQy(3)]);
    set(a,'LineWidth',1,'Color','k');
    hold off;
elseif (MODE.BC == 2)
    figure('Color','w');
    % PLOT FULL BAND DIAGRAM
    for b = 1 : 7               % Just plot some bands
        a = surf(reshape(WN(b,:),SYM.NP,SYM.NP));
        shading interp;
        hold on;
        drawnow;
        view([-40,10]);
    end
    hold off;
    a = get(a,'Parent');
    set(a,'FontSize',12);
    title([MODE.EM '$\textrm{ Mode Full Band Diagram}$'],...
        'Interpreter','LaTex','FontSize',18);
    xlabel('\beta_x','FontSize',15);
    ylabel('\beta_y','FontSize',15);
    zlabel('Normalized Frequency $\frac{\omega a}{2\pi c_0}$',...
        'Interpreter','LaTex','FontSize',15);
    box on;
    colormap(jet);
    T = [1 20 40];
    L = {'-\pi','0','\pi'};
    set(gca,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
    
    % ISOFREQUENCY CONTOURS
    figure('Color','w');
    [C,H] = contour(reshape(WN(1,:),SYM.NP,SYM.NP),8,'LineWidth',1.5);
    set(gca,'FontSize',13);
    title([MODE.EM '$\textrm{ Mode 1st Band Isofrequency Contour}$'],...
        'FontSize',15,'Interpreter','LaTex');
    xlabel('\beta_x','FontSize',12);
    ylabel('\beta_y','FontSize',12);
    clabel(C,H,'FontSize',10);
    colormap(jet);
    T = [1 20 40];
    L = {'-\pi','0','\pi'};
    set(gca,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
    colorbar;
    axis equal tight;
    
    figure('Color','w');
    [C,H] = contour(reshape(WN(2,:),SYM.NP,SYM.NP),8,'LineWidth',1.5);
    set(gca,'FontSize',13);
    title([MODE.EM '$\textrm{ Mode 2nd Band Isofrequency Contour}$'],...
        'FontSize',15,'Interpreter','LaTex');
    xlabel('\beta_x','FontSize',12);
    ylabel('\beta_y','FontSize',12);
    clabel(C,H,'FontSize',10);
    colormap(jet);
    T = [1 20 40];
    L = {'-\pi','0','\pi'};
    set(gca,'XTick',T,'XTickLabel',L,'YTick',T,'YTickLabel',L);
    colorbar;
    axis equal tight;
    
end
    