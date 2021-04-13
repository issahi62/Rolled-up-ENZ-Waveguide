%%%% Developed by : Ibrahim Issah 
%%%% Supervisor   : Assoc. Prof. Humeyra Caglayan
%%%% Faculty of Engineering and Natural Sciences, Photonics, 
%%%% Tampere University
%%%% https://research.tuni.fi/metaplasmonics/

%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE %%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clc; 
clear;


%%%%%%%%%%%%%%%%%%%%%%% DATA IMPORTATION %%%%%%%%%%%%%%%%%%%%%%%%%
load dipole_dipole_normalized.mat; %% Save the normalized decay rate and the dipole-dipole interactions into a file name of your choice 
%%%%%%%%%%%%%% Here we used the above name


r = lum.x0; %% interatomic distance 
gamma = lum.y1; %% decay rate gamma_12
gparam_norm =lum.y0; %dipole-dipole interactions -- g12 

%% Selection of some section of concurrence to plot as a function of time
indt1 = 200; %% INDEX PLOT FOR TIME CONCURRENCE NOT WITH POSITION
indt2 = 413; 


%%%%%%%%%%%%%%%%%%%%%%% DASHBOARD %%%%%%%%%%%%%%%%%%%%%%%%%
time_value = 20; %% time -- normalized time with gamma11 ==1 
t = linspace(0, time_value,413);  
gg = gamma(ceil(length(gamma)/2)); %% self interactions g11  

%% concurrence section 
B = exp(-(gg+gamma).*t');  %% evoluation of sub and superradiant density matrix
C = exp(-(gg-gamma).*t'); 
D = (B-C).^2; 
E = 4*exp(-(2.*gg.*t)).*sin(2.*gparam_norm'.*t).^2;
Conc = 0.5*sqrt(D + E.'); %% concurrece formulations

gamma_0 = gamma(ceil(length(gamma)/2)); 
tt = linspace (0,2,431);
tau = tt;


%%%%%%%%%%%%%%%% PLOTTING SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% POWER PLOT %%%%%%%%%%%%%%%
% 
rho_ee = exp(-tau');
rho_aa = exp(-(gamma_0-gamma).*tau');
rho_ss =exp(-(gamma_0+gamma).*tau');
 
 P = 2*gamma_0*rho_ee + (gamma_0 + gamma).*rho_ss + (gamma_0 - gamma).*rho_aa;


%%%%%%%%%%%%%%%%%%%%%%% PLOTTING SECTION %%%%%%%%%%%%%%%%%%%%%%%%%
fig2 = figure('color', 'w'); 
ax1 = axes;
plot(lum.x0, lum.y0,lum.x1, lum.y1, 'Linewidth', 4)
hold on 
plot(lum.x0, 0*lum.y0, 'k--', 'Linewidth', 2); 
set(ax1, 'XLim', [0 2.09351])
set(ax1, 'YLim', [-1 1])
set(gca, 'Fontsize', 20); 
xlabel('Normalized inter-emitter positon r_{12}/\lambda_0', 'FontSize',21);
ylabel('Coupling and decay rate', 'Fontsize', 21); 
legend('\Gamma_{12}/ \gamma','U_{12}/ \gamma')



% Create figure
figure1 = figure('color', 'w');
% Create axes
axes1 = axes('Parent',figure1,'FontSize',31);
view(axes1,[90 -90]);
grid(axes1,'on');
hold(axes1,'all');
%%%%%%%%%%%%%%%%%%%%%%% PLOTTING SECTION %%%%%%%%%%%%%%%%%%%%%%%%%
[x, y]=meshgrid(t,r);
surf(x,y,Conc.','Parent',axes1,'LineStyle','none',...
    'FaceColor','interp');
% Create colorbar
colorbar('peer',axes1);
colormap('hot'); 
caxis([0 .5])
ylim([0.01 2])

set(gca,'Fontsize', 31); 
%xlabel('Normalized interatomic distance r_{12}/\lambda_{0}'); 
%ylabel('Normalized time t\gamma'); 
xlabel('Normalized time t\gamma','VerticalAlignment','bottom','Rotation',90,...
    'HorizontalAlignment','center',...
    'FontSize',31);

ylabel('Normalized interatomic distance r_{12}/\lambda_{0}','VerticalAlignment','cap',...
    'HorizontalAlignment','center',...
    'FontSize',31);
%title('Concurrence (C)');


%% pOWER PLOT 
% Create figure
figure1 = figure('color', 'w');

% Create axes
axes1 = axes('Parent',figure1,'FontSize',31);
view(axes1,[90 -90]);
grid(axes1,'on');
hold(axes1,'all');
%%%%%%%%%%%%%%%%%%%%%%% PLOTTING SECTION %%%%%%%%%%%%%%%%%%%%%%%%%
[x, y]=meshgrid(tt,r);
surf(x,y,P.','Parent',axes1,'LineStyle','none',...
    'FaceColor','interp');

% Create colorbar
colorbar('peer',axes1);
colormap('hot'); 
%caxis([0.7 1.2])
 ylim([0.01 2])
%xlim([0, .5]); 

xlabel('Normalized time t\gamma','VerticalAlignment','bottom','Rotation',90,...
    'HorizontalAlignment','center',...
    'FontSize',31);

ylabel('Normalized interatomic distance r_{12}/\lambda_{0}','VerticalAlignment','cap',...
    'HorizontalAlignment','center',...
    'FontSize',31);
%title('Power enhancemet')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig4 = figure('color', 'w'); 
ax1 = axes;
plot(lum.x0, Conc(indt1, :), 'Linewidth', 4)
hold on 
plot(lum.x0, Conc(indt2, :), 'Linewidth', 4)
%plot(lum.x0, 0*lum.y0, 'k--', 'Linewidth', 2); 
xlim([0.01 2])
set(gca, 'Fontsize', 20); 
xlabel('Normalized inter-emitter positon r_{12}/\lambda_0', 'FontSize',21);
ylabel('C(\infty)', 'Fontsize', 21); 
legend(strcat('t = ',num2str(ceil(t(indt1))), 's'), strcat('t = ',num2str(ceil(t(indt2))), 's'))