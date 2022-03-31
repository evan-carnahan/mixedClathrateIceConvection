% code for plotting contours of convective shutdown
close all; clear all;
set(groot,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultTextFontSize',8)
fs = 8;
yLab = 'mixture viscosity, $J$';
%% load in relavent matricies
load('titanParamSweepContourPlot_all.mat','dGr','etaGr','JGr','tShut','cLid','botFlux')
tST = tShut;
cLT = cLid;
bFT = botFlux;

load('plutoParamSweepContourPlot_all.mat','tShut','cLid','botFlux')
tSP = tShut;
cLP = cLid;
bFP = botFlux;

% etaUni = etaGr(1,:);
% lEtaUni = log10(etaUni);
gVisc_func = @(d) 1/2*(d*1e-3).^2*264/(9.1e-8)*exp(59e3/(264*8.314)); % d in mm
dLab = 0.5:0.5:3.5;
etaLab = gVisc_func(dLab);
lEtaLab = round(log10(etaLab),1);

% just markers at 1e14, 1e15, and 1e16
dUniLab = [0.38,1.2085,3.83];
etaWhole = [14,15,16];

%% set up plot
fig = figure;%('units','normalized','Position',[0.1,0.1,0.9,0.9]);

%% Titan
%% Convective shutdown
subplot(2,3,1)
ax1 = gca;
ax1.FontSize = fs;
hold on
contourf(dGr,JGr,tST,30,'linestyle','none')
c1 = colorbar;
c1.Location = 'EastOutside';
c1.Label.String = 'Time convection stops, Gyrs';
caxis([0 4])
% xlabel('Grain size, mm','FontSize',fs)
ylabel(yLab,'FontSize',fs);
ax1.XTickLabel = [];
box on

%% Conductive lid thickness
subplot(2,3,2)
ax2 = gca;
ax2.FontSize = fs;
hold on
contourf(dGr,JGr,cLT,30,'linestyle','none')
c2 = colorbar;
c2.Location = 'EastOutside';
c2.Label.String = 'Conductive lid thickness, km';
caxis([0 100])
[cn2,h2] = contour(dGr,JGr,cLT,'w','LevelList',40,'LineWidth',1);
% xlabel('Grain size, mm','FontSize',fs);
ylabel(yLab,'FontSize',fs);
ax2.XTickLabel = [];
box on

%% Basal Heat Flux
subplot(2,3,3)
ax3 = gca;
ax3.FontSize = fs;
hold on
contourf(dGr,JGr,bFT,30,'linestyle','none')
[cn3,h3] = contour(dGr,JGr,bFT,'w','LevelList',[3.2 7.1],'LineWidth',1);
c3 = colorbar;
c3.Location = 'EastOutside';
c3.Label.String = 'Geothermal heat flow, mW/m$^2$';
% xlabel('Grain size, mm','FontSize',fs)
ylabel(yLab,'FontSize',fs);
ax3.XTickLabel = [];
box on

%% Pluto
%% Time till convective shutdown
subplot(2,3,4)
ax4 = gca;
ax4.FontSize = fs;
hold on
[con,hCon] = contourf(dGr,JGr,tSP,30,'linestyle','none');
c4 = colorbar;
c4.Location = 'EastOutside';
c4.Label.String = 'Time convection stops, Gyrs';
caxis([0 4])
xlabel('Ice grain size, mm','FontSize',fs)
ylabel(yLab,'FontSize',fs);
box on

%% Conductive lid thickness
subplot(2,3,5)
ax5 = gca;
ax5.FontSize = fs;
hold on
contourf(dGr,JGr,cLP,30,'linestyle','none')
c5 = colorbar;
c5.Location = 'EastOutside';
c5.Label.String = 'Conductive lid thickness, km';
caxis([0 200])
% [cn2,h2] = contour(dGr,JGr,cLP,'r','LevelList',200,'LineWidth',1);
xlabel('Ice grain size, mm','FontSize',fs);
ylabel(yLab,'FontSize',fs);
box on

%% Basal heat flux
subplot(2,3,6)
ax6 = gca;
ax6.FontSize = fs;
hold on
contourf(dGr,JGr,bFP,30,'linestyle','none')
[cn6,h6] = contour(dGr,JGr,bFP,'w','LevelList',3.35,'LineWidth',1);
c6 = colorbar;
c6.Location = 'EastOutside';
c6.Label.String = 'Geothermal heat flow, mW/m$^2$';
xlabel('Ice grain size, mm','FontSize',fs)
ylabel(yLab,'FontSize',fs);
box on

%% area that satisfies both elastic and thermal constraints

% % titan
% subplot(2,3,1)
% hT = patch([cn2(1,2:end),3.83,cn3(1,2:13)],[cn2(2,2:end),-1,cn3(2,2:13)],linspace(0,1,size([cn2(2,2:end),-1,cn3(2,2:13)],2)),'EdgeColor','none');
% hatchfill2(hT,'single','HatchAngle',45,'HatchDensity',30,'HatchColor','r','HatchLineWidth',0.1);
% hT.FaceAlpha = 0;
% plot(cn2(1,2:end),cn2(2,2:end),'r','LineWidth',1)
% plot(cn3(1,2:13),cn3(2,2:13),'r','LineWidth',1)
% 
% % pluto
% subplot(2,3,4)
% hP = patch([cn6(1,2:10),3.83,3.83],[cn6(2,2:10),1,-1],'r','FaceAlpha',0.25,'EdgeAlpha',0);
% hatchfill2(hP,'single','HatchAngle',45,'HatchDensity',30,'HatchColor','r','HatchLineWidth',0.1);
% hP.FaceAlpha = 0;
% plot(cn6(1,2:10),cn6(2,2:10),'r','LineWidth',1)

%% figure orientation
f = gcf;
f.Units = 'centimeters';
width = 19;
height = 10.7;
f.Position(3:4) = [width height];
drawnow
%% [left bottom width height]
ax1.Units = 'centimeters';
ax2.Units = 'centimeters';
ax3.Units = 'centimeters';
ax4.Units = 'centimeters';
ax5.Units = 'centimeters';
ax6.Units = 'centimeters';

c1.TickLabelInterpreter = 'latex';
c1.Label.Interpreter = 'latex';
c1.Units = 'centimeters';
c2.TickLabelInterpreter = 'latex';
c2.Label.Interpreter = 'latex';
c2.Units = 'centimeters';
c3.TickLabelInterpreter = 'latex';
c3.Label.Interpreter = 'latex';
c3.Units = 'centimeters';
c4.TickLabelInterpreter = 'latex';
c4.Label.Interpreter = 'latex';
c4.Units = 'centimeters';
c5.TickLabelInterpreter = 'latex';
c5.Label.Interpreter = 'latex';
c5.Units = 'centimeters';
c6.TickLabelInterpreter = 'latex';
c6.Label.Interpreter = 'latex';
c6.Units = 'centimeters';


c1.TickLabels(end) = {'4+'};
c4.TickLabels(end) = {'4+'};

%% add text descriptions to contours
subplot(2,3,2)
hold on
text(2,-0.08, 'min elastic thickness','HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontSize',7,'Rotation',-29,'Color','w','FontWeight','bold')

subplot(2,3,3)
hold on
text(2,0.26, 'min heat flow','HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontSize',7,'Rotation',-25,'Color','w','FontWeight','bold')
text(0.6,-0.33, 'max heat flow','HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontSize',7,'Rotation',-65,'Color','w','FontWeight','bold')

% subplot(2,3,5)
% hold on
% text(2.5,-0.4, 'conductive at present','HorizontalAlignment','center',...
%     'VerticalAlignment','Bottom','FontSize',7,'Rotation',-90,'Color','r')

subplot(2,3,6)
hold on
text(2.5,-0.4, 'nominal heat flow','HorizontalAlignment','center',...
    'VerticalAlignment','Bottom','FontSize',7,'Rotation',-90,'Color','w','FontWeight','bold')

%% axis positioning
bMarg = 0.7;
mMarg = 0.5;
lMarg = 1.1;
cMarg1 = 2;
cMarg2 = 2.2;
cSz = 4.2;
cWid = 0.2;

ax4.Position = [lMarg bMarg cSz cSz];
ax1.Position = [lMarg bMarg+cSz+mMarg cSz cSz];
ax5.Position = [lMarg+cSz+cMarg1 bMarg cSz cSz];
ax6.Position = [lMarg+2*cSz+cMarg1+cMarg2 bMarg cSz cSz];
ax2.Position = [lMarg+cSz+cMarg1 bMarg+cSz+mMarg cSz cSz];
ax3.Position = [lMarg+2*cSz+cMarg1+cMarg2 bMarg+cSz+mMarg cSz cSz];

c1.Position(1) = ax1.Position(1) + ax1.Position(3) + 0.1;
c2.Position(1) = ax2.Position(1) + ax2.Position(3) + 0.1;
c3.Position(1) = ax3.Position(1) + ax3.Position(3) + 0.1;
c4.Position(1) = ax4.Position(1) + ax4.Position(3) + 0.1;
c5.Position(1) = ax5.Position(1) + ax5.Position(3) + 0.1;
c6.Position(1) = ax6.Position(1) + ax6.Position(3) + 0.1;

c1.Position(3) = cWid;
c2.Position(3) = cWid;
c3.Position(3) = cWid;
c4.Position(3) = cWid;
c5.Position(3) = cWid;
c6.Position(3) = cWid;

c1.Label.Position(1) = 2.6;
c4.Label.Position(1) = 2.6;
c2.Label.Position(1) = 4;
c5.Label.Position(1) = 4;


% set super x-axis
axX1 = axes('Units','centimeters','position',get(ax1,'position'),'color','none','XAxisLocation','top');
axX1.YTick = [];
axX1.XLim = [0.38,3.83];
axX1.XTick = dUniLab;
axX1.XTickLabel = etaWhole;
axX1.XLabel = ax1.XLabel;
axX1.XLabel.String = 'log$_{10}$ viscosity, Pa s';
axX1.Position(3:4) = ax1.Position(3:4);
axX1.FontSize = fs;
axX1.TickLength = [0 0];
% ax1.XLabel.String = 'Grain size, mm';
% ax1.XLabel.FontSize = 8;


axX2 = axes('Units','centimeters','position',get(ax2,'position'),'color','none','XAxisLocation','top');
axX2.YTick = [];
axX2.XLim = [0.38,3.83];
axX2.XTick = dUniLab;
axX2.XTickLabel = etaWhole;
axX2.XLabel = ax2.XLabel;
axX2.XLabel.String = 'log$_{10}$ viscosity, Pa s';
axX2.Position(3:4) = ax2.Position(3:4);
axX2.FontSize = fs;
axX2.TickLength = [0 0];
% ax2.XLabel.String = 'Grain size, mm';
% ax2.XLabel.FontSize = 8;

axX3 = axes('Units','centimeters','position',get(ax3,'position'),'color','none','XAxisLocation','top');
axX3.YTick = [];
axX3.XLim = [0.38,3.83];
axX3.XTick = dUniLab;
axX3.XTickLabel = etaWhole;
axX3.XLabel = ax3.XLabel;
axX3.XLabel.String = 'log$_{10}$ viscosity, Pa s';
axX3.Position(3:4) = ax3.Position(3:4);
axX3.FontSize = fs;
axX3.TickLength = [0 0];
% ax3.XLabel.String = 'Grain size, mm';
% ax3.XLabel.FontSize = 8;


axX4 = axes('Units','centimeters','position',get(ax4,'position'),'color','none','XAxisLocation','top');
axX4.XTick = [];
axX4.YTick = [];
% axX4.XLim = [0.38,3.83];
% axX4.XTick = dUniLab;
% axX4.XTickLabel = etaWhole;
% axX4.XLabel = ax4.XLabel;
% axX4.XLabel.String = 'log$_{10}$ viscosity, Pa s';
% axX4.Position(3:4) = ax4.Position(3:4);
% axX4.FontSize = fs;
% axX4.TickLength = [0 0];
% ax4.XLabel.String = 'Grain size, mm';
% ax4.XLabel.FontSize = 8;
% 
axX5 = axes('Units','centimeters','position',get(ax5,'position'),'color','none','XAxisLocation','top');
axX5.XTick = [];
axX5.YTick = [];
% axX5.XLim = [0.38,3.83];
% axX5.XTick = dUniLab;
% axX5.XTickLabel = etaWhole;
% axX5.XLabel = ax5.XLabel;
% axX5.XLabel.String = 'log$_{10}$ viscosity, Pa s';
% axX5.Position(3:4) = ax5.Position(3:4);
% axX5.FontSize = fs;
% axX5.TickLength = [0 0];
% ax5.XLabel.String = 'Grain size, mm';
% ax5.XLabel.FontSize = 8;
% 
axX6 = axes('Units','centimeters','position',get(ax6,'position'),'color','none','XAxisLocation','top');
axX6.XTick = [];
axX6.YTick = [];
% axX6.XLim = [0.38,3.83];
% axX6.XTick = dUniLab;
% axX6.XTickLabel = etaWhole;
% axX6.XLabel = ax6.XLabel;
% axX6.XLabel.String = 'log$_{10}$ viscosity, Pa s';
% axX6.Position(3:4) = ax6.Position(3:4);
% axX6.FontSize = fs;
% axX6.TickLength = [0 0];
% ax6.XLabel.String = 'Grain size, mm';
% ax6.XLabel.FontSize = 8;

%% annotate isostrain/isostress labels
% fsIso = ax1.YLabel.FontSize;
ax1.YTickLabel(1) = {'isostress'};
ax1.YTickLabel(end) = {'isostrain'};
ax1.YLabel.Position(1) = -0.1;


ax4.YTickLabel(1) = {'isostress'};
ax4.YTickLabel(end) = {'isostrain'};
ax4.YLabel.Position(1) = -0.1;



%% annotate figure labels
fsM = 1;
lShif = 0.8;
uShif = 0.09;
annotation('textbox','Units','centimeters','Position',[ax1.Position(1)-lShif ax1.Position(2)+ax1.Position(4)+uShif 0.1 0.1],...
    'string','a','fontsize',fsM*fs,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox','Units','centimeters','Position',[ax2.Position(1)-lShif ax2.Position(2)+ax2.Position(4)+uShif 0.1 0.1],...
    'string','b','fontsize',fsM*fs,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox','Units','centimeters','Position',[ax3.Position(1)-lShif ax3.Position(2)+ax3.Position(4)+uShif 0.1 0.1],...
    'string','c','fontsize',fsM*fs,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox','Units','centimeters','Position',[ax4.Position(1)-lShif ax4.Position(2)+ax4.Position(4)+uShif 0.1 0.1],...
    'string','d','fontsize',fsM*fs,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox','Units','centimeters','Position',[ax5.Position(1)-lShif ax5.Position(2)+ax5.Position(4)+uShif 0.1 0.1],...
    'string','e','fontsize',fsM*fs,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')
annotation('textbox','Units','centimeters','Position',[ax6.Position(1)-lShif ax6.Position(2)+ax6.Position(4)+uShif 0.1 0.1],...
    'string','f','fontsize',fsM*fs,...
    'fontweight','bold','interpreter','tex',...
    'HorizontalAlignment','left','VerticalAlignment','bottom',...
    'Margin',0,'LineStyle','none')

