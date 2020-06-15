% close all;

% width = 8;
% height = 6;
% VISIBLE=1;
% The properties we've been using in the figures
% load('./plotting/colormapATDD.mat');
% load('./plotting/colormapMASK.mat');
% 
% load('./plotting/colormap_Phase.mat');
% load('.\plotting\colormapDEM.mat');
% colormapATDD = colormap('hot');
% colormapATDD(1,:)=[1 1 1];
set(0,'defaultLineLineWidth',3);   % set the default line width
%set(0,'defaultLineMarkerSize',4);  % set the default line marker size
set(0,'defaultAxesFontName','Arial');
set(0,'defaultAxesFontSize',23); %30
set(0,'defaultTextFontName','Arial'); 
set(0,'defaultTextFontSize',23); %30

% % Make ALL the plots invisible
% if VISIBLE == 1
%     set(0,'defaultFigureVisible','on');
% else
%     set(0,'defaultFigureVisible','off');
% end

% Set the default Size for display
mida = get(0,'ScreenSize');
mida(3:4)=[1920,1080];
% set(0,'defaultFigurePosition', [defpos(1) defpos(2)-50 width*100, height*100]);
set(0,'defaultFigurePosition',mida);
colormapATDD=colormap(hot(128));
colormapATDD(1,:)=[1 1 1];
set(0,'defaultFigureColormap',colormapATDD);
% close(0);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
% set(0,'defaultFigurePaperUnits','centimeters');
% defsize = get(gcf, 'PaperSize');
% left = (defsize(1)- width)/2;
% bottom = (defsize(2)- height)/2;
% defsize = [left, bottom, width, height];
% set(0, 'defaultFigurePaperPosition', defsize);
set(0,'defaultFigurePaperUnits','points');
% set(0,'defaultFigurePaperPosition', mida);

% plot(...)
% print('figName','-r100','-dpng')