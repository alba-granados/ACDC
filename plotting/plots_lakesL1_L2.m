global wv_length_ku pi_cst

[lakeTiff, R] = geotiffread(['./inputs/mask/lake2.tif']);
AspectRatio   = R.RasterSize(2)/R.RasterSize(1);
Reduced_Size  = round(R.RasterSize/5);
latMargin    = Reduced_Size/2*R.DeltaLat;
lonMargin    = Reduced_Size/2*R.DeltaLon;

h=figure;
[~,sample_retracked] = max(dataL1_single.Wvf(dataL2_single.index_inside(i_surf),:));
D = 1.167;
X_dist =  (dataL1_single.alt(dataL2_single.index_inside(i_surf))*1e-3- dataL2_single.fineH(dataL2_single.index_inside(i_surf)))* dataL1_single.PhD(dataL2_single.index_inside(i_surf),sample_retracked)/1e6*wv_length_ku/2/pi_cst/D;
% subplot(2,2,1); mapshow(lakeTiff, R);hold on;axis image off; axis image on;hold on
subplot(2,3,1); mapshow(lakeTiff, R);hold on;axis image off; axis image on;hold on
figlabels('Longitude','Latitude','','Woerba Lake',12);
max_error = 10; %[m]
min_error = 0.3; %[m]
dataL2_single.s = zeros(1,length(dataL2_single.fineH))+40;
mean_lake_height = 5195; % Woerba
for iPoint=1:length(dataL2_single.index_inside)
    if(abs(dataL2_single.fineH(dataL2_single.index_inside(iPoint))-mean_lake_height)>max_error)
        dataL2_single.c(iPoint,:) = [1 0 0];
    elseif(abs(dataL2_single.fineH(dataL2_single.index_inside(iPoint))-mean_lake_height)<min_error)
        dataL2_single.c(iPoint,:) = [0 1 0];
    else
        dataL2_single.c(iPoint,:) = [(abs(dataL2_single.fineH(dataL2_single.index_inside(iPoint))-mean_lake_height)/max_error) 1-(abs(dataL2_single.fineH(dataL2_single.index_inside(iPoint))-mean_lake_height)/max_error) 0];
    end
    
end

for i_surf=1:length(dataL2_single.index_inside)
%     subplot(2,3,1); hold on
%     scatter(dataL1_single.lon(dataL2_single.index_inside(i_surf)),dataL1_single.lat(dataL2_single.index_inside(i_surf)),dataL2_single.s(dataL2_single.index_inside(i_surf)),dataL2_single.c(i_surf,:),'fill','LineWidth',1);
    subplot(2,3,1); hold on
    scatter(dataL2_single.lon(dataL2_single.index_inside(i_surf)),dataL2_single.lat(dataL2_single.index_inside(i_surf)),dataL2_single.s(dataL2_single.index_inside(i_surf)),dataL2_single.c(i_surf,:),'fill','LineWidth',1);
    set(gca,'YLim',R.Latlim);
    set(gca,'XLim',R.Lonlim);
    subplot(2,3,4);
    plot(dataL1_single.Wvf(dataL2_single.index_inside(i_surf),:));set(gca,'XLim',[1 512]);
    figlabels('Samples','L1 Power Waveform','','',12);
    subplot(2,3,5);
    plot(dataL1_single.PhD(dataL2_single.index_inside(i_surf),:).*1e-6*180/pi);set(gca,'XLim',[1 512]);
    figlabels('Samples','L1 Phase Difference','','',12);
    subplot(2,3,6);
    plot(dataL1_single.Coh(dataL2_single.index_inside(i_surf),:));set(gca,'XLim',[1 512]);
    figlabels('Samples','L1 Coherence','','',12);
    subplot(2,3,2:3);
    scatter(dataL2_single.lat(dataL2_single.index_inside(1:i_surf)),dataL2_single.fineH(dataL2_single.index_inside(1:i_surf)),dataL2_single.s(dataL2_single.index_inside(1:i_surf)),dataL2_single.c(1:i_surf,:),'fill','LineWidth',1);
    figlabels('Latitude','Lake Height L2','','',12);
    set(gca,'XLim',[min(dataL2_single.lat(dataL2_single.index_inside)) max(dataL2_single.lat(dataL2_single.index_inside))]);
    set(gca,'YLim',[min(dataL2_single.fineH(dataL2_single.index_inside)) max(dataL2_single.fineH(dataL2_single.index_inside))]);
    saveas (h, ['./results/CR2 L2 vs L1 vs L1BS/',dataL2_single.date,'_',num2str(i_surf, '%03d'),'_L1vsL2.jpg']);
    subplot(2,3,1);
    set(gca,'YLim',[dataL2_single.lat(dataL2_single.index_inside(i_surf))+latMargin(1) dataL2_single.lat(dataL2_single.index_inside(i_surf))-latMargin(1)]);
    set(gca,'XLim',[dataL2_single.lon(dataL2_single.index_inside(i_surf))-lonMargin(1) dataL2_single.lon(dataL2_single.index_inside(i_surf))+lonMargin(1)]);
    saveas (h, ['./results/CR2 L2 vs L1 vs L1BS/',dataL2_single.date,'_',num2str(i_surf, '%03d'),'_L1vsL2_z.jpg']);
end