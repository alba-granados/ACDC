Color=demcmap([min(alt_surf) max(alt_surf)],length(alt_surf));
alt_aux=alt_surf;
[~,color_indexes]=sort(alt_aux);

figure; scatter(lat(color_indexes),alt_aux(color_indexes),90,Color,'fill','MarkerEdgeColor','k')
hold all;
scatter(lat(loss_track-1),alt_surf(loss_track-1),80, 'fill','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0 0])
scatter(lat(loss_track(1:end-1)+1),alt_surf(loss_track(1:end-1)+1),80, 'fill','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 1 0])
set(gca,'YLim',[min(alt_surf) max(alt_surf)]);
set(gca,'XLim',[min(lat(lat~=0)) max(lat)]);
set(gca,'XLim',[35.6 36.4]);

figlabels('Latitude [degrees]','Elevation [m]','', 'Tracking window' ,12);

lla2kml_demcmap('tracking_window.kml',lat,lon,alt_surf,'.');

                        