%Geoid

for kk=1:1:361
        latitude(kk,1:361) = 90:-0.5:-90;
        longitude(kk,1:361) = kk-1;
        h(kk,1:361) = geoidheight(latitude(kk,1:361),longitude(kk,1:361),'EGM2008');
    end
    load coast
    figure; k=surf((180:-1:-180),90:-0.5:-90,h'); colormap('jet'); set(k, 'edgecolor','none');grid on; colorbar;
    geoshow(flipud(lat),flipud(long),'DisplayType','polygon','FaceColor','none')
    figlabels('Longitude','Latitude','','Geoid EGM 2008',16);