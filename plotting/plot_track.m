function [product_mask] = plot_track(lat,lon,alt,options,product_mask)
global mask_flag
    
    WaterColor = [143 226 255]/255;
    degrees=max(lat)-min(lat);
%     axes(options.axes)
%     options.axes=axesm('MapProjection','lambert','MapLatLimit',[min(lat) max(lat)], 'MapLonLimit',[min(lon)-degrees/2 max(lon)+degrees/2]);
    options.axes = worldmap([min(lat) max(lat)],[min(lon)-degrees/2 max(lon)+degrees/2]);
    npoints = 4;
    lon_length = max(lon)+degrees-min(lon);
    lat_length = max(lat)-min(lat);
    init_lon = min(lon)- degrees/2;
    end_lon  = max(lon)+ degrees/2;
    lat_axis = min(lat):lat_length/npoints:max(lat);
    lon_axis = init_lon:lon_length/npoints:end_lon;
    setm(options.axes,'mlabellocation',lon_axis);
    setm(options.axes,'plabellocation',lat_axis);
    setm(options.axes,'Grid','on');
    p = findobj(options.axes,'type','patch'); % Find background
    set(p,'FaceColor',WaterColor); % Change background to white
    set(p,'LineWidth',1); % Change border width
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(options.axes, land, 'FaceColor', [244/255 164/255 96/255])
    lakes = shaperead('worldlakes', 'UseGeoCoords', true);
    geoshow(lakes, 'FaceColor', 'blue');
    rivers = shaperead('worldrivers', 'UseGeoCoords', true);
    geoshow(rivers, 'Color', 'blue');
    cities = shaperead('worldcities', 'UseGeoCoords', true);
    geoshow(cities, 'Marker', '.', 'Color', 'k');
    plotm(lat,lon,alt,'Parent',options.axes);
    figlabels('Longitude','Latitude','','Ground Track',10,options.axes);

    
    plot(alt, lat,'Parent',options.wd_axes);
    hold(options.wd_axes,'on');
    set(options.wd_axes,'YLim',[min(lat) max(lat)]);
    set(options.wd_axes,'XLim',[min(alt) max(alt)]);

    figlabels('Elevation [m]','Latitude','','Window Elevation',10,options.wd_axes);
    set(options.wd_axes,'YAxisLocation','right');
    
    if(options.draw_mask_check.Value)
       product_mask = [];
       mask_flag = 1;
       set(options.text_display,'String','[INFO]: Draw a rectangle on the map to select the records you want to process');
       set(options.text_display,'ForegroundColor',[0 1 0]);
       rect = getrect(options.axes);       
       [lat_mask_min,lon_mask_min] = minvtran(rect(1),rect(2));
       [lat_mask_max,lon_mask_max] = minvtran(rect(1)+rect(3),rect(2)+rect(4));
       product_mask.coord = [lon_mask_min lon_mask_max lon_mask_max lon_mask_min ; lat_mask_min lat_mask_min lat_mask_max lat_mask_max].';
       rect_draw=plotm(               [lat_mask_min lat_mask_min lat_mask_max lat_mask_max lat_mask_min],[lon_mask_min lon_mask_max lon_mask_max lon_mask_min lon_mask_min],'ko-','Parent',options.axes);
       index_inside = find(inpolygon(lon,lat,product_mask.coord(:,1),product_mask.coord(:,2))); % records inside the given mask
						
        if(isempty(index_inside))
            
            set(options.text_display,'String','[WARNING]: No records inside the mask, draw it again');
            set(options.text_display,'ForegroundColor',[1 1 0]);
            rect = getrect(options.axes);
            delete(rect_draw);
            [lat_mask_min,lon_mask_min] = minvtran(rect(1),rect(2));
            [lat_mask_max,lon_mask_max] = minvtran(rect(1)+rect(3),rect(2)+rect(4));
            product_mask.coord = [lon_mask_min lon_mask_max lon_mask_max lon_mask_min ; lat_mask_min lat_mask_min lat_mask_max lat_mask_max].';
            rect_draw = plotm(               [lat_mask_min lat_mask_min lat_mask_max lat_mask_max lat_mask_min],[lon_mask_min lon_mask_max lon_mask_max lon_mask_min lon_mask_min],'ko-','Parent',options.axes);
       
            index_inside = find(inpolygon(lon,lat,product_mask.coord(:,1),product_mask.coord(:,2))); % records inside the given mask
            if(isempty(index_inside))
                set(options.text_display,'String','[ERROR]: 2nd chance missed');
                set(options.text_display,'ForegroundColor',[1 0 0]);
            else
                 
                set(options.text_display,'String','[INFO]: Well done! Filtered records are being porcessed');
            end
        else
            set(options.text_display,'String','[INFO]: Well done! Filtered records are being porcessed');
               
        end
    end
    lat_length=max(product_mask.coord(:,2))-min(product_mask.coord(:,2));
    lon_length=max(product_mask.coord(:,1))-min(product_mask.coord(:,1));
    
    if(lat_length<lon_length)
        degrees = lon_length;
        options.axes=worldmap([min(product_mask.coord(:,2))-degrees/2 max(product_mask.coord(:,2))+degrees/2],[min(product_mask.coord(:,1)) max(product_mask.coord(:,1))]);
        lat_axis = min(product_mask.coord(:,2))-degrees/2:lat_length/npoints:max(product_mask.coord(:,2)+degrees/2);
        lon_axis = min(product_mask.coord(:,1)):lon_length/npoints:max(product_mask.coord(:,1));
    else
        degrees = lat_length;
        options.axes=worldmap([min(lat_axis) max(lat_axis)],[min(lon_axis)-degrees/2 max(lon_axis)+degrees/2]);
        lat_axis = min(product_mask.coord(:,2)):lat_length/npoints:max(product_mask.coord(:,2));
        lon_axis = min(product_mask.coord(:,1))-degrees/2:lon_length/npoints:max(product_mask.coord(:,1)+degrees/2);

    end
    
    
    
    
    
    
%     setm(options.axes,'mlabellocation',lon_axis);
%     setm(options.axes,'plabellocation',lat_axis);
    p = findobj(options.axes,'type','patch'); % Find background
    set(p,'FaceColor',WaterColor); % Change background to white
    set(p,'LineWidth',1); % Change border width
    land = shaperead('landareas', 'UseGeoCoords', true);
    geoshow(options.axes, land, 'FaceColor', [244/255 164/255 96/255])
    lakes = shaperead('worldlakes', 'UseGeoCoords', true);
    geoshow(lakes, 'FaceColor', 'blue');
    rivers = shaperead('worldrivers', 'UseGeoCoords', true);
    geoshow(rivers, 'Color', 'blue');
    cities = shaperead('worldcities', 'UseGeoCoords', true);
    geoshow(cities, 'Marker', '.', 'Color', 'k');
    plotm(lat,lon,alt,'Parent',options.axes);
    figlabels('Longitude','Latitude','','Ground Track',10,options.axes);
    
    set(options.wd_axes,'YLim',[min(product_mask.coord(:,2)) max(product_mask.coord(:,2))]);
        
    
end