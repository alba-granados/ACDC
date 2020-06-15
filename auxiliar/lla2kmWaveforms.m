%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% CryoSat 2 calibration over transponders
%
% ---------------------------------------------------------
%
% Calling
%
%   lla2kml(filename,lat,long,alt, mode, imagenames)
%
% Inputs 
% Filename, Latitude,Longitude, Altitude (optional) 
% imagesPath: Paths of the images
% description
% Outputs
%   
% Comments:

% Two modes of printing,
% mode '.' ---> prints icons
% mode '-' ---> prints a line
%   
% ----------------------------------------------------------
% 
% Author:   Albert García / isardSAT
%
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Mònica Roca / isardSAT (XX/XX/11)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function lla2kmlWaveforms(filename,lat,lon,alt, mode, imagesPath)



if(length(lat)~=length(lon))
    disp('Latitude and longitude with differents lengths');
   return; 
end

if(~alt)
    alt=zeros(1,length(lat)); % initialize altitude to 0;
end

if(~mode)
   mode='.';
end
if(~imagesPath)
   disp('No Path Waveforms ');
   return;
end
WaveformsPath   =   dir(imagesPath);

%% Save Locations in KML format
fidKML = fopen(sprintf('%s.kml', filename) ,'w');
% headers in 1st line
fprintf(fidKML,'<?xml version="1.0" standalone="yes"?>\n');
fprintf(fidKML,'<kml xmlns="http://earth.google.com/kml/2.2">\n');
fprintf(fidKML,'   <Document>\n');


if strcmp(mode,'-')  
    fidTXT = fopen(sprintf('%s.csv', filename) ,'r');
    
    fprintf(fidKML,'   <Style id="multiTrackB_h">\n');
    fprintf(fidKML,' <IconStyle>\n');
    fprintf(fidKML,'   <color>ffffccff</color>\n');
    fprintf(fidKML,'   <scale>1.2</scale>\n');
    fprintf(fidKML,'   <Icon>\n');
	fprintf(fidKML,' <href>http://earth.google.com/images/kml-icons/track-directional/track-0.png</href>\n');
    fprintf(fidKML,'   </Icon>\n');
    fprintf(fidKML,' </IconStyle>\n');
    fprintf(fidKML,' <LineStyle>\n');
    fprintf(fidKML,'   <color>99ffcc00</color>\n');
    fprintf(fidKML,'   <width>8</width>\n');
    fprintf(fidKML,' </LineStyle>\n');
    fprintf(fidKML,'  </Style>\n');
    fprintf(fidKML,'  <Style id="multiTrackB_n">\n');
    fprintf(fidKML,' <IconStyle>\n');
    fprintf(fidKML,'   <color>ffffcc00</color>\n');
    fprintf(fidKML,'   <Icon>\n');
	fprintf(fidKML,' <href>http://earth.google.com/images/kml-icons/track-directional/track-0.png</href>\n');
    fprintf(fidKML,'   </Icon>\n');
    fprintf(fidKML,' </IconStyle>\n');
    fprintf(fidKML,' <LineStyle>\n');
    fprintf(fidKML,'   <color>9999cc00</color>\n');
    fprintf(fidKML,'   <width>6</width>\n');
    fprintf(fidKML,' </LineStyle>\n');
    fprintf(fidKML,' </Style>\n');
    fprintf(fidKML,' <StyleMap id="multiTrackB">\n');
    fprintf(fidKML,'   <Pair>\n');
	fprintf(fidKML,' <key>normal</key>\n');
	fprintf(fidKML,' <styleUrl>#multiTrackB_n</styleUrl>\n');
    fprintf(fidKML,'   </Pair>\n');
    fprintf(fidKML,'   <Pair>\n');
	fprintf(fidKML,' <key>highlight</key>\n');
   fprintf(fidKML,'      <styleUrl>#multiTrackB_h</styleUrl>\n');
   fprintf(fidKML,'    </Pair>\n');
   fprintf(fidKML,' </StyleMap> \n');
   fprintf(fidKML,'       <Folder id="Line">\n');  

%     fprintf(fidKML,'<ExtendedData>\n');
%     fprintf(fidKML,' <Data name="ANX Longitude">\n');
%     fprintf(fidKML,' <value>12.313000 deg</value>\n');
%     fprintf(fidKML,'</Data>\n');
%     fprintf(fidKML,'</ExtendedData>\n');
    orbit=0;
    endPoint=0; 
    while(~feof(fidTXT))

        ascLine = fgetl(fidTXT);
        separator = strfind(ascLine,';');
        if length(separator)==1 && endPoint==0;
                  fprintf(fidKML,'</coordinates>\n');
                  fprintf(fidKML,'</LineString>\n');    
                  fprintf(fidKML,'</Placemark>\n');
                  endPoint=1;
             
        elseif length(separator)==2
            if endPoint==1 || orbit==0
                orbit=orbit+1;
                fprintf(fidKML,'<Placemark>\n');
                fprintf(fidKML,'<name> Orbit Ground Track</name>\n');
                fprintf(fidKML,' <styleUrl>#multiTrackB</styleUrl>\n');
                fprintf(fidKML,' <LineString>\n');
                fprintf(fidKML,'<tessellate>1</tessellate>\n');
                fprintf(fidKML,'<altitudeMode>clampToGround</altitudeMode>\n');
                fprintf(fidKML,'<coordinates>');
                
            end
            lat=str2double(ascLine(1:separator(1)-1));
            lon=str2double(ascLine(separator(1)+1:separator(2)-1));
            fprintf(fidKML,'%d,%d,%d\n',lon,lat,0);
            endPoint=0;
            
        elseif length(separator)==3
            if endPoint==1 || orbit==0
                orbit=orbit+1;
                fprintf(fidKML,'<Placemark>\n');
                fprintf(fidKML,'<name> Orbit Ground Track</name>\n');
                fprintf(fidKML,' <styleUrl>#multiTrackB</styleUrl>\n');
                fprintf(fidKML,' <LineString>\n');
                fprintf(fidKML,'<tessellate>1</tessellate>\n');
                fprintf(fidKML,'<altitudeMode>clampToGround</altitudeMode>\n');
                fprintf(fidKML,'<coordinates>');
            end
            lat=str2double(ascLine(1:separator(1)-1));
            lon=str2double(ascLine(separator(1)+1:separator(2)-1));
            alt=str2double(ascLine(separator(2)+1:separator(3)-1));
            fprintf(fidKML,',%d,%d,%dºn',lon,lat,alt);
            endPoint=0;
            
        end
    end
    fprintf(fidKML,',</coordinates>\n');
    fprintf(fidKML,'</LineString>\n');    
    fprintf(fidKML,'</Placemark>\n');
    fprintf(fidKML,'</Folder>\n');
     fclose(fidTXT);
else
    fprintf(fidKML,'       <Folder id="Waypoints">\n');
    for iPoint = 1:length(lat)
            
            % Progress bar
%             customText = 'Converting lla to kml...';
%             percentageDone =  iPoint / length(lat);
%             stopBar= progressbar(percentageDone, 0, customText);
%             if (stopBar) 
%                 break; 
%             end
        
                % Estimate color to represent point on a map according
                % error  RGB ----> FF+"BGR"
            if( alt(iPoint) < 20 )     % 5 km
                Color = 'FF008000'; %green      
                IconColor = 'FF008000';
            elseif( 20 <= alt(iPoint) && alt(iPoint) < 250 ) % between 5 and 10 km
                Color = 'FFFFFF00';%blue
                IconColor = 'FFFFFF00';
            elseif( 250 <= alt(iPoint) && alt(iPoint) < 500 ) % between 5 and 10 km
                Color = 'FF00A5FF'; %orange
                IconColor = 'FF00A5FF';
            elseif( 500 <= alt(iPoint) && alt(iPoint) < 1000 ) % between 10 and 30 km
                Color = 'FF0000FF'; %red
                IconColor = 'FF0000FF';
            else
                Color = 'FF000000';
                IconColor = 'FF222222';
            end
            DtoPrint = sprintf('%.1f meters',alt(iPoint));

            %Print points in KML
            fprintf(fidKML,'           <Placemark>\n');
            fprintf(fidKML,'               <Point>\n');
            fprintf(fidKML,'                   <altitudeMode>relativeToGround</altitudeMode>\n');
            fprintf(fidKML,sprintf('                   <coordinates>%f,%f,%f</coordinates>\n', lon(iPoint), lat(iPoint), alt(iPoint)));
            fprintf(fidKML,'               </Point>\n');
            fprintf(fidKML,'               <Style>\n');
            fprintf(fidKML,'                   <IconStyle>');
            fprintf(fidKML,sprintf('                       <color>%s</color>\n',IconColor));
            fprintf(fidKML,'                   </IconStyle>\n');
            fprintf(fidKML,'                   <LabelStyle>\n');
            fprintf(fidKML,sprintf('                       <color>%s</color>\n',Color));
            fprintf(fidKML,'                   </LabelStyle>\n');      
            fprintf(fidKML,'               </Style>\n');
            fprintf(fidKML,sprintf('               <description><![CDATA[<img src="%s%s">Lat=%d; Lon=%d]]></description>\n',imagesPath,WaveformsPath(iPoint+2).name,lat(iPoint),lon(iPoint)));
            fprintf(fidKML,sprintf('               <name>%i</name>\n',iPoint));
            fprintf(fidKML,'               <styleUrl>#gv_waypoint</styleUrl>\n');
            fprintf(fidKML,'           </Placemark>\n');
            
            

    end

    
    %Print property bullets in KML
    fprintf(fidKML,'       <name>Waypoints</name>\n');
    fprintf(fidKML,'       <visibility>1</visibility>\n');

    fprintf(fidKML,'   </Folder>\n');
    fprintf(fidKML,'   <Snippet></Snippet>\n');


    fprintf(fidKML,'   <Style id="gv_waypoint_normal">\n');
    fprintf(fidKML,'       <BalloonStyle>');
    fprintf(fidKML,'           <text><![CDATA[<p align="left" style="white-space:nowrap;"><font size="+1"><b>$[name]</b></font></p> <p align="left">$[description]</p>]]></text>\n');
    fprintf(fidKML,'       </BalloonStyle>\n');
    fprintf(fidKML,'       <IconStyle>\n');
    fprintf(fidKML,'           <Icon>\n');
    fprintf(fidKML,'               <href>http://maps.google.ca/mapfiles/kml/pal4/icon57.png</href>\n');
    fprintf(fidKML,'           </Icon>\n');
    fprintf(fidKML,'           <color>FFFFFFFF</color>\n');
    fprintf(fidKML,'           <hotSpot x="0.5" xunits="fraction" y="0.5" yunits="fraction" />\n');
    fprintf(fidKML,'       </IconStyle>\n');
    fprintf(fidKML,'       <LabelStyle>\n');
    fprintf(fidKML,'           <color>FFFFFFFF</color>\n');
    fprintf(fidKML,'           <scale>0</scale>\n');
    fprintf(fidKML,'       </LabelStyle>\n');
    fprintf(fidKML,'   </Style>\n');
    fprintf(fidKML,'   <Style id="gv_waypoint_highlight">\n');
    fprintf(fidKML,'       <BalloonStyle>\n');
    fprintf(fidKML,'           <text><![CDATA[<p align="left" style="white-space:nowrap;"><font size="+1"><b>$[name]</b></font></p> <p align="left">$[description]</p>]]></text>\n');
    fprintf(fidKML,'       </BalloonStyle>\n');
    fprintf(fidKML,'       <IconStyle>\n');
    fprintf(fidKML,'           <Icon>\n');
    fprintf(fidKML,'               <href>http://maps.google.ca/mapfiles/kml/pal4/icon57.png</href>\n');
    fprintf(fidKML,'           </Icon>\n');
    fprintf(fidKML,'           <color>FFFFFFFF</color>\n');
    fprintf(fidKML,'           <hotSpot x="0.5" xunits="fraction" y="0.5" yunits="fraction" />\n');
    fprintf(fidKML,'           <scale>1.2</scale>\n');
    fprintf(fidKML,'       </IconStyle>\n');
    fprintf(fidKML,'       <LabelStyle>\n');
    fprintf(fidKML,'           <color>FFFFFFFF</color>\n');
    fprintf(fidKML,'           <scale>1</scale>\n');
    fprintf(fidKML,'       </LabelStyle>\n');
    fprintf(fidKML,'   </Style>\n');
    fprintf(fidKML,'   <StyleMap id="gv_waypoint">\n');
    fprintf(fidKML,'       <Pair>\n');
    fprintf(fidKML,'           <key>normal</key>\n');
    fprintf(fidKML,'           <styleUrl>#gv_waypoint_normal</styleUrl>\n');
    fprintf(fidKML,'       </Pair>\n');
    fprintf(fidKML,'       <Pair>\n');
    fprintf(fidKML,'           <key>highlight</key>\n');
    fprintf(fidKML,'           <styleUrl>#gv_waypoint_highlight</styleUrl>\n');
    fprintf(fidKML,'       </Pair>\n');
    fprintf(fidKML,'   </StyleMap>\n');
    fprintf(fidKML,'   <StyleMap id="gv_trackpoint">\n');
    fprintf(fidKML,'       <Pair>\n');
    fprintf(fidKML,'           <key>normal</key>\n');
    fprintf(fidKML,'           <styleUrl>#gv_trackpoint_normal</styleUrl>\n');
    fprintf(fidKML,'       </Pair>\n');
    fprintf(fidKML,'       <Pair>\n');
    fprintf(fidKML,'           <key>highlight</key>\n');
    fprintf(fidKML,'           <styleUrl>#gv_trackpoint_highlight</styleUrl>\n');
    fprintf(fidKML,'       </Pair>\n');
    fprintf(fidKML,'   </StyleMap>\n');

    fprintf(fidKML,'   <name><![CDATA[GPS data]]></name>\n');
   
      
    
end

    fprintf(fidKML,'   <open>0</open>\n');
    fprintf(fidKML,'   <visibility>1</visibility>\n');
    fprintf(fidKML,'   </Document>\n');
    fprintf(fidKML,'</kml>\n');
    fclose(fidKML);
    