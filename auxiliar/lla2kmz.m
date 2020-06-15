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
%   lla2kmz(filename,lat,long,alt)
%
% Inputs 
% Filename, Latitude,Longitude, Altitude (optional) 
%   
% Outputs
%   
% Comments:
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




function lla2kml(filename,lat,long,alt)


if(length(lat)~=length(long))
    disp('Latitude and longitude with differents lengths');
   return; 
end

if(~alt)
    alt=zeros(1,length(lat)); % initialize altitude to 0;
end

%% Save Locations in KML format
fidKML = fopen(sprintf('%s.kml', filename) ,'w');
% headers in 1st line
outText(fidKML,'<?xml version="1.0" standalone="yes"?>');
outText(fidKML,'<kml xmlns="http://earth.google.com/kml/2.2">');
outText(fidKML,'   <Document>');
outText(fidKML,'       <Folder id="Waypoints">');

   
    for iPoint = 1:length(lat)
            
                % Estimate color to represent beacon on a map according localization error
            if( alt < 5 )     % 5 km
                Color = 'FF008000'; %green      
                IconColor = 'FF008000';
            elseif( 5 <= alt && alt < 10 ) % between 5 and 10 km
                Color = 'FF00FFFF';%yellow
                IconColor = 'FF00FFFF';
            elseif( 10 <= alt && alt < 20 ) % between 5 and 10 km
                Color = 'FF00A5FF'; %orange
                IconColor = 'FF00A5FF';
            elseif( 20 <= alt && alt < 120 ) % between 10 and 30 km
                Color = 'FF0000FF'; %red
                IconColor = 'FF0000FF';
            else
                Color = 'FF000000';
                IconColor = 'FF222222';
            end
            DtoPrint = sprintf('%.1f km',beacMeanErr2D);

            %Print points in KML
            outText(fidKML,'           <Placemark>');
            outText(fidKML,'               <Point>');
            outText(fidKML,'                   <altitudeMode>clampToGround</altitudeMode>');
            outText(fidKML,sprintf('                   <coordinates>%f,%f,%f</coordinates>', lon, lat, alt));
            outText(fidKML,'               </Point>');
            outText(fidKML,'               <Style>');
            outText(fidKML,'                   <IconStyle>');
            outText(fidKML,sprintf('                       <color>%s</color>',IconColor));
            outText(fidKML,'                   </IconStyle>');
            outText(fidKML,'                   <LabelStyle>');
            outText(fidKML,sprintf('                       <color>%s</color>',Color));
            outText(fidKML,'                   </LabelStyle>');      
            outText(fidKML,'               </Style>');
            outText(fidKML,sprintf('               <description><![CDATA[%s]]></description>',DtoPrint));
            outText(fidKML,sprintf('               <name>%i</name>',iPoint));
            outText(fidKML,'               <styleUrl>#gv_waypoint</styleUrl>');
            outText(fidKML,'           </Placemark>');
    end

    
    %Print property bullets in KML
    outText(fidKML,'       <name>Waypoints</name>');
    outText(fidKML,'       <visibility>1</visibility>');

    outText(fidKML,'   </Folder>');
    outText(fidKML,'   <Snippet></Snippet>');


    outText(fidKML,'   <Style id="gv_waypoint_normal">');
    outText(fidKML,'       <BalloonStyle>');
    outText(fidKML,'           <text><![CDATA[<p align="left" style="white-space:nowrap;"><font size="+1"><b>$[name]</b></font></p> <p align="left">$[description]</p>]]></text>');
    outText(fidKML,'       </BalloonStyle>');
    outText(fidKML,'       <IconStyle>');
    outText(fidKML,'           <Icon>');
    outText(fidKML,'               <href>http://maps.google.ca/mapfiles/kml/pal4/icon57.png</href>');
    outText(fidKML,'           </Icon>');
    outText(fidKML,'           <color>FFFFFFFF</color>');
    outText(fidKML,'           <hotSpot x="0.5" xunits="fraction" y="0.5" yunits="fraction" />');
    outText(fidKML,'       </IconStyle>');
    outText(fidKML,'       <LabelStyle>');
    outText(fidKML,'           <color>FFFFFFFF</color>');
    outText(fidKML,'           <scale>0</scale>');
    outText(fidKML,'       </LabelStyle>');
    outText(fidKML,'   </Style>');
    outText(fidKML,'   <Style id="gv_waypoint_highlight">');
    outText(fidKML,'       <BalloonStyle>');
    outText(fidKML,'           <text><![CDATA[<p align="left" style="white-space:nowrap;"><font size="+1"><b>$[name]</b></font></p> <p align="left">$[description]</p>]]></text>');
    outText(fidKML,'       </BalloonStyle>');
    outText(fidKML,'       <IconStyle>');
    outText(fidKML,'           <Icon>');
    outText(fidKML,'               <href>http://maps.google.ca/mapfiles/kml/pal4/icon57.png</href>');
    outText(fidKML,'           </Icon>');
    outText(fidKML,'           <color>FFFFFFFF</color>');
    outText(fidKML,'           <hotSpot x="0.5" xunits="fraction" y="0.5" yunits="fraction" />');
    outText(fidKML,'           <scale>1.2</scale>');
    outText(fidKML,'       </IconStyle>');
    outText(fidKML,'       <LabelStyle>');
    outText(fidKML,'           <color>FFFFFFFF</color>');
    outText(fidKML,'           <scale>1</scale>');
    outText(fidKML,'       </LabelStyle>');
    outText(fidKML,'   </Style>');
    outText(fidKML,'   <StyleMap id="gv_waypoint">');
    outText(fidKML,'       <Pair>');
    outText(fidKML,'           <key>normal</key>');
    outText(fidKML,'           <styleUrl>#gv_waypoint_normal</styleUrl>');
    outText(fidKML,'       </Pair>');
    outText(fidKML,'       <Pair>');
    outText(fidKML,'           <key>highlight</key>');
    outText(fidKML,'           <styleUrl>#gv_waypoint_highlight</styleUrl>');
    outText(fidKML,'       </Pair>');
    outText(fidKML,'   </StyleMap>');
    outText(fidKML,'   <StyleMap id="gv_trackpoint">');
    outText(fidKML,'       <Pair>');
    outText(fidKML,'           <key>normal</key>');
    outText(fidKML,'           <styleUrl>#gv_trackpoint_normal</styleUrl>');
    outText(fidKML,'       </Pair>');
    outText(fidKML,'       <Pair>');
    outText(fidKML,'           <key>highlight</key>');
    outText(fidKML,'           <styleUrl>#gv_trackpoint_highlight</styleUrl>');
    outText(fidKML,'       </Pair>');
    outText(fidKML,'   </StyleMap>');

    outText(fidKML,'   <name><![CDATA[GPS data]]></name>');
    outText(fidKML,'   <open>0</open>');
    outText(fidKML,'   <visibility>1</visibility>');
    outText(fidKML,'   </Document>');
    outText(fidKML,'</kml>');
    
    
    fclose(fidKML);
      
       
end