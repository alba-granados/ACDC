%%





%%
function lla2kml_tour(filename,L1BS_buff)
global currenfolder

[~,L1BS_size] = size(L1BS_buff);

fidKML = fopen(sprintf('%s.kml', filename) ,'w');

fprintf(fidKML,'<?xml version="1.0" encoding="UTF-8"?>');
fprintf(fidKML,'<kml xmlns="http://www.opengis.net/kml/2.2" xmlns:gx="http://www.google.com/kml/ext/2.2" xmlns:kml="http://www.opengis.net/kml/2.2" xmlns:atom="http://www.w3.org/2005/Atom">');
fprintf(fidKML,'<gx:Tour>');

fprintf(fidKML,['	<name>' currenfolder '_pass_tour</name>']);

fprintf(fidKML,'	<gx:Playlist>');
fprintf(fidKML,'		<gx:FlyTo>');
fprintf(fidKML,'			<LookAt>');
fprintf(fidKML,'				<gx:horizFov>60</gx:horizFov>');

fprintf(fidKML,['				<longitude>' num2str(L1BS_buff(1).lon_surf) '</longitude>']);
fprintf(fidKML,['				<latitude>' num2str(L1BS_buff(1).lat_surf) '</latitude>']);

fprintf(fidKML,'				<altitude>0</altitude>');
fprintf(fidKML,'				<heading>13</heading>');
fprintf(fidKML,'				<tilt>50</tilt>');
fprintf(fidKML,'				<range>14000</range>');
fprintf(fidKML,'				<gx:altitudeMode>absolute</gx:altitudeMode>');
fprintf(fidKML,'			</LookAt>');
fprintf(fidKML,'		</gx:FlyTo>');
fprintf(fidKML,'		<gx:FlyTo>');

fprintf(fidKML,['			<gx:duration>' num2str(L1BS_size/10) '</gx:duration>']);

fprintf(fidKML,'			<gx:flyToMode>smooth</gx:flyToMode>');
fprintf(fidKML,'			<LookAt>');
fprintf(fidKML,'				<gx:horizFov>60</gx:horizFov>');

fprintf(fidKML,['				<longitude>' num2str(L1BS_buff(end).lon_surf) '</longitude>']);
fprintf(fidKML,['				<latitude>' num2str(L1BS_buff(end).lat_surf) '</latitude>']);


fprintf(fidKML,'				<altitude>0</altitude>');
fprintf(fidKML,'				<heading>13</heading>');
fprintf(fidKML,'				<tilt>50</tilt>');
fprintf(fidKML,'				<range>14000</range>');
fprintf(fidKML,'				<gx:altitudeMode>absolute</gx:altitudeMode>');
fprintf(fidKML,'			</LookAt>');
fprintf(fidKML,'		</gx:FlyTo>');
fprintf(fidKML,'	</gx:Playlist>');
fprintf(fidKML,'</gx:Tour>');
fprintf(fidKML,'</kml>');