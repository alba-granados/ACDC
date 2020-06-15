% Objective compute the geophysical correctons to be applied to the measured range
% v1.0 First version
% geophysical_corrections_sea_ice  = dry_tropo_correction_bursts+wet_tropo_correction_bursts+ inverse_baro_correction_bursts+           GIM_iono_correction_bursts+ ocean_equilibrium_tide_bursts+  long_period_tide_height_bursts+   ocean_loading_tide_bursts+ solid_earth_tide_bursts+ geocentric_polar_tide_bursts;
% geophysical_corrections_land_ice = dry_tropo_correction_bursts+wet_tropo_correction_bursts+                                                                         GIM_iono_correction_bursts+                                                                                                                               ocean_loading_tide_bursts+ solid_earth_tide_bursts+ geocentric_polar_tide_bursts;
% geophysical_corrections_ocean    = dry_tropo_correction_bursts+wet_tropo_correction_bursts+ Dynamic_atmospheric_correction_bursts+    GIM_iono_correction_bursts+ ocean_equilibrium_tide_bursts+  long_period_tide_height_bursts+   ocean_loading_tide_bursts+ solid_earth_tide_bursts+ geocentric_polar_tide_bursts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As done in CR2>
% geophysical_corrections_land_ice =    dry_tropo_correction_bursts+
%                                       wet_tropo_correction_bursts+
%                                       GIM_iono_correction_bursts+
%                                       ocean_loading_tide_bursts+
%                                       solid_earth_tide_bursts+ 
%                                       geocentric_polar_tide_bursts;
% 1. Dry Tropo
% 2. Wet Tropo
% 3. Iono Model
% 4. Ocean Loading
% 5. Solid Earth Tide
% 6. Polar Tide

function [solid_earth, geocentric_tide, ocean_loading,dry_tropo_correction,wet_tropo_correction,iono_correction]=compute_L2_corrections(L2,pos,method)
time_1Hz = L2.data.time_01;
time_20Hz = L2.data.time_20_ku;
dry_tropo_correction = double(L2.data.mod_dry_tropo_cor_zero_altitude_01).*double(L2.attributes.mod_dry_tropo_cor_zero_altitude_01.scale_factor);
% model dry tropospheric correction at zero altitude: 1 Hz
wet_tropo_correction=double(L2.data.rad_wet_tropo_cor_01_ku).*double(L2.attributes.rad_wet_tropo_cor_01_ku.scale_factor);
%radiometer wet tropospheric correction: 1 Hz Ku band
iono_correction=double(L2.data.iono_cor_alt_01_ku).*double(L2.attributes.iono_cor_alt_01_ku.scale_factor);
%altimeter ionospheric correction: 1 Hz Ku band
ocean_loading_tide=double(L2.data.load_tide_sol1_01).*L2.attributes.ocean_tide_sol1_01.scale_factor;
%Solution 1 corresponds to GOT00.2 model. Includes the corresponding loading tide (load_tide_sol1_01) and equilibrium long-period ocean tide height (ocean_tide_eq_01). The permanent tide (zero frequency) is not included in this parameter because it is included in the geoid and mean sea surface (geoid_01, mean_sea_surf_sol1_01).
solid_earth_tide=double(L2.data.solid_earth_tide_01).*double(L2.attributes.solid_earth_tide_01.scale_factor);
%'Cartwright and Edden [1973] Corrected tables of tidal harmonics - J. Geophys. J. R. Astr. Soc., 33, 253-264.'
geocentric_polar_tides=double(L2.data.pole_tide_01).*double(L2.attributes.pole_tide_01.scale_factor);
%'Wahr [1985] Deformation of the Earth induced by polar motion - J. Geophys. Res. (Solid Earth), 90, 9363-9368.'

%interpolate fill values
valid_corr=(dry_tropo_correction~=double(L2.attributes.mod_dry_tropo_cor_zero_altitude_01.FillValue)*1e-4);
dry_tropo_correction = interp1(time_1Hz(valid_corr),dry_tropo_correction(valid_corr),time_1Hz,method);
valid_corr=(wet_tropo_correction~=double(L2.attributes.rad_wet_tropo_cor_01_ku.FillValue)*1e-4);
wet_tropo_correction = interp1(time_1Hz(valid_corr),wet_tropo_correction(valid_corr),time_1Hz,method);
valid_corr=(iono_correction~=double(L2.attributes.iono_cor_alt_01_ku.FillValue)*1e-4);
iono_correction = interp1(time_1Hz(valid_corr),iono_correction(valid_corr),time_1Hz,method);
valid_corr=(ocean_loading_tide~=double(L2.attributes.load_tide_sol1_01.FillValue)*1e-4);
ocean_loading_tide = interp1(time_1Hz(valid_corr),ocean_loading_tide(valid_corr),time_1Hz,method);
valid_corr=(solid_earth_tide~=double(L2.attributes.solid_earth_tide_01.FillValue)*1e-4);
solid_earth_tide = interp1(time_1Hz(valid_corr),solid_earth_tide(valid_corr),time_1Hz,method);
valid_corr=(geocentric_polar_tides~=double(L2.attributes.pole_tide_01.FillValue)*1e-4);
geocentric_polar_tides = interp1(time_1Hz(valid_corr),geocentric_polar_tides(valid_corr),time_1Hz,method);




geophysical_corrections = dry_tropo_correction + wet_tropo_correction + iono_correction + iono_correction + ocean_loading_tide +solid_earth_tide+geocentric_polar_tides;
geophysical_correction = geophysical_corrections(pos);
% figure; plot(dry_tropo_correction); hold all; plot(pos,dry_tropo_correction(pos),'ro');
% plot(wet_tropo_correction); plot(pos,wet_tropo_correction(pos),'ro');
% plot(iono_correction); plot(pos,iono_correction(pos),'ro');
% plot(ocean_loading_tide); plot(pos,ocean_loading_tide(pos),'ro');
% plot(solid_earth_tide);plot(pos,solid_earth_tide(pos),'ro');
% plot(geocentric_polar_tides);plot(pos,geocentric_polar_tides(pos),'ro');
% figlabels('records','correction [meters]','',['Geophysical correction TRP record:' num2str(geophysical_corrections(pos)) 'meters'],12);
% disp(['Total=' num2str(geophysical_corrections(pos))...
%     '\ Dry = ' num2str(dry_tropo_correction(pos)) ...
%     '\ Wet = ' num2str(wet_tropo_correction(pos)) ...
%     '\ Iono = ' num2str(iono_correction(pos)) ...
%     '\ Solid Earth = ' num2str(solid_earth_tide(pos)) ...
%     '\ Geocentric = ' num2str(geocentric_polar_tides(pos)) ...
%     '\ Ocean Loading = ' num2str(ocean_loading_tide(pos)) ...
%     ]);

solid_earth     =solid_earth_tide(pos);
geocentric_tide	=geocentric_polar_tides(pos);
ocean_loading   =ocean_loading_tide(pos);

end
