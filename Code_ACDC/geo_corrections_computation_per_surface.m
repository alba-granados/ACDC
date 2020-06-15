function [res] = geo_corrections_computation_per_surface(L1A_buffer,L1BS,cnf_p_ACDC)
% -------------------------------------------------------------------------
% Created by isardSAT S.L. 
% -------------------------------------------------------------------------
% This code organzies and geophysical correctiosn to be applied to the
% retracked range based on the info available in the L1B product and
% depending on the surface being observed
% -------------------------------------------------------------------------
% 
% Author:           Eduard Makhoul / isardSAT
%
% Reviewer:         Mònica Roca / isardSAT
%
% Last revision:    Eduard Makhoul / isardSAT V1 03/10/2016
% This software is built with internal funding 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% INPUT:
%      MANDATORY
%       -data    =   structure of input data used in the L2 processing       
%       -cnf_p = configuration parameters structure for L2 processing
%       
% OUTPUT:
%       res        =   structure of Corrections to be included in the L2
%       product
% RESTRICTIONS: The type dependent surface corrections for each record
% independently is implemented only for CS2 data only for ocean (falg
% surface==0) and land ice (flag surface==2). Application of same
% geophysical corrections for all the records based on a single surface
% type is currently available only for "open ocean", "land ice" and "sea
% ice". No inclusion of the sea state bias corrections is available since
% there is no such available info at L1B.
% -------------------------------------------------------------------------  
% -------------------------------------------------------------------------
% Versions control:
% v1.0: 

global mission

%------------ COMPUTE THE BURST AT NADIR ----------------------------------
%for computation of the geophysical correction
%Compute nadir burst for a given surface
[~,beam_index_nadir]=min(abs(-pi/2+L1BS.beam_ang_surf'));
burst_index_nadir= L1BS.burst_index(beam_index_nadir);

surface_type_flag = L1A_buffer(burst_index_nadir).surface_type_flag_bursts;
dry_trop=L1A_buffer(burst_index_nadir).dry_tropo_correction_bursts;
wet_trop=L1A_buffer(burst_index_nadir).wet_tropo_correction_bursts;
inv_bar=L1A_buffer(burst_index_nadir).inverse_baro_correction_bursts;
dac=L1A_buffer(burst_index_nadir).Dynamic_atmospheric_correction_bursts;
GIM_iono_correction=L1A_buffer(burst_index_nadir).GIM_iono_correction_bursts;
model_iono_correction=L1A_buffer(burst_index_nadir).model_iono_correction_bursts;
ocean_equilibrium_tide=L1A_buffer(burst_index_nadir).ocean_equilibrium_tide_bursts;
ocean_longperiod_tide=L1A_buffer(burst_index_nadir).long_period_tide_height_bursts;
ocean_loading_tide=L1A_buffer(burst_index_nadir).ocean_loading_tide_bursts;
solidearth_tide=L1A_buffer(burst_index_nadir).solid_earth_tide_bursts;
geocentric_polar_tide=L1A_buffer(burst_index_nadir).geocentric_polar_tide_bursts;

%----------------- Define variables ---------------------------------------
res.total_geo_corr        =   0.0;
%------- corrections flags --------------------------------------------
% flags indicating the inclusion 1 or not 0 of the related correction
res.flags.dry_trop                  =  0;
res.flags.wet_trop                  =  0;
res.flags.inv_bar                   =  0;
res.flags.dac                       =  0;
res.flags.gim_ion                   =  0;
res.flags.model_ion                 =  0;
res.flags.ocean_equilibrium_tide    =  0;
res.flags.ocean_longperiod_tide     =  0;
res.flags.ocean_loading_tide        =  0;
res.flags.solidearth_tide           =  0;
res.flags.geocentric_polar_tide     =  0;


if cnf_p_ACDC.force_geocorr_surf_type
    %forcing all product to same corrections
    type_surface=cnf_p_ACDC.product_type_surface;
else
    %apply correction depending on the type of surface
    switch mission
        case {'CR2','CS2','S3'}
            switch surface_type_flag
                case 0
                    %open Ocean
                    type_surface='open_ocean';
                case 1
                    %Enclosed seas or lakes
                    type_surface='sea_ice';
                case 2
                    %Continental ice
                    type_surface='land_ice';
                case 3
                    %land/transponder
                    type_surface='land_ice';
            end
    end
end
%---------------------total as if all where same surface --------------
% Different contributions depending on surface extracted from CS2
% product handbook
% force the same group of corrections to all records
% independently of the surface flag included in the L1B product

% ---------- Checking iono corrections to be applied ------------------
if any(GIM_iono_correction)
    %from product handbook of Cryosat-2 GIM is nominal choice
    %and if not available use Bent model
    iono_corr=GIM_iono_correction;
    res.flags.gim_ion(:)=1;
else
    iono_corr=model_iono_correction;
    res.flags.model_ion(:)=1;
end

switch type_surface
    case 'open_ocean'
        res.total_geo_corr= ocean_equilibrium_tide+... %is referred in cryosat-2 as ocean tide or elastic ocean tide
            ocean_longperiod_tide+...
            ocean_loading_tide+...
            solidearth_tide+...
            geocentric_polar_tide+...
            dry_trop+...
            wet_trop+...
            iono_corr+...
            dac; %missing sea state bias
        
        res.flags.ocean_equilibrium_tide(:)=1;
        res.flags.ocean_longperiod_tide(:)=1;
        res.flags.ocean_loading_tide(:)=1;
        res.flags.solidearth_tide(:)=1;
        res.flags.geocentric_polar_tide(:)=1;
        res.flags.dry_trop(:)=1;
        res.flags.wet_trop(:)=1;
        res.flags.dac(:)=1;
    case 'sea_ice'
        res.total_geo_corr=ocean_equilibrium_tide+...
            ocean_longperiod_tide+...
            ocean_loading_tide+...
            solidearth_tide+...
            geocentric_polar_tide+...
            dry_trop+...
            wet_trop+...
            iono_corr+...
            inv_bar;
        res.flags.ocean_equilibrium_tide(:)=1;
        res.flags.ocean_longperiod_tide(:)=1;
        res.flags.ocean_loading_tide(:)=1;
        res.flags.solidearth_tide(:)=1;
        res.flags.geocentric_polar_tide(:)=1;
        res.flags.dry_trop(:)=1;
        res.flags.wet_trop(:)=1;
        res.flags.inv_bar(:)=1;
    case 'land_ice'
        res.total_geo_corr=res.ocean_loading_tide+...
            res.solidearth_tide+...
            res.geocentric_polar_tide+...
            res.dry_trop+...
            res.wet_trop+...
            iono_corr;
        res.flags.ocean_loading_tide(:)=1;
        res.flags.solidearth_tide(:)=1;
        res.flags.geocentric_polar_tide(:)=1;
        res.flags.dry_trop(:)=1;
        res.flags.wet_trop(:)=1;
        
end



end

