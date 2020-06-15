ref_alt= 5060;

for i_surf=1:L1BS.N_total_surf_loc
    
            diff=abs(L1BS.alt_surf(i_surf)- ref_alt);
            if( diff < 60 )     % 5 km
                Color(i_surf,:) = [0 0.5 0]; %green      
            elseif( 60 <= diff && diff < 250 ) % between 5 and 10 km
                Color(i_surf,:) = [0 1 1];%blue
            elseif( 250 <= diff && diff < 500 ) % between 5 and 10 km
                Color(i_surf,:) = [1 0.6 0.2]; %orange
            elseif( 500 <= diff && diff < 1000 ) % between 10 and 30 km
                Color(i_surf,:) = [0.847 0.161 0]; %red
            else
                Color(i_surf,:) = [0 0 0];
            end
  
end

figure; scatter3(L1BS.lon_surf,L1BS.lat_surf,L1BS.alt_surf,40,Color,'fill','MarkerEdgeColor','k');
figlabels('Longitude [degrees]','Latitude [degrees]','Window Position [m]','',16);