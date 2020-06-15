% Delete CALIBRATION and DEGRADED bursts


function [data] 			= delete_bursts(data,indexes)
    
    names = fieldnames(data);
    N_parameters = length(names);
    
    for i_parameter = 1:N_parameters
        
        if(length(size(data.(char(names(i_parameter)))))==3)
            data.(char(names(i_parameter)))(indexes,:,:)= [];  
        elseif(length(size(data.(char(names(i_parameter)))))==2)
            if(size(data.(char(names(i_parameter))),2)==data.N_total_bursts_sar_ku)
                data.(char(names(i_parameter)))(indexes)= [];
            end
        
        end
     
    end
    data.N_total_bursts_sar_ku = length(data.lat_sar_sat);
end