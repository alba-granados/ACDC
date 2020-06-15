function L2_met = L2_metrics (L2)
% L2_metrics
% This algorithm computes metrics that caracterize the isardSAT L2 outputs.
% This metrics are computed with respect to a fitting of the data,
% computed with a sliding window.
% In the future, outliers could be discarded using another sliding window.


L2_met = L2;

%%
% area = 'DomeC'; %DomeC,Spirit,Vostok
% main_path = ['C:\Users\roger.ISARDSAT\Documents\Feina\SEOM\SPICE\Data\CryoSat-2\SAR\',area,'-SAR\L2\isardSAT\'];
% files = dir([main_path,'*.nc']);

% Flags
RANGE = 1;
SIG0  = 1;
ELEV  = 1;

%%
N_surf = L2_met.ISD_num_surfaces_filtered;

% Corrected range
if RANGE
    L2_met.WINDOW_RANGE = 5; %number of samples of each side of the sliding window
                      %total_number_samples = 2* N_samples_window + 1
    
    range = zeros(1,N_surf);
    for i_surf = 1:N_surf
        start               = max(i_surf-L2_met.WINDOW_RANGE,1);
        finish              = min(i_surf+L2_met.WINDOW_RANGE,N_surf);
        range(i_surf)       = sum(L2_met.ISD_range(start:finish))/(finish-start+1);
    end
%         figure; plot(rang,'LineWidth',2); hold on; plot(range(i_file,:),'r','LineWidth',2)

    % Metrics
    L2_met.range_mean_error = mean(L2_met.ISD_range-range);
    L2_met.range_rmse       = mean(sqrt((L2_met.ISD_range-range).^2));
    L2_met.range_fitting    = range;
end


% Sigma0
if SIG0
    L2_met.WINDOW_SIG0 = 5; %number of samples of each side of the sliding window
                            %total_number_samples = 2* N_samples_window + 1
    
    sigma0 = zeros(1,N_surf);
    for i_surf = 1:N_surf
        start               = max(i_surf-L2_met.WINDOW_SIG0,1);
        finish              = min(i_surf+L2_met.WINDOW_SIG0,N_surf);
        sigma0(i_surf)      = sum(L2_met.ISD_sigma0(start:finish))/(finish-start+1);
    end
%         figure; plot(sig0,'LineWidth',2); hold on; plot(sigma0(i_file,:),'r','LineWidth',2)

    % Metrics
    L2_met.sig0_mean_error  = mean(L2_met.ISD_sigma0-sigma0);
    L2_met.sig0_rmse        = sqrt(mean((L2_met.ISD_sigma0-sigma0).^2));
    L2_met.sig0_fitting     = sigma0;
end


% Elevation
if ELEV
    L2_met.WINDOW_ELEV = 5; %number of samples of each side of the sliding window
                    %total_number_samples = 2* N_samples_window + 1
    
    elevation = zeros(1,N_surf);
    for i_surf = 1:N_surf
        start               = max(i_surf-L2_met.WINDOW_ELEV,1);
        finish              = min(i_surf+L2_met.WINDOW_ELEV,N_surf);
        elevation(i_surf)   = sum(L2_met.ISD_SSH(start:finish))/(finish-start+1);
    end
%         figure; plot(elev,'LineWidth',2); hold on; plot(elevation(i_file,:),'r','LineWidth',2)

    % Metrics
    L2_met.elev_mean_error  = mean(L2_met.ISD_SSH-elevation);
    L2_met.elev_rmse        = sqrt(mean((L2_met.ISD_SSH-elevation).^2));
    L2_met.elev_fitting     = elevation;
end




