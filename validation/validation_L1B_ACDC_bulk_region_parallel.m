function validation_L1B_ACDC_bulk_region_parallel(input_path_L1B_ISR,path_comparison_results,varargin)
time_init=tic;
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off','MATLAB:DELETE:FileNotFound');
version_matlab=version;
%==========================================================================
%==========================HANDLING input argument=========================
%==========================================================================
if(nargin<2 || nargin>(2+9*2))
    error('Wrong number of input parameters');
end
p = inputParser;
p.addParamValue('input_path_L2_ESA',{''},@(x)ischar(x));
p.addParamValue('figures_visible',0);
p.addParamValue('num_pools',1);
p.addParamValue('flag_outliers_removal',0);
p.addParamValue('type_outliers_removal','percentiles');
p.addParamValue('sh_name_nc','ssh');
p.addParamValue('active_comparison_tracks',1);
p.addParamValue('active_validation_tracks',1);
p.addParamValue('filename_mask_KML',{''},@(x)ischar(x));
p.parse(varargin{:});
input_path_L2_ESA=char(p.Results.input_path_L2_ESA);
figures_visible=p.Results.figures_visible;
num_pools=p.Results.num_pools;
flag_outliers_removal=p.Results.flag_outliers_removal;
type_outliers_removal=p.Results.type_outliers_removal;
sh_name_nc=p.Results.sh_name_nc;
active_comparison_tracks=p.Results.active_comparison_tracks;
active_validation_tracks=p.Results.active_validation_tracks;
filename_mask_KML=p.Results.filename_mask_KML;
clear p;

%----------------- Define linstyles for bulk comparison -------------------
linestyle_ESA='or';
linestyle_conventional='*b';
linestyle_ACDC='or';
fontsize_xlabel_tracks=8;
size_marker=12;
thick_marker=3.0;
title_name_conventional='Conventional';
title_name_ACDC='ACDC';

%-------------- Minimum Number of Surfaces within geo mask ----------------
min_num_surf_validation=50;



filesBulk.inputPath       =   input_path_L1B_ISR;
mkdir(path_comparison_results);

filesBulk.inputFiles      =   dir(filesBulk.inputPath);
filesBulk.indexaDirs      =   find(([filesBulk.inputFiles.isdir]));
filesBulk.indexFiles      =   find(not([filesBulk.inputFiles.isdir]));
filesBulk.nFiles          =   length(filesBulk.indexFiles);             % number of input files
aux=struct2cell(filesBulk.inputFiles); aux=aux(1,:); %Keep the
filesBulk.indexFilesNC=find(~cellfun(@isempty,strfind(aux,'.nc')));
filesBulk.nFilesNC=length(filesBulk.indexFilesNC);
filesBulk.NCFiles=filesBulk.inputFiles(filesBulk.indexFilesNC);

i_files_valid=0;
%% -------------- check available files -----------------------------------
if isempty(filename_mask_KML)
    within_geo_mask=1;
    geo_mask  = [];
else
    geo_mask  = kml2lla(filename_mask_KML);
end
%--------------------------------------------------------------------------
for i_file=1:filesBulk.nFilesNC
    filename_L2_ISR=char(filesBulk.inputFiles(filesBulk.indexFilesNC(i_file)).name);
    data_string=filename_L2_ISR(17:17+30);
    file_char_name=filename_L2_ISR(17+30+1+17:end-3);
    
    % ---- Checking whether the track within geomask if available ---------
    if ~isempty(filename_mask_KML)
        LAT_ISR_L2=double(ncread([input_path_L1B_ISR filename_L2_ISR],'lat_20_ku')).';        
        LON_ISR_L2=double(ncread([input_path_L1B_ISR filename_L2_ISR],'lon_20_ku')).';
        
        idx_lt_0= LON_ISR_L2<0;
        %longitudes +- values (-180,180)
        if any(idx_lt_0)
            LON_ISR_L2(idx_lt_0)=LON_ISR_L2(idx_lt_0)+360.0;
        end

        clear idx_lt_0;
                
        idx_int=inpolygon(LON_ISR_L2,LAT_ISR_L2,geo_mask.coord(:,1),geo_mask.coord(:,2));
        if ~any(idx_int)
            disp(strcat('Track,',{' '},filename_L2_ISR,{' '},'outside the limits of the geographical mask'))
            within_geo_mask=0;
        else
            if (length(find(idx_int)) < min_num_surf_validation)
                disp(strcat('Track,',{' '},filename_L2_ISR,{' '},'has a # surfaces within geographical mask below threshold'))
                within_geo_mask=0;
            else
                within_geo_mask=1;
            end            
        end
        
        clear LAT_ISR_L2 LON_ISR_L2 idx_int idx_lt_180 idx_lt_90;
        
    end
    
    if isempty(input_path_L2_ESA)
        inputL2ESAFiles=1;
    else
        inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' data_string(1:15) '*_C001.DBL']));
        if isempty(inputL2ESAFiles)
            %add one second to initial time acquisition
            init_acq_time=datestr(datenum((data_string(1:15)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
            inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,data_string(16:31)) '*_C001.DBL']));
            if isempty(inputL2ESAFiles)
                %add one second to initial time acquisition
                end_acq_time=datestr(datenum((data_string(25:31)),'yyyymmddTHHMMSS')+1/24/60/60,'yyyymmddTHHMMSS');
                inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' strcat(data_string(1:16),end_acq_time) '*_C001.DBL']));
                if isempty(inputL2ESAFiles)
                    inputL2ESAFiles   = dir(fullfile(input_path_L2_ESA,['*' strcat(init_acq_time,'_',end_acq_time) '*_C001.DBL']));
                end
            end
            
        end
    end
       
    if ~isempty(inputL2ESAFiles) && within_geo_mask
        i_files_valid=i_files_valid+1;
        filesBulk.indexFilesNC_valid(i_files_valid)=filesBulk.indexFilesNC(i_file);
        date_file_id(i_files_valid,1:(length(data_string)+length(file_char_name)+1))=strcat(strrep(data_string(1:31),'_','-'),'_',file_char_name);        
    end

end


%% ------------------ RUN the validation for each available file ----------
%--------------------------------------------------------------------------
filesBulk.nFilesNC_valid=length(filesBulk.indexFilesNC_valid);
disp(strcat('Total number of input L1B ISR for evaluation: ',num2str(filesBulk.nFilesNC)));
disp(strcat('Total number of valid files for evaluation: ',num2str(filesBulk.nFilesNC_valid)));
if active_validation_tracks
    if num_pools~=1
        %create pools
        if str2double(version_matlab(end-5:end-2))>2013
            parpool(num_pools);
        else
            matlabpool('open',num_pools);
        end
        %% ------------- Loop per file to be processed ------------------------
        parfor i_files_valid=1:filesBulk.nFilesNC_valid
            try
                run_L1B_ACDC_validation(filesBulk,i_files_valid,input_path_L1B_ISR,...
                    path_comparison_results,...
                    'input_path_L2_ESA',input_path_L2_ESA,....
                    'figures_visible',figures_visible,...
                    'flag_outliers_removal',flag_outliers_removal,...
                    'type_outliers_removal',type_outliers_removal,...
                    'sh_name_nc',sh_name_nc,...
                    'geo_mask',geo_mask)
            catch
                disp(strcat('Some error in validation for file',{''},char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name)));
                continue;
            end
        end
        %close pools
        if str2double(version_matlab(end-5:end-2))>2013
            poolobj = gcp('nocreate');
            delete(poolobj);
        else
            matlabpool('close');
        end
    else
        for i_files_valid=1:filesBulk.nFilesNC_valid
            try
                run_L1B_ACDC_validation(filesBulk,i_files_valid,input_path_L1B_ISR,...
                    path_comparison_results,...
                    'input_path_L2_ESA',input_path_L2_ESA,....
                    'figures_visible',figures_visible,...
                    'flag_outliers_removal',flag_outliers_removal,...
                    'type_outliers_removal',type_outliers_removal,...
                    'sh_name_nc',sh_name_nc,...
                    'geo_mask',geo_mask)
            catch
                disp(strcat('Some error in validation for file',{''},char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name)));
                continue;
            end
        end
    end
end
%% -------------COMPARISON OF TRACKS --------------------------------------
%--------------------------------------------------------------------------
if active_comparison_tracks
    %loading and reordering the data into a single array of structures
    %filesBulk.inputFilesEvaluation      =   dir(fullfile(path_comparison_results,'*_L2_Evaluation.mat'));
    for i_files_valid=1:filesBulk.nFilesNC_valid
        %load(strcat(path_comparison_results,char(filesBulk.inputFilesEvaluation(i_files_valid).name)))
        name_file_L1B_ISR=char(filesBulk.inputFiles(filesBulk.indexFilesNC_valid(i_files_valid)).name);
        name_file_L1B_ISR=name_file_L1B_ISR(17:17+30);
        %    aux=name_file_L2_ISR;
        %    aux2=strcat(strrep(char(name_file_L2_ISR),'.nc','_'),'L2_Evaluation.mat');
        load(strcat(path_comparison_results,name_file_L1B_ISR,'/',strcat(name_file_L1B_ISR,'_L2_Evaluation.mat')));
        SIGMA0(i_files_valid)=res.SIGMA0;
        SSH(i_files_valid)=res.SSH;
        SWH(i_files_valid)=res.SWH;
        COR(i_files_valid)=res.COR;
    end
    save(strcat(path_comparison_results,'L1B_ACDC_Bulk_validation_information.mat'),'SIGMA0','SSH','SWH','COR');
    
    
    %% ----------  Ploting ----------------------------------------------------
    xlabels=date_file_id;
    if figures_visible
        set(0, 'DefaultFigureVisible', 'on');
    else
        set(0, 'DefaultFigureVisible', 'off');
    end
    set(0,'defaultLineMarkerSize',size_marker);  % set the default line marker size
    set(0,'defaultLineLineWidth',thick_marker);
    %% ---------------------------  SSH ---------------------------------------
    if ~isempty(input_path_L2_ESA)
        %--------------------------------------------------------------------------
        %$$$$$$$$$$$$$$$$$ Comparison retracker w.r.t ESA $$$$$$$$$$$$$$$$$$$$$$$$$
        %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        figure;
        legend_text={''};        
        results=[SSH(:).CONVENTIONAL];
        plot([results.RMSE_error_L2],linestyle_conventional)
        legend_text=[legend_text,{strcat('ESA-ISR (',title_name_conventional,')')}];
        hold on;
        grid on;
        results=[SSH(:).ACDC];
        plot([results.RMSE_error_L2],linestyle_ACDC)
        legend_text=[legend_text,{strcat('ESA-ISR (',title_name_ACDC,')')}];                      
        title(strcat('RMSE error on',{' '},upper(sh_name_nc),': ESA & isardSAT'))
        legend(legend_text(~cellfun(@isempty,legend_text)));
        xlabel_str='Track'; ylabel(strcat('RMSE_{',upper(sh_name_nc),{'} '},'[m]'),'Interpreter','Tex');
        xlabel(xlabel_str)        
        print('-dpng',strcat(path_comparison_results,'RMSE_SSH_comp_ESA_ISR.png'))
        %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        figure;
        legend_text={''};        
        results=[SSH(:).CONVENTIONAL];
        plot([results.mean_error_L2],linestyle_conventional)
        legend_text=[legend_text,{strcat('ESA-ISR (',title_name_conventional,')')}];
        hold on;
        grid on;        
        results=[SSH(:).ACDC];
        plot([results.mean_error_L2],linestyle_ACDC)
        legend_text=[legend_text,{strcat('ESA-ISR (',title_name_ACDC,')')}];        
        title(strcat('Mean error on',{' '},upper(sh_name_nc),': ESA & isardSAT'))
        legend(legend_text(~cellfun(@isempty,legend_text)));
        xlabel_str='Track'; ylabel(strcat('\epsilon_{',upper(sh_name_nc),'} [m]'),'Interpreter','Tex');
        xlabel(xlabel_str)
        print('-dpng',strcat(path_comparison_results,'Mean_error_SSH_comp_ESA_ISR.png'))
    end

    %--------------------------------------------------------------------------
    %$$$$$$$$$$$$$$$$$ Comparison retrackers ISD $$$$$$$$$$$$$$$$$$$$$$$$$
    %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[SSH(:).CONVENTIONAL];
    plot([results.RMSE_error_ACDC],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR conventional-ACDC')}];
    hold on;    
    grid on;
    title(strcat('RMSE error on',{' '},upper(sh_name_nc),': isardSAT conventional-ACDC'))
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel(strcat('RMSE_{',upper(sh_name_nc),{'} '},'[m]'),'Interpreter','Tex');
    xlabel(xlabel_str)    
    print('-dpng',strcat(path_comparison_results,'RMSE_SSH_comp_ISR.png'))
    %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[SSH(:).CONVENTIONAL];
    plot([results.mean_error_ACDC],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR conventional-ACDC')}];
    hold on;
    grid on;    
    title(strcat('Mean error on',{' '},upper(sh_name_nc),': isardSAT conventional-ACDC'))
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel(strcat('\epsilon_{',upper(sh_name_nc),'} [m]'),'Interpreter','Tex');
    xlabel(xlabel_str)
    print('-dpng',strcat(path_comparison_results,'Mean_error_SSH_comp_ISR.png'))
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&    
    figure;
    legend_text={''};    
    results=[SSH(:).CONVENTIONAL];
    plot([results.rmse_fitting],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR (',title_name_conventional,')')}];
    hold on;
    grid on;
    
    results=[SSH(:).ACDC];
    plot([results.rmse_fitting],linestyle_ACDC)
    legend_text=[legend_text,{strcat('ISR (',title_name_ACDC,')')}];
    
    if ~isempty(input_path_L2_ESA)
        results=[SSH(:).ESA_L2];
        plot([results.rmse_fitting],linestyle_ESA)
        legend_text=[legend_text,{'ESA'}];
        title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': ESA & isardSAT (ISR)'))
    else
        title(strcat('RMSE error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR)'))
    end
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('RMSE_{SSH} [m]','Interpreter','Tex');
    xlabel(xlabel_str)
    % %-------------- set the x-tick locations ----------------------------------
    % % Reduce the size of the axis so that all the labels fit in the figure.
    % pos = get(gca,'Position');
    % set(gca,'Position',[pos(1), .2, pos(3) .65])
    % Xt=1:1:filesBulk.nFilesNC_valid;
    % Xl=[1 filesBulk.nFilesNC_valid];
    % Yl=get(gca,'ylim');
    % set(gca,'XTick',Xt,'XLim',Xl);
    % ax = axis;    % Current axis limits
    % axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
    % % Place the text labels
    % t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
    % set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    %       'Rotation',45);
    % % Remove the default labels
    % set(gca,'XTickLabel','')
    % % Get the Extent of each text object.  This
    % % loop is unavoidable.
    % for i = 1:length(t)
    %   ext(i,:) = get(t(i),'Extent');
    % end
    % % Determine the lowest point.  The X-label will be
    % % placed so that the top is aligned with this point.
    % LowYPoint = min(ext(:,2));
    % % Place the axis label at this point
    % XMidPoint = Xl(1)+abs(diff(Xl))/2;
    % tl = text(XMidPoint,LowYPoint,xlabel_str, ...
    %           'VerticalAlignment','top', ...
    %           'HorizontalAlignment','center');
    print('-dpng',strcat(path_comparison_results,'RMSE_fitted_SSH.png'))
    %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    

    results=[SSH(:).CONVENTIONAL];
    plot([results.mean_error_fitting],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR (',title_name_conventional,')')}];
    hold on; 
    grid on;
    
    results=[SSH(:).ACDC];
    plot([results.mean_error_fitting],linestyle_ACDC)
    legend_text=[legend_text,{strcat('ISR (',title_name_ACDC,')')}];
    
    if ~isempty(input_path_L2_ESA)
        results=[SSH(:).ESA_L2];
        plot([results.mean_error_fitting],linestyle_ESA)
        legend_text=[legend_text,{'ESA'}];
        title(strcat('Mean error on fitted',{' '},upper(sh_name_nc),': ESA & isardSAT (ISR)'))
    else
        title(strcat('Mean error on fitted',{' '},upper(sh_name_nc),': isardSAT (ISR)'))
    end
    
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('\epsilon_{SSH} [m]','Interpreter','Tex');
    xlabel(xlabel_str)
    % %-------------- set the x-tick locations ----------------------------------
    % % Reduce the size of the axis so that all the labels fit in the figure.
    % pos = get(gca,'Position');
    % set(gca,'Position',[pos(1), .2, pos(3) .65])
    % Xt=1:1:filesBulk.nFilesNC_valid;
    % Xl=[1 filesBulk.nFilesNC_valid];
    % Yl=get(gca,'ylim');
    % set(gca,'XTick',Xt,'XLim',Xl);
    % ax = axis;    % Current axis limits
    % axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
    % % Place the text labels
    % t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
    % set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    %       'Rotation',45);
    % % Remove the default labels
    % set(gca,'XTickLabel','')
    % % Get the Extent of each text object.  This
    % % loop is unavoidable.
    % for i = 1:length(t)
    %   ext(i,:) = get(t(i),'Extent');
    % end
    % % Determine the lowest point.  The X-label will be
    % % placed so that the top is aligned with this point.
    % LowYPoint = min(ext(:,2));
    % % Place the axis label at this point
    % XMidPoint = Xl(1)+abs(diff(Xl))/2;
    % tl = text(XMidPoint,LowYPoint,xlabel_str, ...
    %           'VerticalAlignment','top', ...
    %           'HorizontalAlignment','center');
    print('-dpng',strcat(path_comparison_results,'Mean_error_fitted_SSH.png'))
    %% ---------------------------  SIGMA0 ------------------------------------
    %$$$$$$$$$$$$$$$$$ Comparison retracker w.r.t ESA $$$$$$$$$$$$$$$$$$$$$$$$$
    %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    if ~isempty(input_path_L2_ESA)
        %--------------------------------------------------------------------------
        %$$$$$$$$$$$$$$$$$ Comparison retracker w.r.t ESA $$$$$$$$$$$$$$$$$$$$$$$$$
        %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        figure;
        legend_text={''};
        results=[SIGMA0(:).CONVENTIONAL];
        plot([results.RMSE_error_L2],linestyle_conventional)
        legend_text=[legend_text,{strcat('ESA-ISR (',title_name_conventional,')')}];
        hold on;
        grid on;
        results=[SIGMA0(:).ACDC];
        plot([results.RMSE_error_L2],linestyle_ACDC)
        legend_text=[legend_text,{strcat('ESA-ISR (',title_name_ACDC,')')}];
        title(strcat('RMSE error on',{' '},upper(sh_name_nc),': ESA & isardSAT'))
        legend(legend_text(~cellfun(@isempty,legend_text)));
        xlabel_str='Track'; ylabel('RMSE_{\sigma^0} [dB]','Interpreter','Tex');
        xlabel(xlabel_str)
        print('-dpng',strcat(path_comparison_results,'RMSE_SIG0_comp_ESA_ISR.png'))
        %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        figure;
        legend_text={''};
        results=[SIGMA0(:).CONVENTIONAL];
        plot([results.mean_error_L2],linestyle_conventional)
        legend_text=[legend_text,{strcat('ESA-ISR (',title_name_conventional,')')}];
        hold on;
        grid on;
        results=[SIGMA0(:).ACDC];
        plot([results.mean_error_L2],linestyle_ACDC)
        legend_text=[legend_text,{strcat('ESA-ISR (',title_name_ACDC,')')}];
        title('Mean error on \sigma^0: ESA & isardSAT','Interpreter','Tex')
        legend(legend_text(~cellfun(@isempty,legend_text)));
        ylabel('\epsilon_{\sigma^0} [dB]','Interpreter','Tex');
        xlabel(xlabel_str)
        print('-dpng',strcat(path_comparison_results,'Mean_error_SIG0_comp_ESA_ISR.png'))
    end
%--------------------------------------------------------------------------
    %$$$$$$$$$$$$$$$$$ Comparison retrackers ISD $$$$$$$$$$$$$$$$$$$$$$$$$
    %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[SIGMA0(:).CONVENTIONAL];
    plot([results.RMSE_error_ACDC],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR conventional-ACDC')}];
    hold on;  
    grid on;
    title(strcat('RMSE error on \sigma^0: isardSAT conventional-ACDC'))
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('RMSE_{\sigma^0} [dB]','Interpreter','Tex');
    xlabel(xlabel_str)    
    print('-dpng',strcat(path_comparison_results,'RMSE_SIG0_comp_ISR.png'))
    %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[SIGMA0(:).CONVENTIONAL];
    plot([results.mean_error_ACDC],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR conventional-ACDC')}];
    hold on;
    grid on;    
    title(strcat('Mean error on \sigma^0: isardSAT conventional-ACDC'))
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('\epsilon_{\sigma^0} [dB]','Interpreter','Tex');
    xlabel(xlabel_str)
    print('-dpng',strcat(path_comparison_results,'Mean_error_SIG0_comp_ISR.png'))    
        
    %$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[SIGMA0(:).CONVENTIONAL];
    plot([results.rmse_fitting],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR (',title_name_conventional,')')}];
    hold on;
    grid on;
    
    results=[SIGMA0(:).ACDC];
    plot([results.rmse_fitting],linestyle_ACDC)
    legend_text=[legend_text,{strcat('ISR (',title_name_ACDC,')')}];
    
    if ~isempty(input_path_L2_ESA)
        results=[SIGMA0(:).ESA_L2];
        plot([results.rmse_fitting],linestyle_ESA)
        legend_text=[legend_text,{'ESA'}];
        title(strcat('RMSE error on fitted \sigma^0: ESA & isardSAT (ISR)'))
    else
        title(strcat('RMSE error on fitted \sigma^0: isardSAT (ISR)'))
    end
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('RMSE_{\sigma^0} [dB]','Interpreter','Tex');
    xlabel(xlabel_str)
    % %-------------- set the x-tick locations ----------------------------------
    % % Reduce the size of the axis so that all the labels fit in the figure.
    % pos = get(gca,'Position');
    % set(gca,'Position',[pos(1), .2, pos(3) .65])
    % Xt=1:1:filesBulk.nFilesNC_valid;
    % Xl=[1 filesBulk.nFilesNC_valid];
    % Yl=get(gca,'ylim');
    % set(gca,'XTick',Xt,'XLim',Xl);
    % ax = axis;    % Current axis limits
    % axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
    % % Place the text labels
    % t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
    % set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    %       'Rotation',45);
    % % Remove the default labels
    % set(gca,'XTickLabel','')
    % % Get the Extent of each text object.  This
    % % loop is unavoidable.
    % for i = 1:length(t)
    %   ext(i,:) = get(t(i),'Extent');
    % end
    % % Determine the lowest point.  The X-label will be
    % % placed so that the top is aligned with this point.
    % LowYPoint = min(ext(:,2));
    % % Place the axis label at this point
    % XMidPoint = Xl(1)+abs(diff(Xl))/2;
    % tl = text(XMidPoint,LowYPoint,xlabel_str, ...
    %           'VerticalAlignment','top', ...
    %           'HorizontalAlignment','center');
    print('-dpng',strcat(path_comparison_results,'RMSE_fitted_SIGMA0.png'))
    %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    

    results=[SIGMA0(:).CONVENTIONAL];
    plot([results.mean_error_fitting],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR (',title_name_conventional,')')}];
    hold on; 
    grid on;
    
    results=[SIGMA0(:).ACDC];
    plot([results.mean_error_fitting],linestyle_ACDC)
    legend_text=[legend_text,{strcat('ISR (',title_name_ACDC,')')}];
    
    if ~isempty(input_path_L2_ESA)
        results=[SIGMA0(:).ESA_L2];
        plot([results.mean_error_fitting],linestyle_ESA)
        legend_text=[legend_text,{'ESA'}];
        title(strcat('Mean error on fitted \sigma^0: ESA & isardSAT (ISR)'))
    else
        title(strcat('Mean error on fitted \sigma^0: isardSAT (ISR)'))
    end
    
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('\epsilon_{\sigma^0} [dB]','Interpreter','Tex');
    xlabel(xlabel_str)
    % %-------------- set the x-tick locations ----------------------------------
    % % Reduce the size of the axis so that all the labels fit in the figure.
    % pos = get(gca,'Position');
    % set(gca,'Position',[pos(1), .2, pos(3) .65])
    % Xt=1:1:filesBulk.nFilesNC_valid;
    % Xl=[1 filesBulk.nFilesNC_valid];
    % Yl=get(gca,'ylim');
    % set(gca,'XTick',Xt,'XLim',Xl);
    % ax = axis;    % Current axis limits
    % axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
    % % Place the text labels
    % t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
    % set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    %       'Rotation',45);
    % % Remove the default labels
    % set(gca,'XTickLabel','')
    % % Get the Extent of each text object.  This
    % % loop is unavoidable.
    % for i = 1:length(t)
    %   ext(i,:) = get(t(i),'Extent');
    % end
    % % Determine the lowest point.  The X-label will be
    % % placed so that the top is aligned with this point.
    % LowYPoint = min(ext(:,2));
    % % Place the axis label at this point
    % XMidPoint = Xl(1)+abs(diff(Xl))/2;
    % tl = text(XMidPoint,LowYPoint,xlabel_str, ...
    %           'VerticalAlignment','top', ...
    %           'HorizontalAlignment','center');
    print('-dpng',strcat(path_comparison_results,'Mean_error_fitted_SIGMA0.png'))
    
    %% ---------------------------  SWH ---------------------------------------
    %$$$$$$$$$$$$$$$$$ Comparison retrackers ISD $$$$$$$$$$$$$$$$$$$$$$$$$
    %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[SWH(:).CONVENTIONAL];
    plot([results.RMSE_error_ACDC],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR conventional-ACDC')}];
    hold on;   
    grid on;
    title(strcat('RMSE error on SWH: isardSAT conventional-ACDC'))
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('RMSE_{SWH} [m]','Interpreter','Tex');
    xlabel(xlabel_str)    
    print('-dpng',strcat(path_comparison_results,'RMSE_SWH_comp_ISR.png'))
    %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[SWH(:).CONVENTIONAL];
    plot([results.mean_error_ACDC],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR conventional-ACDC')}];
    hold on;
    grid on;    
    title(strcat('Mean error on SWH: isardSAT conventional-ACDC'))
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('\epsilon_{SWH} [m]','Interpreter','Tex');
    xlabel(xlabel_str)
    print('-dpng',strcat(path_comparison_results,'Mean_error_SWH_comp_ISR.png'))    
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$ Fitting $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[SWH(:).CONVENTIONAL];
    plot([results.rmse_fitting],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR (',title_name_conventional,')')}];
    hold on;
    grid on;
    
    results=[SWH(:).ACDC];
    plot([results.rmse_fitting],linestyle_ACDC)
    legend_text=[legend_text,{strcat('ISR (',title_name_ACDC,')')}];
    
    title(strcat('RMSE error on fitted SWH: isardSAT (ISR)'))
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('RMSE_{SWH} [m]','Interpreter','Tex');
    xlabel(xlabel_str)
    % %-------------- set the x-tick locations ----------------------------------
    % % Reduce the size of the axis so that all the labels fit in the figure.
    % pos = get(gca,'Position');
    % set(gca,'Position',[pos(1), .2, pos(3) .65])
    % Xt=1:1:filesBulk.nFilesNC_valid;
    % Xl=[1 filesBulk.nFilesNC_valid];
    % Yl=get(gca,'ylim');
    % set(gca,'XTick',Xt,'XLim',Xl);
    % ax = axis;    % Current axis limits
    % axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
    % % Place the text labels
    % t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
    % set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    %       'Rotation',45);
    % % Remove the default labels
    % set(gca,'XTickLabel','')
    % % Get the Extent of each text object.  This
    % % loop is unavoidable.
    % for i = 1:length(t)
    %   ext(i,:) = get(t(i),'Extent');
    % end
    % % Determine the lowest point.  The X-label will be
    % % placed so that the top is aligned with this point.
    % LowYPoint = min(ext(:,2));
    % % Place the axis label at this point
    % XMidPoint = Xl(1)+abs(diff(Xl))/2;
    % tl = text(XMidPoint,LowYPoint,xlabel_str, ...
    %           'VerticalAlignment','top', ...
    %           'HorizontalAlignment','center');
    print('-dpng',strcat(path_comparison_results,'RMSE_fitted_SWH.png'))
    %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    

    results=[SWH(:).CONVENTIONAL];
    plot([results.mean_error_fitting],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR (',title_name_conventional,')')}];
    hold on; 
    grid on;
    
    results=[SWH(:).ACDC];
    plot([results.mean_error_fitting],linestyle_ACDC)
    legend_text=[legend_text,{strcat('ISR (',title_name_ACDC,')')}];
    
    title(strcat('Mean error on fitted SWH: isardSAT (ISR)'))
    
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('\epsilon_{SWH} [m]','Interpreter','Tex');
    xlabel(xlabel_str)
    % %-------------- set the x-tick locations ----------------------------------
    % % Reduce the size of the axis so that all the labels fit in the figure.
    % pos = get(gca,'Position');
    % set(gca,'Position',[pos(1), .2, pos(3) .65])
    % Xt=1:1:filesBulk.nFilesNC_valid;
    % Xl=[1 filesBulk.nFilesNC_valid];
    % Yl=get(gca,'ylim');
    % set(gca,'XTick',Xt,'XLim',Xl);
    % ax = axis;    % Current axis limits
    % axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
    % % Place the text labels
    % t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
    % set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    %       'Rotation',45);
    % % Remove the default labels
    % set(gca,'XTickLabel','')
    % % Get the Extent of each text object.  This
    % % loop is unavoidable.
    % for i = 1:length(t)
    %   ext(i,:) = get(t(i),'Extent');
    % end
    % % Determine the lowest point.  The X-label will be
    % % placed so that the top is aligned with this point.
    % LowYPoint = min(ext(:,2));
    % % Place the axis label at this point
    % XMidPoint = Xl(1)+abs(diff(Xl))/2;
    % tl = text(XMidPoint,LowYPoint,xlabel_str, ...
    %           'VerticalAlignment','top', ...
    %           'HorizontalAlignment','center');
    print('-dpng',strcat(path_comparison_results,'Mean_error_fitted_SWH.png'))
    
    %% ----------------------  PEARSON CORR COEF. -----------------------------
        %$$$$$$$$$$$$$$$$$ Comparison retrackers ISD $$$$$$$$$$$$$$$$$$$$$$$$$
    %&&&&&&&&&&&&&&&&&&&&& RMSE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[COR(:).CONVENTIONAL];
    plot([results.RMSE_error_ACDC],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR conventional-ACDC')}];
    hold on;        
    title(strcat('RMSE error on corr. coeff. \rho_{pearson}: isardSAT conventional-ACDC'))
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('RMSE_{\rho_{pearson}} [%]','Interpreter','Tex');
    xlabel(xlabel_str)    
    print('-dpng',strcat(path_comparison_results,'RMSE_COR_comp_ISR.png'))
    %&&&&&&&&&&&&&&&&&&&&& MEAN ERROR &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    figure;
    legend_text={''};    
    results=[COR(:).CONVENTIONAL];
    plot([results.mean_error_ACDC],linestyle_conventional)
    legend_text=[legend_text,{strcat('ISR conventional-ACDC')}];
    hold on;
    grid on;    
    title(strcat('Mean error on corr. coeff. \rho_{pearson}: isardSAT conventional-ACDC'))
    legend(legend_text(~cellfun(@isempty,legend_text)));
    xlabel_str='Track'; ylabel('\epsilon_{\rho_{pearson}} [%]','Interpreter','Tex');
    xlabel(xlabel_str)
    print('-dpng',strcat(path_comparison_results,'Mean_error_COR_comp_ISR.png'))    
    
    %$$$$$$$$$$$$$$$$$$$$$$$$$ MEAN VALUE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    figure;
    legend_text={''};
    results=[COR(:).CONVENTIONAL];
    plot([results.mean],linestyle_conventional)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_conventional,')')}];
    hold on;
    grid on;
    results=[COR(:).ACDC];
    plot([results.mean],linestyle_ACDC)
    legend_text=[legend_text,{strcat('L2 ISR (',title_name_ACDC,')')}];
    hold on;
    grid on;
    title('Mean value pearson corr. coeff. \rho_{pearson}: isardSAT','Interpreter','Tex')
    legend(legend_text(~cellfun(@isempty,legend_text)));
    ylabel('\rho_{pearson} [%]','Interpreter','Tex');
    xlabel(xlabel_str)
    % %-------------- set the x-tick locations ----------------------------------
    % % Reduce the size of the axis so that all the labels fit in the figure.
    % pos = get(gca,'Position');
    % set(gca,'Position',[pos(1), .2, pos(3) .65])
    % Xt=1:1:filesBulk.nFilesNC_valid;
    % Xl=[1 filesBulk.nFilesNC_valid];
    % Yl=get(gca,'ylim');
    % set(gca,'XTick',Xt,'XLim',Xl);
    % ax = axis;    % Current axis limits
    % axis(axis);    % Set the axis limit modes (e.g. XLimMode) to manual
    % % Place the text labels
    % t = text(Xt,Yl(1)*ones(1,length(Xt)),xlabels(1:1:filesBulk.nFilesNC_valid,:),'FontSize',fontsize_xlabel_tracks);
    % set(t,'HorizontalAlignment','right','VerticalAlignment','top', ...
    %       'Rotation',45);
    % % Remove the default labels
    % set(gca,'XTickLabel','')
    % % Get the Extent of each text object.  This
    % % loop is unavoidable.
    % for i = 1:length(t)
    %   ext(i,:) = get(t(i),'Extent');
    % end
    % % Determine the lowest point.  The X-label will be
    % % placed so that the top is aligned with this point.
    % LowYPoint = min(ext(:,2));
    % % Place the axis label at this point
    % XMidPoint = Xl(1)+abs(diff(Xl))/2;
    % tl = text(XMidPoint,LowYPoint,xlabel_str, ...
    %           'VerticalAlignment','top', ...
    %           'HorizontalAlignment','center');
    print('-dpng',strcat(path_comparison_results,'Mean_COR.png'))
    close all
end
time_end=toc(time_init);
minutes_processing = floor(time_end/60);
secs_processing = time_end - minutes_processing*60;
disp(['Validation/processing time for ',input_path_L1B_ISR,': ',num2str(minutes_processing),' minutes and ',num2str(secs_processing),' seconds']);    
end

