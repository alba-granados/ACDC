% read_all_netcdf
% this script uses readanyNETCDF_V1.m function to read all the files in the
% specified directory
clear

PLOT = 1;
SAVE = 0;

ALL = 0;
if ALL == 0
    start_cycle = 200;
    stop_cycle = 213;
end

path =  'C:\Users\roger.ISARDSAT\Documents\Feina\Jason-2\Data\Inputs\L2\GDR\';
% path =  'C:\Users\roger.ISARDSAT\Documents\Feina\Jason-2\Data\Inputs\L2\IGDR\';


list_aux = dir(path);
if ALL == 0
    for i_list = 1:length(list_aux)
        if strcmp(list_aux(i_list).name,num2str(start_cycle)) == 1
            START = i_list;
        elseif strcmp(list_aux(i_list).name,num2str(stop_cycle)) == 1
            STOP = i_list;
            break;
        end
    end
    if start_cycle == stop_cycle
        STOP = START;
    end
else
    START = 3;
    STOP = length(list_aux);
end
list_cycle = list_aux(START:STOP);

track_num = [15,91,117,193];
% track_num = [15,91,117,193,26,154,167,230];
NTRACKS = length(track_num);
lat_lim1 = 18; %18N
lat_lim2 = 30; %30N
% lon_lim1 = 83+180; %83W
% lon_lim2 = 97+180; %97W


for i_cycle = 1:length(list_cycle)
    pathname = [path,list_cycle(i_cycle).name];
    list_file = dir([pathname,'\*.nc']);
    
    cycle_name = ['cycle',list_cycle(i_cycle).name];
    
    data.(cycle_name).ssha = [];
    data.(cycle_name).mss = [];
    data.(cycle_name).alt = [];
    data.(cycle_name).lat = [];
    data.(cycle_name).lon = [];

    for i_file = 1:NTRACKS
%         progressbar(i_file/NTRACKS)
        
        filename = [pathname,'\',list_file(i_file).name];
        file_data = readanyNETCDF_V2(filename);
        
        %putting the data in a vector
        ssha    = double(file_data.data.ssha)*1e-3;
        mss     = double(file_data.data.mean_sea_surface)*1e-4;
        alt     = double(file_data.data.alt)*1e-4 + 1.3e6;
        lat     = double(file_data.data.lat)*1e-6;
        lon     = double(file_data.data.lon)*1e-6;
%         flag_acq_mode   = file_data.data.alt_state_flag_acq_mode_20hz;
        
%         num_rec = length(lat);
        
        %taking only the information for the desired area
        indexes = lat>=lat_lim1 & lat<=lat_lim2;
        num_rec_filt = sum(indexes);
        
        %storing the data
        data.(cycle_name).num_rec_filt(i_file) = num_rec_filt;
%         data.(cycle_name).flag_acq_mode(i_file,1:num_rec_filt) = flag_acq_mode(indexes);
        
        data.(cycle_name).ssha_mtx(i_file,1:num_rec_filt) = ssha(indexes);
        data.(cycle_name).mss_mtx(i_file,1:num_rec_filt)  = mss(indexes);
        data.(cycle_name).elev_mtx(i_file,1:num_rec_filt) = data.(cycle_name).ssha_mtx(i_file,1:num_rec_filt) + data.(cycle_name).mss_mtx(i_file,1:num_rec_filt);
        data.(cycle_name).alt_mtx(i_file,1:num_rec_filt)  = alt(indexes);
        data.(cycle_name).lat_mtx(i_file,1:num_rec_filt)  = lat(indexes);
        data.(cycle_name).lon_mtx(i_file,1:num_rec_filt)  = lon(indexes);
        
        data.(cycle_name).ssha  = [data.(cycle_name).ssha, double(data.(cycle_name).ssha_mtx(i_file,1:num_rec_filt))];
        data.(cycle_name).mss   = [data.(cycle_name).mss, double(data.(cycle_name).mss_mtx(i_file,1:num_rec_filt))];
        data.(cycle_name).elev  = data.(cycle_name).ssha + data.(cycle_name).mss;
        data.(cycle_name).alt   = [data.(cycle_name).alt, double(data.(cycle_name).alt_mtx(i_file,1:num_rec_filt))];
        data.(cycle_name).lat   = [data.(cycle_name).lat, double(data.(cycle_name).lat_mtx(i_file,1:num_rec_filt))];
        data.(cycle_name).lon   = [data.(cycle_name).lon, double(data.(cycle_name).lon_mtx(i_file,1:num_rec_filt))];
    end
    progressbar(i_cycle/length(list_cycle))
end

data.track_num = track_num;



%% ---------- PLOTTING -----------
if PLOT == 1
    NTRACKS = length(data.track_num);

    %Legend
    mylegend = [];
    for i_cycle = 1:length(list_cycle)
        mylegend = [mylegend;list_cycle(i_cycle).name];
    end

    for i_track = 1:NTRACKS
        h = figure; hold all
        for i_cycle = 1:length(list_cycle)
            cycle_name = ['cycle',list_cycle(i_cycle).name];        
            x_pre = data.(cycle_name).lat_mtx(i_track,1:data.(cycle_name).num_rec_filt(i_track));
            y_pre = data.(cycle_name).ssha_mtx(i_track,1:data.(cycle_name).num_rec_filt(i_track)) + ...
                        data.(cycle_name).mss_mtx(i_track,1:data.(cycle_name).num_rec_filt(i_track));
%             plot(x_pre,y_pre,'.','MarkerSize',10)
            
            y_pre_mean = mean(y_pre);
            ind = y_pre < y_pre_mean;
            x = x_pre(ind);
            y = y_pre(ind);
            plot(x,y,'LineWidth',2)
%             plot(x,y,'.','MarkerSize',10)
            x_mean(i_track,i_cycle) = mean(x);
            y_mean(i_track,i_cycle) = mean(y);
        end
        title(['Track #',num2str(data.track_num(i_track))])
        legend(mylegend)
        if SAVE == 1
            saveas(h,['track_num',num2str(data.track_num(i_track)),'.fig'])
        end
    end
end


if SAVE == 1
    save('J2_elevations','data')
end


clear track_num alt cycle_name cycle_num file_data filename...
    flag_acq_mode i_cycle i_file i_track indexes lat lat_lim1...
    lat_lim2 list lon num_rec num_rec_filt path pathname ssha mss...
    start_cycle stop stop_cycle h ind list_aux list_cycle list_file...
    mylegend x_pre y_pre y_pre_mean NTRACKS START STOP...
    i_list ALL PLOT SAVE x x_mean y y_mean