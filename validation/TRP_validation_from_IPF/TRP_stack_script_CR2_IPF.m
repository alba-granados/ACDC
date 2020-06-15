mode='SAR';
c_cst =299792458;
filename_FBR='./inputs/CS_OFFL_SIR1SAR_FR_20151107T112552_20151107T112605_C001.DBL';
antenna=1;
zpf=128;
str_selec='STR_ATT_REF';
plots_flag=1;
STACK_FILE = './inputs/CS_OPER_BKPSARSSSS_20151107T112551_20151107T112627_0001.DBL';
L1B_FILE   = './inputs/CS_OFFL_SIR_SAR_1B_20151107T112552_20151107T112605_C001.DBL';
passDate= '20151107';
config_selec='A12';
TRP_delay   = 9.88;
TRP_delay   = 0;
[out_FBR,record_numFBR] = auto_readFBR(filename_FBR, true,mode);
[pos]=stack_detection_C(STACK_FILE,record_numFBR, mode);
[wvf_stack, times] = read_bkpsinssss2_isard_C(STACK_FILE, passDate,record_numFBR, pos,antenna,plots_flag,config_selec,str_selec,'./results/'); 
[out_L1B, record_num_L1B] = auto_readL1b_all_modes_C(L1B_FILE,1,0); 
nbeams = wvf_stack.t_header.CONTRIB_BEAMS;
nburst = nbeams;
nstack = pos+1;    % It is added 1 because it has been enumerated since cero
nsamples_L1 = length(out_L1B(1).waveform_group(1).averaged_power_echo_waveform); 

% ----- ZERO-PADDING -----
     
    complex_wvf_A1 = complex(wvf_stack.A1_Re,wvf_stack.A1_Img);
    nsamples_stack = length(complex_wvf_A1);
    nwvf = 0;
    start = 1;
    final = 0;
    Sy2_A1=zeros(1,nsamples_L1*zpf);
    wvf_A1=zeros(1,nsamples_L1*zpf);
    maxpos_A1=zeros(1,nbeams);
    init_sample = nsamples_stack-nsamples_L1+1;
    end_sample  = nsamples_stack;

    
    %--------------ANT-1 ---------------------------------
    for count = 1:nbeams
%         customText = 'ANT1 AZ FFT ...';
%         percentageDone =  count / nbeams;
%         stopBar= progressbar(percentageDone, 0, customText);
%         if (stopBar) 
%             break; 
%         end

        nwvf = nwvf +1;
        wvf_A1          = fft(complex_wvf_A1(count,:));         
        Sy_A1           = CalcSpectraPulse2(wvf_A1', zpf);
        Sy2_A1          = fftshift(Sy_A1);
%         Sy2_A1          = Sy_A1(1:(nsamples_stack*zpf)); % Zero-pading interpolates from 1 to nsamples_stack*zpf. I just need to keep the values from 1 to 512
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [maxval_A1(nwvf),b_A1]        = max(abs(Sy2_A1));
        maxpos_A1(nwvf) = b_A1/zpf-1/zpf+1;
%         modSy2_A1(:,nwvf) = abs(Sy2_A1(:,nwvf));
    end
    [~,central_beam_A1]=max(maxval_A1);
    nwvf = 0;
    
    
    DB=nstack;
   while (DB>20) 
      DB =DB-20;
   end   
    record = ceil(nstack/20);
    bandwidth_value = 320;
    offset_inst_corr = out_L1B(record).measurements_group(DB).instrument_range_correction_tx_rx; % [millimiters]
    win_del_stack      = out_FBR(record).measurements_group(DB).window_delay ; % [ps] 
    win_del_stack_L1   = out_L1B(record).measurements_group(DB).window_delay ; % [ps] 
    [~,max_pos_L1]= max(out_L1B(record).waveform_group(DB).averaged_power_echo_waveform);
    % remove TRP_Delay
    win_del_stack      = win_del_stack - TRP_delay/c_cst*1e12;
    win_del_stack_L1   = win_del_stack_L1 - TRP_delay/c_cst*1e12;

    temps_A1           = win_del_stack_L1+(maxpos_A1'- (nsamples_stack-nsamples_L1)/2 - (nsamples_stack/2+1))*(1000000/(2*bandwidth_value)) ;%%Resolution on Reprocessed Data Half than old processor
    temps_A1_L1        = win_del_stack_L1+(max_pos_L1'- (nsamples_L1/2+1))*(1000000/(2*bandwidth_value)) ;
    range_meas_stack_A1  = temps_A1*1e-9 * c_cst/2.  ; % [millimiters]
    range_meas_L1_A1    = temps_A1_L1*1e-9 * c_cst/2.  ; % [millimiters]
    
   range_without_slant_A1 = slantcorr(wvf_stack, range_meas_stack_A1, nbeams, start, final);  % remove the slant correction (para q se vea como parabola)
   subplot(2,2,2); hold all;plot(range_without_slant_A1'*1e-3);
   subplot(2,2,4);hold all; plot(range_without_slant_A1'-theo_range(1:length(range_without_slant_A1)));
   set(gca,'YLim',[-120 120]);
    