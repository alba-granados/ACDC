% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% CryoSat 2 calibration over transponders
%
% ---------------------------------------------------------
% Stack Detection: reads the stack data files and finds TRP stack
% 
% Calling
%   [out_stack, dates] = stack_detection(filename, record_num)
%
% Input
%   filename    : stack file to be read
%   pass        : number of the dataset to be analysed
%
% Outputs
%   [pos_of_max]   : stack that corresponds to the Transponder
%   
%
% ----------------------------------------------------------
% 
% Author:   Albert García / isardSAT
% Reviewer: Mònica Roca / isardSAT
% Last revision: Mònica Roca / isardSAT --/--/13)
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [pos_of_max] = stack_detection_C(filename, record_num, mode)
  val_of_max =0;
  iStack=0;
  fid = fopen(filename,'r', 'b');

    for i = 1:record_num*20-1 %we will read from start to position of the maximum
       
%   Progress bar
        customText = 'Looking for the TRP Stack...';
        percentageDone =  (i+1) / (record_num*20);
        
%         stopBar= progressbar(percentageDone, 0, customText);
%                 if (stopBar) 
%                     break; 
%                 end
%         
        
	    t_header = struct('day', fread(fid,1,'int32'),...
			  'seconds',fread(fid,1,'float64'),...
			  'SAMPLE',fread(fid,1,'int32'),...
			  'RW_SAMPLES', fread(fid,1,'int32'),...
			  'EXTENDED_RW_SAMPLES', fread(fid,1,'int32'),...
              'PULSES_PER_BURST', fread(fid,1,'int32'), ...  
              'ANTENNAE', fread(fid,1,'int32'),...            
              'MODE', fread(fid,1,'int32'),...                
              'CONTRIB_BEAMS', fread(fid,1,'int32'),...       
              'SUB_STACK_BEAMS', fread(fid,1,'int32'),...
              'CONTRIB_BEAMS_BEFORE_WEIGHTING', fread(fid,1,'int32') );
          
        if(isempty(t_header.day))
            
            break;  
        elseif(t_header.day < 0)
            break; 
        elseif(abs(t_header.EXTENDED_RW_SAMPLES) > 20000)
            break; 
        end
        
        hour_st = fix(t_header.seconds/3600.);
        minute_st = t_header.seconds-hour_st*3600;
        min_st = fix(minute_st / 60.);
        sec_st = minute_st(1)-min_st*60;
        dates(i+1) = t_header.seconds;

          fread(fid, 1, 'uint32');
%         fread(fid, (13+t_header.CONTRIB_BEAMS+t_header.RW_SAMPLES+t_header.SUB_STACK_BEAMS+t_header.SUB_STACK_BEAMS+100+t_header.SUB_STACK_BEAMS*t_header.RW_SAMPLES), 'float32');
          fseek(fid,(13+t_header.CONTRIB_BEAMS+t_header.RW_SAMPLES+t_header.SUB_STACK_BEAMS+t_header.SUB_STACK_BEAMS+100+t_header.SUB_STACK_BEAMS*t_header.RW_SAMPLES)*4,'cof');
          
              A1_Re = zeros(t_header.CONTRIB_BEAMS,t_header.EXTENDED_RW_SAMPLES);
              A1_Img = zeros(t_header.CONTRIB_BEAMS,t_header.EXTENDED_RW_SAMPLES);
              A2_Re = zeros(t_header.CONTRIB_BEAMS,t_header.EXTENDED_RW_SAMPLES);
              A2_Img = zeros(t_header.CONTRIB_BEAMS,t_header.EXTENDED_RW_SAMPLES);
              complex_wvf= zeros(t_header.CONTRIB_BEAMS,t_header.EXTENDED_RW_SAMPLES);
              wvf= zeros(t_header.CONTRIB_BEAMS,t_header.EXTENDED_RW_SAMPLES);
              wvf_t=zeros(t_header.CONTRIB_BEAMS,t_header.EXTENDED_RW_SAMPLES);
              
          
          for k = 1:t_header.CONTRIB_BEAMS
%             fread(fid,3,'uint32');
%             fread(fid,2,'int16');
%             fread(fid,2,'float32');
%             fread(fid,1,'int32');
%             fread(fid,5,'float32');
            fseek(fid, 3*4+2*2+3*4+1*4+5*4,'cof'); % Field 32 Jumps until waveforms
            if(strcmp(mode,'SIN'))
                t_sample_store = struct('ANTENNA1',fread(fid, [2,t_header.EXTENDED_RW_SAMPLES],'float32'),...
                                        'ANTENNA2',fread(fid, [2,t_header.EXTENDED_RW_SAMPLES],'float32' ));
    %             for j = 1:t_header.EXTENDED_RW_SAMPLES   
    %                     avg_wvf_real (k,j) = (t_sample_store.ANTENNA1(1,j) + t_sample_store.ANTENNA2(1,j))/2.;
    %                     avg_wvf_img (k,j) = (t_sample_store.ANTENNA1(2,j) + t_sample_store.ANTENNA2(2,j))/2.;
    %                     avg_wvf (k,j) = module([avg_wvf_real(k,j) ,  avg_wvf_img(k,j)]);
    %                     waveform(k,j) = module(t_sample_store.ANTENNA1(:,j));
    %                     waveform2(k,j) = module(t_sample_store.ANTENNA2(:,j));                  
    %             end
                A1_Re(k,:) = t_sample_store.ANTENNA1(1,:);
                A1_Img(k,:) = t_sample_store.ANTENNA1(2,:);
                A2_Re(k,:) = t_sample_store.ANTENNA2(1,:);
                A2_Img(k,:) = t_sample_store.ANTENNA2(2,:);
            elseif(strcmp(mode,'SAR')) 
               
                t_sample_store = struct('ANTENNA1',fread(fid, [2,t_header.EXTENDED_RW_SAMPLES],'float32'));
                A1_Re(k,:) = t_sample_store.ANTENNA1(1,:);
                A1_Img(k,:) = t_sample_store.ANTENNA1(2,:);
            end
            complex_wvf(k,:) = complex(A1_Re(k,:),A1_Img(k,:));
            
            wvf_t(k,:)=fft(complex_wvf(k,:));
            wvf(k,:)         = CalcSpectraPulse2(wvf_t(k,:)', 1);
            wvf(k,:)          = fftshift(wvf(k,:));
            avg_Power(k)=max(abs(wvf(k,:)));
          end
          stack_power(i+1) = sum(avg_Power);
        
    end
    
     fclose(fid);
     
     if(strcmp(mode,'SIN'))
        [~,pos_of_max]=max(stack_power(150:200)); %Do not go far than stack 200, false stack detection
         pos_of_max=pos_of_max+149;
     elseif(strcmp(mode,'SAR'))
        [~,pos_of_max]=max(stack_power(:)); %Do not go far than stack 200, false stack detection
     end
     pos_of_max=pos_of_max-2;
  end
   
   
  
  
   

   




