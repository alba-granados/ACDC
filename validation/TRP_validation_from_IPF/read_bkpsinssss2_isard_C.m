% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% CryoSat 2 calibration over transponders
% 
% This code implements the algorithm as described in the
% ISARD_ESA_CR2_TRP_CAL_DPM_030 2.b of 26/05/2011
%
% ---------------------------------------------------------
% READ_BKPSINSSSS2: reads the stack data files provided by ESA
% 
% Calling
%   [out_stack, dates] = read_bkpsinssss2_isard(filename)
%
% Input
%   filename    : stack file to be read
%   pass        : number of the dataset to be analysed
%
% Outputs
%   out_stack   : stack that corresponds to the Transponder
%   dates       : all the times from the stack
%   pos         : stack position of TRP
%
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
%           Josep Montolio / Pildo Labs
%           Pablo García / isardSAT
%           Albert García / isardSAT
%
% Reviewer: Mònica Roca / isardSAT
%
% Last revision: Mònica Roca / isardSAT (26/05/11)
%
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [out_stack, dates] = read_bkpsinssss2_isard_C(filename, pass, record_num, pos_of_max, antenna, plots_flag, config_selec,str_selec, resultPath)
% val_of_max =0;
% iStack=0;
  fid = fopen(filename,'r', 'b')	;


  status = fclose(fid);
  fid = fopen(filename,'r', 'b');

for i = 0:record_num*20 %we will read from start to position of the maximum
    
   
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
              'CONTRIB_BEAMS_BEFORE_WEIGHTING', fread(fid,1,'int32'));
        if(isempty(t_header.day))
            
            break;  
        end
        hour_st = fix(t_header.seconds/3600.);
        minute_st = t_header.seconds-hour_st*3600;
        min_st = fix(minute_st / 60.);
        sec_st = minute_st(1)-min_st*60;
%         fprintf(fid2,'%d\t\t %d\t\t %d\t\t %f\t\t\n', i, hour_st, min_st, sec_st) ;   

        dates(i+1) = t_header.seconds;
    
    if (i ~= pos_of_max)%%&&(i ~= pos_of_max+1)&&(i ~= pos_of_max-1)
     fread(fid, 1, 'uint32');
     %seek of the t_ocog struck
     fseek(fid,(13+t_header.CONTRIB_BEAMS+t_header.RW_SAMPLES+t_header.SUB_STACK_BEAMS+t_header.SUB_STACK_BEAMS+100+t_header.SUB_STACK_BEAMS*t_header.RW_SAMPLES)*4,'cof');
    %seek of the t_sample_store struck
    if(t_header.MODE==0)% SAR MODE
        fseek(fid, t_header.CONTRIB_BEAMS*(3*4+2*2+3*4+1*4+5*4+2*t_header.EXTENDED_RW_SAMPLES*4),'cof');     
    elseif(t_header.MODE==1)% SARIN MODE
        fseek(fid, t_header.CONTRIB_BEAMS*(3*4+2*2+3*4+1*4+5*4+4*t_header.EXTENDED_RW_SAMPLES*4),'cof');       
    end

     
    else                 
        t_ocog= struct(  'FLAG', fread(fid,1,'uint32'), ...
                         'SUM_OF_SQUARES',fread(fid,1,'float32'), ...
                         'SUM_OF_QUADS',fread(fid,1,'float32'), ...
                         'AMPLITUDE',fread(fid,1,'float32'), ...
                         'WIDTH',fread(fid,1,'float32'), ...
                         'TWICE_WIDTH_SQ',fread(fid,1,'float32'), ...
                         'HALF_WIDTH',fread(fid,1,'float32'),  ...
                         'CENTRE',fread(fid,1,'float32'), ...
                         'SKEW',fread(fid,1,'float32'), ...
                         'KURT',fread(fid,1,'float32') , ...
                         'THRESHOLD',fread(fid,1,'float32'),...
                         'SUM_SQUARED',fread(fid,1,'float32'),...
                         'MEAN',fread(fid,1,'float32'),...
                         'STDDEV',fread(fid,1,'float32') ,...
                         'REAL_BEAM_WEIGHTS',fread(fid,t_header.CONTRIB_BEAMS, 'float32'), ...
                         'SUM_OF_WEIGHTS',fread(fid,t_header.RW_SAMPLES, 'float32'), ...
                         'RANGE_STACKED_POWER',fread(fid,t_header.SUB_STACK_BEAMS, 'float32'), ...
                         'CUMULATIVE_POWER',fread(fid,t_header.SUB_STACK_BEAMS, 'float32'), ...
                         'PERCENTILE',fread(fid,100,'float32'), ...
                         'SUB_STACK_POWER',fread(fid,t_header.SUB_STACK_BEAMS*t_header.RW_SAMPLES,'float32') );
   

          for k = 1:t_header.CONTRIB_BEAMS
              %   Progress bar
                customText = 'Reading TRP Stack...';
                percentageDone =  k / t_header.CONTRIB_BEAMS;
%                 stopBar= progressbar(percentageDone, 0, customText);
%                 if (stopBar) 
%                     break; 
%                 end
              if(t_header.MODE==0)% SAR MODE
                  t_sample_store = struct('FBR_DAYS', fread(fid,1,'uint32'), ...
                                         'FBR_SECONDS', fread(fid,1,'uint32'), ...
                                         'FBR_MICROSECONDS',fread(fid,1,'uint32'), ...
                                         'BURST',fread(fid, 1, 'int16'),...
                                         'BEAM_NUM',fread(fid, 1, 'int16'), ...
                                         'ANGLE',fread(fid, 1, 'float32'), ...
                                         'BORESIGHT_ANGLE',fread(fid, 1, 'float32'), ...
                                         'LOOK_ANGLE',fread(fid, 1, 'float32'), ... % new field baseline C
                                         'DELTA_COUNT',fread(fid, 1, 'int32'), ...
                                         'CENTRE_ANGLE',fread(fid, 1, 'float32'), ...
                                         'DOPPLER_CORR',fread(fid, 1, 'float32'), ...
                                         'DOPPLER_COARSE_CORR',fread(fid, 1, 'float32'), ...
                                         'DOPPLER_FINE_CORR',fread(fid, 1, 'float32'), ... 
                                         'SLANT_R',fread(fid, 1, 'float32'), ...
                                         'ANTENNA1',fread(fid, [2,t_header.EXTENDED_RW_SAMPLES],'float32'));
                     for j = 1:t_header.EXTENDED_RW_SAMPLES   
                            avg_wvf_real (k,j) = (t_sample_store.ANTENNA1(1,j));
                            avg_wvf_img (k,j) = (t_sample_store.ANTENNA1(2,j));
                            avg_wvf (k,j) = module([avg_wvf_real(k,j) ,  avg_wvf_img(k,j)]);
                            waveform(k,j) = module(t_sample_store.ANTENNA1(:,j));

                      end
                  
              elseif(t_header.MODE==1) % SARIN MODE
                  
                  t_sample_store = struct('FBR_DAYS', fread(fid,1,'uint32'), ...
                                     'FBR_SECONDS', fread(fid,1,'uint32'), ...
                                     'FBR_MICROSECONDS',fread(fid,1,'uint32'), ...
                                     'BURST',fread(fid, 1, 'int16'),...
                                     'BEAM_NUM',fread(fid, 1, 'int16'), ...
                                     'ANGLE',fread(fid, 1, 'float32'), ...
                                     'BORESIGHT_ANGLE',fread(fid, 1, 'float32'), ...
                                     'LOOK_ANGLE',fread(fid, 1, 'float32'), ... % new field baseline C
                                     'DELTA_COUNT',fread(fid, 1, 'int32'), ...
                                     'CENTRE_ANGLE',fread(fid, 1, 'float32'), ...
                                     'DOPPLER_CORR',fread(fid, 1, 'float32'), ...
                                     'DOPPLER_COARSE_CORR',fread(fid, 1, 'float32'), ...
                                     'DOPPLER_FINE_CORR',fread(fid, 1, 'float32'), ... 
                                     'SLANT_R',fread(fid, 1, 'float32'), ...
                                     'ANTENNA1',fread(fid, [2,t_header.EXTENDED_RW_SAMPLES],'float32'),...
                                     'ANTENNA2',fread(fid, [2,t_header.EXTENDED_RW_SAMPLES],'float32' ));
                                 
                  for j = 1:t_header.EXTENDED_RW_SAMPLES   
                        avg_wvf_real (k,j) = (t_sample_store.ANTENNA1(1,j) + t_sample_store.ANTENNA2(1,j))/2.;
                        avg_wvf_img (k,j) = (t_sample_store.ANTENNA1(2,j) + t_sample_store.ANTENNA2(2,j))/2.;
                        avg_wvf (k,j) = module([avg_wvf_real(k,j) ,  avg_wvf_img(k,j)]);
                        waveform(k,j) = module(t_sample_store.ANTENNA1(:,j));
                        waveform2(k,j) = module(t_sample_store.ANTENNA2(:,j));                  
                  end
              end
              
              slantrange(k) =  t_sample_store.SLANT_R;
              dopplercoarsecorr(k) = t_sample_store.DOPPLER_COARSE_CORR; 
              dopplercorr(k) = t_sample_store.DOPPLER_CORR;
              dopplerfinecorr(k) = t_sample_store.DOPPLER_FINE_CORR;
              A1_Re(k,:) = t_sample_store.ANTENNA1(1,:);
              A1_Img(k,:) = t_sample_store.ANTENNA1(2,:);
              
              if(t_header.MODE==1)
                A2_Re(k,:) = t_sample_store.ANTENNA2(1,:);
                A2_Img(k,:) = t_sample_store.ANTENNA2(2,:);
              
              end
              
              days_FBR(k) = t_sample_store.FBR_DAYS;
              seconds_FBR(k) = t_sample_store.FBR_SECONDS;
              microseconds_FBR(k) = t_sample_store.FBR_MICROSECONDS;
              boresight_angle(k) = t_sample_store.BORESIGHT_ANGLE;
              angle(k) = t_sample_store.ANGLE;
              centre_angle(k) = t_sample_store.CENTRE_ANGLE;
              iburst(k)=t_sample_store.BURST;
              ibeam(k)=t_sample_store.BEAM_NUM;
              
              
              


              
              end
         end
   
     
   if (i == pos_of_max)%%||(i == pos_of_max+1)||(i == pos_of_max-1)
        
         if(t_header.MODE==1)
       out_stack = struct('t_header', t_header,...
                      'waveform', waveform, ...
                      'waveform2', waveform2, ...
                      'avg_waveform', avg_wvf, ...
                      'avg_wvf_real', avg_wvf_real, ...
                      'avg_wvf_img', avg_wvf_img, ...
                      'slant_range', slantrange, ...
                      'doppler_coarse_corr', dopplercoarsecorr, ... 
                      'doppler_corr', dopplercorr, ...
                      'doppler_fine_corr', dopplerfinecorr, ...
                      'A1_Re', A1_Re, ...
                      'A1_Img', A1_Img, ...
                      'A2_Re', A2_Re, ...
                      'A2_Img', A2_Img, ...
                      'days_FBR', days_FBR, ...
                      'seconds_FBR', seconds_FBR, ...
                      'microseconds_FBR', microseconds_FBR, ...
                      'boresight_angle',boresight_angle, ...
                      'centre_angle',centre_angle, ...
                      'angle',angle,...
                      'iburst',iburst,...
                      'ibeam',ibeam);
                  
         elseif(t_header.MODE==0)
            out_stack = struct('t_header', t_header,...
                      'waveform', waveform, ...
                      'avg_waveform', avg_wvf, ...
                      'avg_wvf_real', avg_wvf_real, ...
                      'avg_wvf_img', avg_wvf_img, ...
                      'slant_range', slantrange, ...
                      'doppler_coarse_corr', dopplercoarsecorr, ... 
                      'doppler_corr', dopplercorr, ...
                      'doppler_fine_corr', dopplerfinecorr, ...
                      'A1_Re', A1_Re, ...
                      'A1_Img', A1_Img, ...
                      'days_FBR', days_FBR, ...
                      'seconds_FBR', seconds_FBR, ...
                      'microseconds_FBR', microseconds_FBR, ...
                      'boresight_angle',boresight_angle, ...
                      'centre_angle',centre_angle, ...
                      'angle',angle,...
                      'iburst',iburst,...
                      'ibeam',ibeam); 
         end
         if(plots_flag)         
             
            figure; mesh(waveform);

             colormap('hot');
%              set(kk, 'edgecolor','none');axis xy;%view(20,40);
             set(gca,'XLim',[1 size(waveform,2)], 'YLim',[1 size(waveform,1)],'FontSize',20) % --> eixos
             xlab = get(gca,'XLabel'); set(xlab,'String','Samples','FontSize',18,'FontName','Arial')
             ylab = get(gca,'YLabel'); set(ylab,'String','Beams','FontSize',18,'FontName','Arial')
             zlab = get(gca,'ZLabel'); set(zlab,'String','FFT p.u.','FontSize',18,'FontName','Arial')

             title(sprintf('Waveform Antenna %i. Pass: %i Pos of max: %i ',1, pass,  i ))

            %          saveas(gcf, sprintf('%s%i_%s%s_Stack_%i_ant%i.fig',resultPath,pass,config_selec,str_selec,i,antenna) , 'fig');
%              saveas(gcf, sprintf('%s%i_%s%s_Stack_%i_ant%i.png',resultPath,pass,config_selec,str_selec,i,1) , 'png');
             close all
             figure; imagesc(waveform);
             
             colormap('hot');
             set(gca,'XLim',[1 size(waveform,2)], 'YLim',[1 size(waveform,1)],'FontSize',20) % --> eixos
             xlab = get(gca,'XLabel'); set(xlab,'String','Samples','FontSize',18,'FontName','Arial')
             ylab = get(gca,'YLabel'); set(ylab,'String','Beams','FontSize',18,'FontName','Arial')
             zlab = get(gca,'ZLabel'); set(zlab,'String','FFT p.u.','FontSize',18,'FontName','Arial')
            title(sprintf('Waveform Antenna %i. Pass: %i Pos of max: %i ',antenna, pass,  i ))
            %          saveas(gcf, sprintf('%s%i_%s%s_StackXY_%i_ant%i.fig',resultPath,pass,config_selec,str_selec,i,antenna) , 'fig');
%              saveas(gcf, sprintf('%s%i_%s%s_StackXY_%i_ant%i.png',resultPath,pass,config_selec,str_selec,i,1) , 'png');

             close all
             if(antenna==2)
                    figure; mesh(waveform2);

                    colormap('hot');
%                     set(kk, 'edgecolor','none');axis xy;%view(20,40);
                    set(gca,'XLim',[1 size(waveform,2)], 'YLim',[1 size(waveform,1)],'FontSize',20) % --> eixos
                    xlab = get(gca,'XLabel'); set(xlab,'String','Samples','FontSize',18,'FontName','Arial')
                    ylab = get(gca,'YLabel'); set(ylab,'String','Beams','FontSize',18,'FontName','Arial')
                    zlab = get(gca,'ZLabel'); set(zlab,'String','FFT p.u.','FontSize',18,'FontName','Arial')

                    title(sprintf('Waveform Antenna %i. Pass: %i Pos of max: %i ',antenna, pass,  i ))

                    %          saveas(gcf, sprintf('%s%i_%s%s_Stack_%i_ant%i.fig',resultPath,pass,config_selec,str_selec,i,antenna) , 'fig');
                    saveas(gcf, sprintf('%s%i_%s%s_Stack_%i_ant%i.png',resultPath,pass,config_selec,str_selec,i,antenna) , 'png');
                    close all
                    figure; imagesc(waveform2);
                    colormap('hot');
                   
                    set(gca,'XLim',[1 size(waveform,2)], 'YLim',[1 size(waveform,1)],'FontSize',20) % --> eixos
                    xlab = get(gca,'XLabel'); set(xlab,'String','Samples','FontSize',18,'FontName','Arial')
                    ylab = get(gca,'YLabel'); set(ylab,'String','Beams','FontSize',18,'FontName','Arial')
                    zlab = get(gca,'ZLabel'); set(zlab,'String','FFT p.u.','FontSize',18,'FontName','Arial')
                    title(sprintf('Waveform Antenna %i. Pass: %i Pos of max: %i ',antenna, pass,  i ))

                    %          saveas(gcf, sprintf('%s%i_%s%s_StackXY_%i_ant%i.fig',resultPath,pass,config_selec,str_selec,i,antenna) , 'fig');
                    saveas(gcf, sprintf('%s%i_%s%s_StackXY_%i_ant%i.png',resultPath,pass,config_selec,str_selec,i,antenna) , 'png');

                    close all 
             end
         
         end
%         stack_max(i+1) = max(max(waveform,[],2)); 
    end
  end
   
   
%    fclose(fid2);
   fclose(fid);
   
%     for i=1:length(wvf_stack)
%         sumatory(i)= sum(sum(wvf_stack(i).waveform(:,291:295)));
%     end
%     plot(sumatory)
%     [~,pos]=max(sumatory(i))
   
end



