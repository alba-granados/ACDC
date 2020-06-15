global N_samples_sar_chd zp_fact_range_cnf T0_chd c_cst
% lla2kml('./results/CS_OFFL_SIR_SIN_FR_20140818T053651_20140818T053937_B001.kml', L1BS.lat_surf, L1BS.lon_surf, L1BS.alt_surf,'.')
i_surf=139;
%% L1B-S
beams=find(L1BS.ProcessID_beam(i_surf,:)==58);
zp_fact_range_cnf=2;
for i_beam = 1:L1BS.N_beams_stack(i_surf)
    for i_sample=1:N_samples_sar_chd*zp_fact_range_cnf
        shift(i_beam,i_sample) = L1BS.shift_coarse(i_surf,i_beam)+i_sample; 
        shift_mod(i_beam,i_sample) =  ceil(shift(i_beam,i_sample)/N_samples_sar_chd/zp_fact_range_cnf);
     end
end
aux = 1./zeros(L1BS.N_beams_stack(i_surf),N_samples_sar_chd*zp_fact_range_cnf*max(max(shift_mod))+N_samples_sar_chd*zp_fact_range_cnf/2);
aux2 = aux;
aux3 = aux;

aux(:,end-N_samples_sar_chd*zp_fact_range_cnf+1:end) =(abs(fftshift(fft(squeeze(L1BS.beams_surf(i_surf,1:L1BS.N_beams_stack(i_surf),:))',N_samples_sar_chd*zp_fact_range_cnf),1))');
aux(:,end-N_samples_sar_chd*zp_fact_range_cnf+1:end) =squeeze(L1BS.beams_rng_cmpr(i_surf,1:L1BS.N_beams_stack(i_surf),:));

aux2(:,end-N_samples_sar_chd*zp_fact_range_cnf+1:end) = squeeze(L1BS.phase_diff(i_surf,1:L1BS.N_beams_stack(i_surf),:));
stack_wd_corr = 1./zeros(L1BS.N_beams_stack(i_surf),N_samples_sar_chd*zp_fact_range_cnf*max(max(shift_mod))+N_samples_sar_chd*zp_fact_range_cnf/2);
phase_wd_corr = stack_wd_corr;
coh_wd_corr = stack_wd_corr;
for i_beam = 1:L1BS.N_beams_stack(i_surf)
    shift(i_beam) = (-L1BS.shift(i_surf,i_beam));%-mod(L1BS.shift(i_surf,i_beam),N_samples_sar_chd*zp_fact_range_cnf));
    stack_wd_corr(i_beam,:) = circshift(aux(i_beam,:),[0,round(shift(i_beam)*zp_fact_range_cnf)]);
%     phase_wd_corr(i_beam,:) = circshift(aux2(i_beam,:),[0,round(shift(i_beam)*zp_fact_range_cnf)]);
%     coh_wd_corr(i_beam,:)   = circshift(aux3(i_beam,:),[0,round(shift(i_beam)*zp_fact_range_cnf)]);
    
end

for i_beam = 1:L1BS.N_beams_stack(i_surf)
    aux_(i_beam) = sum(squeeze(L1BS.beams_rng_cmpr(i_surf,i_beam,:)));
end

figure; imagesc(fliplr(stack_wd_corr(beams,:))); colormap('jet');
figlabels('Samples','Beams','','Stack at the border of the lake',16);
figure; imagesc(fliplr(phase_wd_corr(beams,:))); colormap('jet');


aux2=fliplr(stack_wd_corr(beams,:));
for i_sample = 1:length(stack_wd_corr)
    wfmL1(i_sample)= mean((aux2(isfinite((aux2(:,i_sample))),i_sample)));
end

for i_sample = 1:N_samples_sar_chd*zp_fact_range_cnf
    wfmL2(i_sample)=mean((aux2(isfinite((aux2(:,i_sample))),i_sample)));
end
for i_sample = 1+N_samples_sar_chd*zp_fact_range_cnf:N_samples_sar_chd*zp_fact_range_cnf*2
    wfmL3(i_sample)=mean((aux2(isfinite((aux2(:,i_sample))),i_sample)));
end
for i_sample = 1+N_samples_sar_chd*zp_fact_range_cnf*2:N_samples_sar_chd*zp_fact_range_cnf*3
    wfmL4(i_sample)=mean((aux2(isfinite((aux2(:,i_sample))),i_sample)));
end

figure; subplot(2,1,1);plot(wfmL1); 
set(gca,'XLim',[1 length(stack_wd_corr)]);set(gca,'YLim',[1 600]);
figlabels('Samples','','','Possible L1B waveforms at the border of the lake',16);
subplot(2,3,4);plot(wfmL2,'r'); 
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 600]);
figlabels('Samples','','','',16);
subplot(2,3,5);plot(1:N_samples_sar_chd*zp_fact_range_cnf,wfmL3(1+N_samples_sar_chd*zp_fact_range_cnf:N_samples_sar_chd*zp_fact_range_cnf*2),'r');
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 600]);
figlabels('Samples','','','',16);
subplot(2,3,6);plot(1:N_samples_sar_chd*zp_fact_range_cnf,wfmL4(1+N_samples_sar_chd*zp_fact_range_cnf*2:N_samples_sar_chd*zp_fact_range_cnf*3),'g');
set(gca,'XLim',[1 N_samples_sar_chd*zp_fact_range_cnf]);set(gca,'YLim',[1 600]);
figlabels('Samples','','','',16);












for i_beam = 1:L1BS.N_beams_stack(i_surf)
    for i_sample=1:N_samples_sar_chd*zp_fact_range_cnf
        stack_wd_corr(i_beam,i_sample+(N_samples_sar_chd*zp_fact_range_cnf/2)*(shift_mod(i_beam,i_sample)-1)) = (L1BS.beams_rng_cmpr(i_surf,i_beam,i_sample));
    end
    
end

max_wd=max(win_del_surf_beam);
min_wd=min(win_del_surf_beam);
N_windows = ceil((max_wd-min_wd)/N_samples_sar_chd/zp_fact_range_cnf)+1.5;
[~,ref_surf] = min(win_del_surf_beam);
stack_wd_corr = 1./zeros(L1BS.N_beams_stack(i_surf),N_samples_sar_chd*zp_fact_range_cnf*N_windows);
aux =  1./zeros(L1BS.N_beams_stack(i_surf),N_samples_sar_chd*zp_fact_range_cnf*N_windows);

for i_beam = 1:L1BS.N_beams_stack(i_surf)
    aux(i_beam,1:N_samples_sar_chd*zp_fact_range_cnf)  = squeeze(L1BS.beams_rng_cmpr(i_surf,i_beam,:))./max(squeeze(L1BS.beams_rng_cmpr(i_surf,i_beam,:)));
    wd_shift(i_beam) = ((win_del_surf_beam(i_beam)-win_del_surf_beam(ref_surf)));
    wd_shift2(i_beam) = wd_shift(i_beam)-mod(wd_shift(i_beam),N_samples_sar_chd*zp_fact_range_cnf);
    stack_wd_corr(i_beam,:) = circshift(aux(i_beam,:),[0,round(wd_shift2(i_beam)*zp_fact_range_cnf)]);
    
end


figure;imagesc(1:N_samples_sar_chd*zp_fact_range_cnf,1:60,(stack_wd_corr)); colormap('jet');
figlabels('Samples','Beams','','Stack at the border of the lake',16);




window_init=L1BS.alt_surf(i_surf)-(min(win_del_surf_beam)*T0_chd)*c_cst/2;    
x_axis = window_init-T0_chd*c_cst/2:-T0_chd*c_cst/2:window_init-512*zp_fact_range_cnf*N_windows*T0_chd*c_cst/2;
y_axis = L1A.lat_sar_sat(L1BS.burst_index(i_surf,1:L1BS.N_beams_stack(i_surf)));



figure; k=surf(y_axis,x_axis,stack_wd_corr');set(k, 'edgecolor','none');view(90,-90);


subplot(2,2,3);


%% L1B
N_windows = ceil((max(L1BS.win_delay_surf)-min(L1BS.win_delay_surf))./T0_chd/N_samples_sar_chd/zp_fact_range_cnf);
L1B.wfm_cor_i2q2_sar_ku_wdcorr = zeros(L1BS.N_total_surf_loc,N_samples_sar_chd*zp_fact_range_cnf*N_windows);
aux =  zeros(L1BS.N_total_surf_loc,N_samples_sar_chd*zp_fact_range_cnf*N_windows);
[~,ref_surf] = min(L1BS.win_delay_surf);

for i_surf = 1:L1BS.N_total_surf_loc
    aux(i_surf,1:N_samples_sar_chd*zp_fact_range_cnf)  = L1B.wfm_cor_i2q2_sar_ku(i_surf,:);
    wd_shift(i_surf) = ((L1BS.win_delay_surf(i_surf)-L1BS.win_delay_surf(ref_surf))-(L1BS.alt_sat(i_surf)-L1BS.alt_sat(ref_surf))*2/c_cst) / T0_chd;
    L1B.wfm_cor_i2q2_sar_ku_wdcorr(i_surf,:) = circshift(aux(i_surf,:),[0,round(wd_shift(i_surf)*zp_fact_range_cnf)])./max(aux(i_surf,:));
end

L1BS.win_delay_surf_aligned = L1BS.win_delay_surf-wd_shift*T0_chd;
figure;  k=surf(L1B.wfm_cor_i2q2_sar_ku_wdcorr);set(k, 'edgecolor','none');view(90,-90);
