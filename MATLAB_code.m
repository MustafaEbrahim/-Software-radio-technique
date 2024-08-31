clear all
close all
if(exist('a6','var'))
clf(figure(1),figure(2));
delete([a1 a2 a3 a4  a5 a6]);
end
%% General Notes
frame_length =100; %Required no. of tranmitted bits
Ensemble_width = 500; %Number of realizations
DAC_samples = 7; % Dac samples = Dac sampling freq / Bit rate = (70 ms)/(10 ms)
cond=1;
while(cond)
A = input(" Enter the amplitude of the line code : "); % Rails
if(A>0 && A<10)
    cond=0;
    else 
        fprintf(" Wrong Input A in ]0,10] , please try again ..... \n");
    end
end
final_frame_length = frame_length * DAC_samples; %The frame length after the DAC
Rs = (10^3)/70; %(1/Ts) = 1/70ms --> Symbol rate 
Rb = (10^3)/10; %(1/Tb = n/Ts) = 1/10ms --> bit rate
%% Generation of the general random Delay & Tx stream (not mapped yet)
% Add extra bit is used to introduce the delay given by the last transmitted frame
tx = randi([0 1],Ensemble_width,frame_length+1);
tx_DAC_out = repelem(tx,1,DAC_samples);  
% Delay the realizations
delay_period = randi([0 DAC_samples-1],1,Ensemble_width);
tx_delayed = zeros(Ensemble_width,size(tx,2)*DAC_samples);
for i=1:Ensemble_width
    tx_delayed(i,:) = circular_shift(tx_DAC_out(i,:),delay_period(i));
end
% pickup the shifted frame
for i=1:Ensemble_width
    tx_final = tx_delayed(:,delay_period(i)+1:delay_period(i)+final_frame_length);
end
%% Choose the required Signaling
Signaling_line_codes = {'Unipolar NRZ','Polar NRZ','Polar RZ'};  
[indicies,check] = listdlg('ListString',Signaling_line_codes,'promptstring','Select the line code to continue (Unipolar Signaling by default)','SelectionMode','single');
if (~check)
    indicies = 1 ;
end
switch cell2mat(Signaling_line_codes(indicies))
    case 'Unipolar NRZ'
       Signaling_mapped = tx_final*A;  % Unipolar NRZ Signaling Line Code
    case 'Polar NRZ'
       Signaling_mapped = (2*tx_final-1)*A;  % polar NRZ Signaling Line Code
    case 'Polar RZ'
       Signaling_mapped = (2*tx_final-1)*A; % polar RZ (Take NRZ and replace some bits with Zeros)
       for k=5:7:final_frame_length-2
            for l = k:k+2
                Signaling_mapped(:,l)=0;   
            end
       end
end 
%% Checking the sationarity using statistical (Mean/ACF)
Signaling_stat_mean = Stat_mean(Signaling_mapped,final_frame_length); 
Signaling_stat_ACF_single_sided_positive = stat_acf(Signaling_mapped,final_frame_length);
Signaling_stat_ACF_single_sided_negative = fliplr(Signaling_stat_ACF_single_sided_positive);
Signaling_stat_ACF = [Signaling_stat_ACF_single_sided_negative(:,1:end-1) Signaling_stat_ACF_single_sided_positive];
tau_domain_length = length(Signaling_stat_ACF);
tau_domain = (-((tau_domain_length-1)/2):((tau_domain_length-1)/2)).*(Rs/Rb);
time_domain = 0:final_frame_length-1;
%{ 
       1- To see if the mean is constant , we can check that :
         * Unipolar :   
           Average of the mean across the time instants ~ A/2 & Std  ~ 0
         * Polar NRZ:
           Average of the mean across the time instants ~ 0 & Std  ~ 0
         * Polar RZ:
           Average of the mean across the time instants ~ 0 & Std ~ 0
       2- We can check the statistic ACF by :
        * Unipolar :
            ACF @(tau=0) ~ A^2/2 & DC ~ A^2/4 
        *Polar NRZ: 
            ACF @(tau=0) ~ A^2 & DC ~ 0 
        *Polar RZ : 
            ACF @(tau=0) ~ A^2 *4/7 & DC ~ 0 (4/7 = Duty Cycle)
%}
str_stat_mean_Signaling = sprintf(" Average = %.3f \nStd = %.3f",mean(Signaling_stat_mean),std(Signaling_stat_mean));
str_stat_ACF_Signaling = sprintf(" DC = %.3f \nR(0) = %.3f",mean(Signaling_stat_ACF),Signaling_stat_ACF(tau_domain==0));
%% plot the statistic mean & ACF
figure(1);
plot(time_domain,Signaling_stat_mean);
xlabel('Time (msec) ');
ylabel('E\{X(t)\}')
set(gca,'XAxisLocation','origin');
title(sprintf("Stat %s Signaling Mean",cell2mat(Signaling_line_codes(indicies))));
a1=annotation('textbox',[0.4 0.62 0.3 0.3],'String',str_stat_mean_Signaling,'FitBoxToText','on');
figure(2);
plot(tau_domain,Signaling_stat_ACF);
xlabel('\tau *(Rs/Rb) (ms)','FontSize',12);
ylabel('R(\tau)','FontSize',14);
set(gca,'XAxisLocation','origin');
title(sprintf("Stat %s Signaling ACF",cell2mat(Signaling_line_codes(indicies))));
a2=annotation('textbox',[0.7 0.62 0.3 0.3],'String',str_stat_ACF_Signaling,'FitBoxToText','on');
figure();
zoom = final_frame_length+[-50 50];
plot(tau_domain(zoom(1):zoom(2)),Signaling_stat_ACF(zoom(1):zoom(2)));
a2=annotation('textbox',[0.7 0.62 0.3 0.3],'String',str_stat_ACF_Signaling,'FitBoxToText','on');
title(sprintf("Stat %s Signaling ACF (zoomed version)",cell2mat(Signaling_line_codes(indicies))));
xlabel('\tau *(Rs/Rb) (ms)','FontSize',12);
ylabel('R(\tau)','FontSize',14);
set(gca,'XAxisLocation','origin');
%% Checking the ergodicity using time(mean/ACF):
R_domain = 0:Ensemble_width-1;
Signaling_time_mean = time_mean(Signaling_mapped,Ensemble_width);
Signaling_time_ACF_singlesided_positive =  time_acf(Signaling_mapped,Ensemble_width);
Signaling_time_ACF_singlesided_negative = fliplr(Signaling_time_ACF_singlesided_positive);
Signaling_time_ACF = [Signaling_time_ACF_singlesided_negative(:,1:end-1) Signaling_time_ACF_singlesided_positive];
r_length = length(Signaling_time_ACF);
r_domain = -(r_length-1)/2:(r_length-1)/2; 
str_time_mean_Signaling = sprintf("Average = %.3f \nStd = %.3f",mean(Signaling_time_mean),std(Signaling_time_mean));
str_time_ACF_Signaling = sprintf("DC = %.3f \nR(0) = %.3f",mean(Signaling_time_ACF),Signaling_time_ACF(r_domain==0));
%% plot the time mean & ACF
figure();
plot(R_domain,Signaling_time_mean,'LineWidth',0.25);
xlabel('Realization');
ylabel('<X(t)>')
set(gca,'XAxisLocation','origin');
title(sprintf("Time %s Signaling Mean",cell2mat(Signaling_line_codes(indicies))));
a3=annotation('textbox',[0.4 0.62 0.3 0.3],'String',str_time_mean_Signaling,'FitBoxToText','on');
figure();
plot(r_domain,Signaling_time_ACF);
xlabel('\tau','FontSize',16);
ylabel('<X(r).X(r+\tau)>','FontSize',14);
set(gca,'XAxisLocation','origin');
title(sprintf("Time %s Signaling Autocorrelation Function",cell2mat(Signaling_line_codes(indicies))));
a4=annotation('textbox',[0.7 0.62 0.3 0.3],'String',str_time_ACF_Signaling,'FitBoxToText','on');
figure();
zoom= round(Ensemble_width + ([-50 50]*(Rs/Rb)));
plot(r_domain(1,zoom(1):zoom(2)),Signaling_time_ACF(1,zoom(1):zoom(2)));
xlabel('\tau','FontSize',16);
ylabel('<X(r).X(r+\tau)>','FontSize',14);
set(gca,'XAxisLocation','origin');
a4=annotation('textbox',[0.7 0.62 0.3 0.3],'String',str_time_ACF_Signaling,'FitBoxToText','on');
title(sprintf("Time %s Signaling Autocorrelation Function (zoomed version)",cell2mat(Signaling_line_codes(indicies))));
%% PSD
freq_domain = (-((tau_domain_length-1)/2):((tau_domain_length-1)/2)).*(Rb/Rs/(tau_domain_length-1)); % Normalized frequency scale (k) wrt to the sampling
switch cell2mat(Signaling_line_codes(indicies))
    case 'Unipolar NRZ'
       Signaling_psd = abs(fftshift(fft(Signaling_stat_ACF)))/tau_domain_length;
       zero_crossing = psd_zero_cross(Signaling_psd,Rs*freq_domain,cell2mat(Signaling_line_codes(indicies)),Rs);
       str_psd_Signaling=sprintf('S(0) = %.3f \nZero crossing @ f = %.3f',Signaling_psd(freq_domain==0),zero_crossing);
       figure();
       plot(freq_domain,Signaling_psd,'LineWidth',0.2);
       xlabel('Normalized frequency scale : k = (Rs/Rb) * F ');
       ylabel('|S(f)|');
       set(gca,'YAxisLocation','origin');
       a5 = annotation('textbox',[0.55 0.52 0.3 0.3],'QString',str_psd_Signaling,'FitBoxToText','on');
    case {'Polar NRZ','Polar RZ'} 
             Signaling_psd = abs(fftshift(fft(Signaling_stat_ACF)))/final_frame_length;
             Mag_Norm = [1 Rb/Rs];
             zero_crossing = psd_zero_cross(Signaling_psd,Rs*freq_domain,cell2mat(Signaling_line_codes(indicies)),Rs);
       for p=1:2
           if(p==1)
             psd_title = sprintf('%s Signaling PSD Un-normalized magnitude',cell2mat(Signaling_line_codes(indicies)));
             psd_ylabel = sprintf('|S(f)|');Q
           else
            psd_title = sprintf('%s Signaling PSD Normalized magnitude ',cell2mat(Signaling_line_codes(indicies)));              
            psd_ylabel = sprintf('|S(f)|*(Rb/Rs)|');
           end
           figure(); 
           plot(freq_domain,Signaling_psd*Mag_Norm(p),'LineWidth',0.2);
           xlabel('Normalized frequency scale : k = (Rs/Rb) * F ');
           ylabel(psd_ylabel);
           set(gca,'YAxisLocation','origin');
           title(psd_title);
           str_psd_Signaling=sprintf('S(0) = %.3f \nZero crossing @ f = %.3f',Signaling_psd(freq_domain==0)*Mag_Norm(p),zero_crossing);
           a5 = annotation('textbox',[0.55 0.52 0.3 0.3],'String',str_psd_Signaling,'FitBoxToText','on');
                if(indicies==2)
                    a6=annotation('ellipse','Position',[0.59 0.1 0.05 0.05],'Color','r');
                else 
                    a6=annotation('ellipse','Position',[0.69 0.1 0.05 0.05],'Color','r');
                end
      end
end
%% Aiding Functions
function z = circular_shift(x,shift_val) %circular shift the realizations with the delay
len = length(x);
z= zeros(1,len);
    if(shift_val~=0)
        start = len-(shift_val-1);
        i = start;
        j=1;
        while (j<=len)
            z(j)=x(i);
            j=j+1;
            i=i+1;
            if(i>len)
               i=1;
            end
        end
    else
        z = x;
    end
    
end
function z = average(x)
len = length(x);
sum = 0;
for i= 1:len
    sum = sum + x(i);
end 
z = sum/len;
end
function [y]= Stat_mean(x,frame_length)
y = zeros(1,frame_length);
for t=1:frame_length
y(t)= average(x(:,t));
end
end
function [z] = stat_acf(x,frame_length) %The code exercises the +ve & -ve sides of the ACF so it implicitly checks if it is even or not
y=zeros(frame_length,frame_length);
z = zeros(1,frame_length);   
tau_hash = zeros(1,frame_length);
for t1 = 1 : frame_length
  tau=0;
    while(1)
        t2=t1+tau;
        y(t1,abs(tau)+1) = average(x(:,t1).*x(:,t2));
        tau_hash(abs(tau)+1) = tau_hash(abs(tau)+1) +1;
        if(t1<=(frame_length/2))
            tau=tau+1;
        else
            tau = tau-1;
        end
        if((t1+tau)>frame_length || (t1+tau) <1)
            break;
        end
        
    end
end
        sum_instants = sum(y);
        for i=1:frame_length
        z(i) = sum_instants(i)/tau_hash(i);
        end
end
function [y] = time_mean(x,realization_count)
y = zeros(realization_count,1);
for r=1:realization_count
y(r)=average(x(r,:));
end
end
function [z] = time_acf(x,realization_count)
y=zeros(realization_count,realization_count);
z = zeros(1,realization_count);   
tau_hash = zeros(1,realization_count);
for n1 = 1 : realization_count
  tau=0;
    while(1)
        n2=n1+tau;
        y(n1,abs(tau)+1) = average(x(n1,:).*x(n2,:));
        tau_hash(abs(tau)+1) = tau_hash(abs(tau)+1) +1;
        if(n1<=(realization_count/2))
            tau=tau+1;
        else
            tau = tau-1;
        end
        
        if((n1+tau)>realization_count || (n1+tau) <1)
            break;
        end
        
    end
end
    sum_instants = sum(y);
    for i=1:realization_count
        z(i) = sum_instants(i)/tau_hash(i);
    end
end
function z = psd_zero_cross(x,f_domain,linecode,threshold) % get the first zero crossing  
%The threshold sets the theoritical zero cross of the psd to scope it without
%being affected by the error due to the practical data set
error = threshold/12.5; % 8% acceptable error in the zero cross 
switch linecode
case {'Unipolar NRZ','Polar NRZ'}
       l = f_domain(round(x)==0 & abs(f_domain-threshold)<=error);
       z = l(abs((l-threshold))==min(sort(abs(l-threshold))));
case 'Polar RZ'
       l = f_domain(round(x)==0 & abs(f_domain-2*threshold)<=error);
       z = l(abs((l-2*threshold))==min(sort(abs(l-2*threshold))));
end
end