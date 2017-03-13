%% DemodulateSignal

% Copyright (c) 2017, Felix Schlegel
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


% Section 1: Data Import

clearvars data
data = lvm_import() % available on: mathworks.com/matlabcentral/fileexchange/19913
modF=input('Modulation Frequency in Hz? ');


% check for nans in ch1
clearvars ch1 ch2 t_raw y k y1 y21
ch1 = data.Segment1.data(:,2);                                                
ch2 = data.Segment1.data(:,4);
t_raw = data.Segment1.data(:,1);

y=[];
k=1;
for l = 24:length(ch1)
    if isnan(ch1(l))
    else
        y1(k)= (-1)*ch1(l);
        y2(k)= ch2(l);  
        t(k)= t_raw(l);
        k = k+1;
    end
end

% digital lock-in for the reference channel

clearvars NFFT Y f YSig_sub demodSig
Fs=1/data.Segment1.Delta_X(1);  %readout sampling rate
L = length(y1);   % Length of signal
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(y1,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% power spectral density
figure;
hold on
[psdestx,Fxx] = periodogram(y1,rectwin(length(y1)),length(y1),Fs);
plot(Fxx,10*log10(psdestx)); grid on;
xlabel('Hz'); ylabel('Power/Frequency (dB/Hz)');
title('Periodogram Power Spectral Density Estimate');

%take the power in a frequency band
k = 1;
powbandRef=[];
demodSig=[];
N_samples = 50e-3*Fs;

t_downs=[];
t_downs=downsample(t_raw,N_samples);
t_downs=t_downs(1:end-1)-min(t_downs);

for l = 1:(length(y2)-N_samples)
    ySig_sub = y1(k:k + N_samples);
    demodSig(l) = bandpower(ySig_sub,Fs,[modF-0.5 modF+0.5]);
%   demodSig(l) = bandpower(ySig_sub,Fs,[1099.5 1100.5]); 
    k= k + N_samples;
end


%% Section 2: Detrending and plotting

% Simple detrending using a polynomial fit:
demodSig_detrended=[];
opol = 2;
[p,s,mu] = polyfit(t_downs',demodSig,opol);
f_y = polyval(p,t_downs,[],mu);
filtered_signal=demodSig-f_y'+mean(demodSig);


% Find stimulus events:
z=diff(ch2);
k=1;
stim_events=[];
for i=2:length(z)
    if z(i)>1 && z(i-1)<1
    stim_events(k)=i;
    k=k+1;
    end
end
stim_events=stim_events/(20*N_samples);


% Plot signals before/after detrending, add stimulus events.
figure('Color','w','Renderer','painters');
hold on
plot(t_downs,demodSig,'k');
plot(t_downs,filtered_signal); 

for i=1:length(stim_events)
vline(stim_events(i),'r')  % The function vline is available on: mathworks.com/matlabcentral/fileexchange/1039
end


%% Section 3: Calcium Regressor

cutoff_time=20; % initial time in s discarded in the fMRI analysis (to reach magnetization equilibrium)
bin_size=20; % By default 20 (=sampling rate of calcium signal) for 1s TR.

filtered_signal_cropped=filtered_signal((stim_events(1)+cutoff_time)*sr:end);  % first stim_event = start of EPI scan (in case of resting state)
t_downs_cropped=t_downs((stim_events(1)+cutoff_time)*sr:end); % vice versa for time axis
t_downs_cropped=t_downs_cropped-min(t_downs_cropped)+1;

% smoothing/downsampling of the data
calcium_regressor=downsamplebin(filtered_signal_cropped,bin_size); %averaging the data in bins of 1s (EPI temporal resolution)
% scaling data to prevent rounding errors (values between 0 and 1)
calcium_regressor=calcium_regressor-min(calcium_regressor);
scaling_factor=1/max(calcium_regressor);
calcium_regressor=calcium_regressor.*scaling_factor; 


figure
subplot(2,1,1)
plot(t_downs_cropped,filtered_signal_cropped); 
subplot(2,1,2)
plot(calcium_regressor)

%write .1D files for AFNI
dlmwrite('C:\calcium_regressor.1D',calcium_regressor,'delimiter','\t');
calcium_regressor_time_axis=0:length(calcium_regressor)-1;
dlmwrite('C:\calcium_regressor_time_axis.1D',calcium_regressor_time_axis,'delimiter','\t');

