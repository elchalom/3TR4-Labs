%This code sets up the time and frequency vectors for all the numerical
%experiments of Lab3

%Junseo Mun and Matthew El Chalouhi
clear 
format long e
tend = 10;
tbeg = -10;
N=100000;
tstep = (tend-tbeg)/N;
sampling_rate = 1/tstep;

%Time window =
tt = tbeg:tstep:tend-tstep;

% load('lab4_num_expt1')
load('lab4_num_expt2')
%load('lab4_num_expt3')

maxlag = 20000;
%Autocorrelation of yt
Ry  = xcorr(yt,yt,maxlag);
%tau vector
tau_vec = -(maxlag*tstep):tstep:maxlag*tstep;
%Abs. PSD corresponding to yt
Sy = abs(fftshift(fft(fftshift(Ry))));
%define the frequency vector corresponding to tau_vec
Ntau = length(tau_vec);
%Nyquist sampling rate
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/Ntau;
%Frequency window
freq = fmin:fstep:fmax-fstep;

%%
%Autocorrelation Function

figure(1)
Hp1 = plot(tau_vec, Ry);
set(Hp1, 'LineWidth', 1)
Ha = gca;
set(Ha, 'FontSize', 16)
Hx = xlabel('\tau (s)');
set(Hx, 'FontWeight', 'bold', 'FontSize', 16)
Hx = ylabel('Ry');
Hx = yline(0);
set(Hx, 'FontWeight', 'bold', 'FontSize', 16)
title('Autocorrelation Plot');
axis([min(tau_vec) max(tau_vec) min(Ry) max(Ry)])
pause(1)

%%
% Plot of power spectrum density
figure(2)
Hp1 = plot(freq, Sy);
set(Hp1, 'LineWidth', 1)
Ha = gca;
set(Ha, 'FontSize', 16)
Hx = xlabel('Frequency (Hz)');
set(Hx, 'FontWeight', 'bold', 'FontSize', 16)
Hx = ylabel('Sy');
Hx = yline(0);
set(Hx, 'FontWeight', 'bold', 'FontSize', 16)
title('Power Spectrum Density Plot');
axis([min(freq) max(freq) min(Sy) max(Sy)])
pause(1)

%%
fig = figure(3);
plot(tt,yt);
xlim([-100*tstep 100*tstep]);
title ("Time Domain Signal");
subtitle("y(t)", 'interpreter', 'latex');
xlabel("Time(s)", 'FontWeight', 'bold');
ylabel("y(t)", "FontWeight", 'bold');


%%
fig = figure(4);
Nyt = length(yt);

fmax = sampling_rate/2;
fmin = -fmax;
fstep = (fmax-fmin)/Nyt;

freq = fmin:fstep:fmax-fstep;

% Plot Magnitude Spectrum

plot(freq, abs(fftshift(fft(fftshift(yt)))));
title("Magnitude Plot");
xlabel("Frequency (Hz)", 'FontWeight', 'bold');
ylabel("Magnitude", 'FontWeight', 'bold');
grid on;