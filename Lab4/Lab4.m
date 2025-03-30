% This code sets up the time and frequency vectors for all the numerical
% experiments of Lab3

% Junseo Mun and Matthew El Chalouhi

clear 
endTime = 10;
beginTime = -10;
N=100000;
timeStep = (endTime-beginTime)/N;
sampling_rate = 1/timeStep;

% Time window =
timeWindow = beginTime:timeStep:endTime-timeStep;

% Data 
%load('lab4_num_expt1')
load('lab4_num_expt2')

maxLag = 20000;
% Autocorrelation of yt
Ry = xcorr(yt,yt,maxLag);
% Vector of Tau values
tauVector = -(maxLag*timeStep):timeStep:maxLag*timeStep;
% Abs. PSD corresponding to yt
Sy = abs(fftshift(fft(fftshift(Ry))));
% define the frequency vector corresponding to tauVector
Ntau = length(tauVector);
% Nyquist sampling rate
fMax = sampling_rate/2; 
fMin = -fMax;
fStep = (fMax-fMin)/Ntau;
% Frequency window
freq = fMin:fStep:fMax-fStep;

%%

% Plot Autocorrelation Function
figure(1)
Hp1 = plot(tauVector, Ry);
set(Hp1, 'LineWidth', 1)
Ha = gca;
set(Ha, 'FontSize', 16)
Hx = xlabel('\tau (s)');
set(Hx, 'FontWeight', 'bold', 'FontSize', 16)
Hx = ylabel('Ry');
Hx = yline(0);
set(Hx, 'FontWeight', 'bold', 'FontSize', 16)
title('Autocorrelation Plot');
axis([min(tauVector) max(tauVector) min(Ry) max(Ry)])

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

%%

% Plot of our time domain signal
figure(3);
plot(timeWindow,yt);
xlim([-100*timeStep 100*timeStep]);
title ("Time Domain Signal");
subtitle("y(t)", 'interpreter', 'latex');
xlabel("Time(s)", 'FontWeight', 'bold');
ylabel("y(t)", "FontWeight", 'bold');


%%
Nyt = length(yt);

fMax = sampling_rate/2;
fMin = -fMax;
fStep = (fMax-fMin)/Nyt;

freq = fMin:fStep:fMax-fStep;

% Plot Magnitude Spectrum
figure(4);
plot(freq, abs(fftshift(fft(fftshift(yt)))));
title("Magnitude Plot");
xlabel("Frequency (Hz)", 'FontWeight', 'bold');
ylabel("Magnitude", 'FontWeight', 'bold');
grid on;