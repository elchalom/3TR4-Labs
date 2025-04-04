clear
hold off
format long e
N = 4096; %No. of FFT samples
sampling_rate = 500e3; %unit Hz
tstep = 1/sampling_rate;
tmax = N*tstep/2;

tmin = -tmax;
tt = tmin:tstep:tmax-tstep;
fmax = sampling_rate/2; 
fmin = -fmax;
fstep = (fmax-fmin)/N;
freq = fmin:fstep:fmax-fstep;

% Signal Parameter
% Carrier
% i <3 ty
fc=64e3;
Ac = 1;
ct=Ac*cos(2*pi*fc*tt);

% Message
Am = 1;
fm = 4e3; 
mt = Am*cos(2*pi*fm*tt);

% modulated signal
st = mt.*ct;

% new carrier power is 0.125 or 0.375 or 0.500 // total sideband power(2*(Ac^2)*Am^2)/8 = 0.25 //half of total sideband power = ((Ac^2)*Am^2)/8 = 0.125
% requires the newly added carrier to have amplitude of 1/2
% AM Signal

% amplitude required is 0.5, 0.866, and 1
new_carrier_amplitude = 1;
% 100% modulation
modulation = 1;

AM_signal = new_carrier_amplitude.* ct + modulation .* st;


%%
figure(1)
Hp1 = plot(tt,ct);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('Carrier c(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Carrier : Time domain');
axis([-1e-3 1e-3 -1.1 1.1])
pause(1)

%%
figure(2)
Hp1 = plot(tt,mt);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('message  m(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('message signal : Time domain');
axis([-0.01 0.01 0 1.1])
pause(1)
%%
% DSB-SC in time domain

% change st to AM_signal for the am signal

figure(3)

Hp1 = plot(tt,st);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('s(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('DSB-SC modulated wave : Time domain');
axis([-10/fc 10/fc min(st) max(st)])
pause(1)

%%
figure(4)
%The amplitude of the spectrum is different from the Fourier transform
%amplitude due to discretization of discrete Fourier transform
Mf = fftshift(fft(fftshift(mt)))/(2*fmax);
Hp1=plot(freq,abs(Mf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|M(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the message signal');
axis ([-15e3 15e3 0 max(abs(Mf))])
pause(1)

%%
% DSB-SC in frequency domain

% change st to AM_Signal for AM
Sf = fftshift(fft(fftshift(st)))/(2*fmax);
figure(5)
Hp1=plot(freq,abs(Sf));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum of the DSB-SC wave S(f)');
axis([-10e4 10e4 0 max(abs(Sf))])

xticks(-10e4:1e4:10e4)

%%
pause(1)

%DSB-SC demodulation
%Local oscillator at the receiver perfectly synchronized
thet=0;
lo = cos(2*pi*fc*tt + thet); 
st1 = st .* lo;
figure(6)
Hp1=plot(tt,st1);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('  s hat(t)  (Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('signal after remodulation at Rx, s hat(t) : Time domain');
axis([-10/fc 10/fc min(st1) max(st1)])
pause(1)

%%
Sf1 = fftshift(fft(fftshift(st1)))/(2*fmax);
figure(7)
Hp1=plot(freq,abs(Sf1));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|S hat(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum S hat(f)');
axis ([-200e3 200e3 0 max(abs(Sf1))])
pause(1)

%%
%Low pass filtering
f_cutoff = 30e3;
%ideal low pass filter
Hf = zeros(1, length(freq));
n=1;
for f = freq
    if abs(f) < f_cutoff
        Hf(n) = 1;
    else
        Hf(n) = 0;
    end
n=n+1;
end

Mf1 = Sf1 .* Hf;
figure(8)
Hp1=plot(freq,abs(Mf1));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Frequency (Hz) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('|M(f)|');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Spectrum output of low pass filter signal');
axis ([-1e4 1e4 0 max(abs(Mf1))])
pause(1)

%%
mt1 = 2*fmax*fftshift(ifft(fftshift(Mf1)));
figure(9)
Hp1=plot(tt,mt1,'r',tt,mt*0.5,'g.');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
Hx=xlabel('Time (sec) ');
set(Hx,'FontWeight','bold','Fontsize',16)
Hx=ylabel('m hat(t)(Volt)');
set(Hx,'FontWeight','bold','Fontsize',16)
title('Output of low pass filter,  m hat(t)  : Time domain');
axis([-0.01e-1 0.01e-1 min(mt*0.5) max(mt*0.5)])
legend('LPF output', 'message sig');
pause(1)





