%% Square wave generator
clc
clear all
hold off

f0=500;     %fundamental freq of input sqaure wave
T0 = 1/f0;  %period 
tstep = 0.005*T0;
no_sample = 3*T0/tstep + 1; %no. of samples  within  3*T0

%tt = -0.5*T0:tstep:0.5*T0;
tt = -1.5*T0:tstep:1.5*T0;

square_in = 2*square(tt*2*pi*f0,50);

figure(1)
Hp1 = plot(tt,square_in);
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('input - time domain')
pause

%% Fourier series representation of signal (Amplitude Spectrum)
      
K=1/(2*pi);
N= 100; %no. of harmonics
nvec = -N:N;
c_in = zeros(size(nvec)); %fourier coefficients
A = 2;
for n = nvec
    m = n+N+1;
    c_in(m) = (2*A / (n * pi)) * sin((n * pi) / 2);
    
    if (n == 0)
      c_in(m) = 0.0;
    end
end


f = nvec*f0; %frequency vector
figure(2)
Hp1=stem(f,abs(c_in));
axis([-8*f0 8*f0 0 max(abs(c_in))])
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('magnitude spectrum of input')
pause

%% Fourier series representation of signal (Phase Spectrum)

figure(3)
Hp1=stem(f,angle(c_in));
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
axis([-0.1e5 0.1e5 -pi pi])
title('phase spectrum of input')
pause

%% Designing the 2nd order Butterworth filter


R=10e3; %10k ohms
C=1e-8; %10 nF
fc=1/(2*pi*R*C);     %cutoff freq of filter
%fc = 5000;

Hf = 1./ (power(1i*f/fc,2) + 1.414*(1i*f/fc) + 1); %filter transfer function

c_out = c_in .* Hf; %Fourier coefficients of the filter output

figure(4)
stem(f,abs(c_in),'r','LineWidth',2);
hold on
stem(f,abs(c_out),'b','LineWidth',2);
hold off
axis([-8*f0 8*f0 0 max(abs(c_in))])
Ha = gca;
set(Ha,'Fontsize',16)
title('magnitude spectrum of filter output and input')
Ha = gca;
set(Ha,'Fontsize',16)
legend('input','output')
% pause

% figure(5)
% Hp1=plot(f,angle(c_out));
% set(Hp1,'LineWidth',2)
% Ha = gca;
% set(Ha,'Fontsize',16)
% title('phase spectrum of output')
% axis([-0.1e4 0.1e4 -pi pi])

figure(5)
stem(f,angle(c_in),'r','LineWidth',2);
hold on
stem(f,angle(c_out),'b','LineWidth',2);
hold off
axis([-0.1e4 0.1e4 -pi pi])
Ha = gca;
set(Ha,'Fontsize',16)
title('phase spectrum of input and output')
Ha = gca;
set(Ha,'Fontsize',16)
legend('input','output')
pause


%% Construct the output signal from the Cout Fourier coefficients

A = zeros(3*N+1,ceil(no_sample));
for n = nvec
    m=n+N+1;
    A(m,:) = c_out(m) .* exp(1i*2*pi*n*f0*tt);
end
gp_out = sum(A);
figure(6)
Hp1 = plot(tt,real(gp_out),'b',tt,square_in,'r');
set(Hp1,'LineWidth',2)
Ha = gca;
set(Ha,'Fontsize',16)
title('filter input and output-time domain')
set(Ha,'Fontsize',16)
legend('output','input')
pause