%% A digital crossover can be designed as shown in Figure 8.60.
% Given audio specifications as:
%   Sampling rate: 44,100 Hz
%   Crossover frequency: 1,000 Hz
%   Highpass filter: third-order Butterworth type at a cutoff frequency of 1,000 Hz
%   Lowpass filter: third-order Butterworth type at a cutoff frequency of 1,000 Hz,
% use the MATLAB BLT design method to determine: 

clear all;
fs=44100; % Sampling Rate
T=1/fs; % Periods
f=1000; % Crossover Frequency
wd=2*pi*f; % Digital frquency in radians per second (Wd=2*pi*f (rads/s))
wa=(2/T)*tan(wd*T/2); % Design Procedure (rads/s)

%% A) Frequency responses for the highpass filter and the lowpass filter

% Lowpass Filter
figure(1)
Bp=1; % Bp: vector containing the numerator coefficients of the lowpass prototype
Ap=[1 2 2 1]; % Ap: vector containing the denominator coefficients of the lowpass prototype -> (S+1)*(S^2+S+1) = S^3 + 2S^2 + 2s + 1 -> Butterworth Order 3

% B_l: vector containing the numerator coefficients of the analog filter.
% A_l : vector containing the denominator coefficients of the analog filter
[B_l,A_l]=lp2lp(Bp,Ap,wa); % Lowpass to lowpass

% bL: vector containing the numerator coefficients of the digital filter for the Lowpass.
% aL: vector containing the denominator coefficients of the digital filter for the Lowpass.
[bL,aL]=bilinear(B_l,A_l,fs); % Plot of the magnitude and phase frequency responses of the digital filter
[hL,ff]=freqz(bL,aL,512,fs);

freqz(bL,aL,512,fs);

% Highpass Filter
figure(2)
Bp=1; % Bp: vector containing the numerator coefficients of the highpass prototype
Ap=[1 2 2 1]; % Ap: vector containing the denominator coefficients of the lowpass prototype -> (S+1)*(S^2+S+1) = S^3 + 2S^2 + 2s + 1

% B_h: vector containing the numerator coefficients of the analog filter.
% A_h : vector containing the denominator coefficients of the analog filter
[B_h,A_h]=lp2hp(Bp,Ap,wa); % Lowpass to highpass

% bH: vector containing the numerator coefficients of the digital filter for the Highpass.
% aH: vector containing the denominator coefficients of the digital filter for the Highpass.
[bH,aH]=bilinear(B_h,A_h,fs) % Plot of the magnitude and phase frequency responses of the digital filter

[hH,ff]=freqz(bH,aH,512,fs);
freqz(bH,aH,512,fs)

%% B) combined frequency response for both filters.
figure(3)
H=abs(hL)+abs(hH);
plot(ff,20*log10(abs(hL)),ff,20*log10(abs(hH)),'-.', ff,20*log10(H)); 
grid on;

%% C) The Transfer Function and Difference Equations for the Lowpass and
%     Highpass

% Lowpass:
% 1) Transfer Function:
tf_l=filt(bL,aL);
tf_l

% 2) Difference Equation:
de_l=idpoly([1 -2.715 2.47 -0.7519],[0.0003151 0.0009452 0.0009452 0.0003151],'NoiseVariance',0)
de_l

% Highpass:
% 1) Transfer Function:
tf_h=filt(bH,aH);
tf_h

% 2) Difference Equation:
de_h=idpoly([1 -2.715 2.47 -0.7519],[0.8671 -2.601 2.601 -0.8671],'NoiseVariance',0)
de_h