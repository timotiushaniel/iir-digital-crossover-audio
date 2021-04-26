%% A digital crossover can be designed as shown in Figure 8.60.
% Given audio specifications as:
%   Sampling rate: 44,100 Hz
%   Crossover frequency: 1,000 Hz
%   Highpass filter: third-order Butterworth type at a cutoff frequency of 1,000 Hz
%   Lowpass filter: third-order Butterworth type at a cutoff frequency of 1,000 Hz,
% Use the MATLAB BLT design method to determine: 

%% INITIAL VALUE OF GIVEN AUDIO SPECIFICATIONS
format long
fs = 44100;
order = 3; % Butterworth Filter Order
T = 1/fs; %Sample Period
fc = 1000; % Cut-off frequency
Wn = fc/(fs/2); % Normalized Cut-off Frequency

%% READ ORIGINAL AUDIO
[raw_aud] = audioread('paul.aac');
figure(1);

subplot(2,1,1);
t = (0:T:(length(raw_aud)-1)/fs);
plot(t, raw_aud);grid;
title('Time Domain Representation - Unfiltered Soundl');
xlabel('Time (seconds)');
ylabel('Amplitude');
xlim([0 t(end)]);

% Fourier Transform
aud_N = length(raw_aud); %Panjang sinyal dari raw_aud
%POW
% Mencari panjang sample asli.
% n adalah panjang vektor hasil perhitungan FFT. 
% idealnya adalah bilangan dari angka 2 yang dipangkatkan, 
% dan lebih panjang dari sinyal informasi (m). Sehingga digunakanlah fungsi nextpow2.

n = pow2(nextpow2(aud_N)); %nextpow2(m) = 17 >> pow2(17) >> 2^17 = 131072
% Mengubah panjang sinyal, jadi jumlah sample = pow2
% Membuat perhitungan FFT menjadi lebih cepat
% Terutama untuk ukuran sample data dengan faktor prima yang lebih besar

% Untuk mendapatkan normalisasi magnitud sinyal, 
% nilai hasil FFT harus dibagi dengan L (panjang data)
y = fft(raw_aud, n);
f = (0:n-1)*(fs/n); %Membuat sumbu x dalam domain frekuensi (hasil ft)
amplitude = abs(y)/n; %Meng-absolute kan nilai dari y (fft dari sample_data)
subplot(2,1,2)
%semilogy(f(1:floor(n/2)),amplitude(1:floor(n/2)))
plot(f(1:floor(n/2)),amplitude(1:floor(n/2))),grid;
title('Frequency Domain Representation - Unfiltered Sound')
xlabel('Frequency')
% ylabel('dB')
ylabel('Amplitude')

%% LOW PASS FILTER OF BUTTERWORTH
[b_low,a_low] = butter(order,Wn,'low');

% Plot the magnitude and phase responses
figure(2);
freqz(b_low,a_low,512,fs);
axis([0 fs/2 -40 5]);

%% HIGHPASS FILTER OF BUTTERWORTH
format long
order = 3;
fc = 1000;
Wn = fc/(fs/2);
[b_high,a_high] = butter(order,Wn,'high');

% Plot the magnitude and phase responses
figure(3);
freqz(b_high,a_high,512,fs);
axis([0 fs/2 -40 5]);

%% LOWPASS FILTERED AUDIO
% Add filter to audio
lowpass_filter = filter(b_low,a_low,raw_aud);
% Fourier Transform
aud_N = length(lowpass_filter); 
aud_f = [0:aud_N/2]*fs/aud_N; 
Axk = 2*abs(fft(lowpass_filter))/aud_N;
Axk(1) = Axk(1)/2;

figure(4);
subplot(2,1,1);
plot(aud_f,Axk(1:aud_N/2+1)); 
axis([0 22050 0 0.2]);
title('Lowpass Filtered Audio Spectrum');
xlabel('Frequency (Hz)');ylabel('Amplitude');grid;

%% HIGHPASS FILTERED AUDIO
% Add filter to audio
highpass_filter = filter(b_high,a_high,raw_aud);
% Fourier Transform
aud_N = length(highpass_filter); 
aud_f = [0:aud_N/2]*fs/aud_N; 
Axk = 2*abs(fft(highpass_filter))/aud_N;
Axk(1) = Axk(1)/2;

subplot(2,1,2);
plot(aud_f,Axk(1:aud_N/2+1));
axis([0 22050 0 0.2]);
title('Highpass Filtered Audio Spectrum');
xlabel('Frequency (Hz)');ylabel('Amplitude');grid;