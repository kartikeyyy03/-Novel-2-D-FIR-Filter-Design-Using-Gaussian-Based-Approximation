% MATLAB code for Parks-McClellan 1D FIR filter extended to 2D using freqz2

% Filter specifications
N = 11; % Filter size (11x11)
Fp = 0.3; % Passband edge (normalized frequency, 0 to 1)
Fst = 0.36; % Stopband edge (normalized frequency, 0 to 1)
Rp = 1; % Passband ripple in dB
Rs = 60; % Stopband attenuation in dB

% Update weight function for Parks–McClellan method
weights = [0.7 1]; % 0.7 for passband, 1 for stopband

% Design the Parks-McClellan FIR filter (1D)
filter_1d = firpm(N-1, [0 Fp Fst 1], [1 1 0 0], weights); % Apply weights
filter_2d = filter_1d' * filter_1d; % Outer product to get 2D filter
% Design the Parks-McClellan FIR filter (1D)

% Compute the frequency response
[H, w1, w2] = freqz2(filter_2d, 256, 256); % 2D frequency response
H_dB = 20 * log10(abs(H) + eps); % Convert magnitude to dB (add eps for log safety)

% Plot the frequency response
figure;
surf(w1/pi, w2/pi, H_dB);
shading interp;
xlabel('Normalized frequency (x)');
ylabel('Normalized frequency (y)');
zlabel('Magnitude response (in dB)');
title('Parks-McClellan Frequency Response (11x11)');
axis([-1 1 -1 1 -60 5]); % Adjust axis limits
colormap jet;
colorbar;