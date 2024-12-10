%%A Novel 2-D FIR Filter Design Methodology Based on a Gaussian-Based
%%Approximation:
% Parameters (from research paper)
wc = 0.3; % Cut-off frequency
ws = 0.36; % Stop-band edge
delta = 1; % Maximum ripple in the pass-band
N = 11; % Filter size (11x11)

% Derived parameters(From research paper)
sigma = 0.8 * (ws - wc); % Transition bandwidth scaling
Delta_mu = 2 * sqrt(-2 * sigma^2 * log(1 - delta)); % Equispacing of Gaussian centers
N1 = ceil(wc / Delta_mu); % Bound for Gaussian centers

% Gaussian center coordinates (mu_x and mu_y)
%mu_x = (-N1:N1) * Delta_mu; % Centers along x-axis
%mu_y = (-N1:N1) * Delta_mu; % Centers along y-axis

% Compute M
M = 1 / (2 * pi * sigma^2) * (1 + ...
    2 * sum((1 - delta).^(4 * (1:N1).^2)) + ...
    2 * sum((1 - delta).^(4 * (1:N1).^2)) + ...
    4 * sum(((1-delta).^(4*(1:N1).^2)).*((1-delta).^(4*(1:N1).^2))));


% Designing the filter by computing the filter impulse response
h = zeros(N, N);
mid = (N - 1) / 2;
for n = -mid:mid
    for m = -mid:mid
        gaussian_product = exp(-(n^2 + m^2) * sigma^4 / 2);
        cosine_sum = 1 + 2 * sum(cos(n * (1:N1)*Delta_mu)) + ...
                        2 * sum(cos(m * (1:N1)*Delta_mu)) + ...
                        4 * sum(cos(n * (1:N1)*Delta_mu) .* cos(m * (1:N1)*Delta_mu));
        h(n + mid + 1, m + mid + 1) = gaussian_product * cosine_sum / (4 * pi^2 * M);
    end
end

% Normalizing the filter so that the total energy of the filter remains
% constant
%h(:) reshapes the @D filter into a single column vector.
%sum(h(:)) computes the sum of all elements in the filter
h = h / sum(h(:));
%h/sum(h(:)) ensures that the sum of all filter coefficients equals to 1.

% Frequency response in dB
[H, w1, w2] = freqz2(h, 256, 256); % computes 2D frequency response. 256 specifies the resolution(size) of the response.
% This means the frequency response is sampled on a 256*256 grid in the normalized frequency domain.
H_dB = 20 * log10(abs(H)); % Convert magnitude to dB

% Plotting the frequency response
figure;
surf(w1/pi, w2/pi, H_dB);
shading interp;
xlabel('Normalized frequency (x)');
ylabel('Normalized frequency (y)');
zlabel('Magnitude response (in dB)');
title('Frequency Response for 11x11 2D FIR Filter');
axis([-1 1 -1 1 -35 5]); % Match graph limits
colormap jet;
colorbar;

%% 


image=imread('library.jpg');
image_gray=rgb2gray(image);


filtered_image=conv2(double(image_gray),h,'same');


imshow(uint8(filtered_image));
title('Filtered Image');


high_pass_image=double(image_gray) - filtered_image;
imshow(uint8(high_pass_image));
title('High-pass Filtered Image');


Sharpened_image = double(image_gray)+high_pass_image;


figure;


subplot(2,2,1);
imshow(image_gray);
title('Original Image');


subplot(2,2,2);
imshow(uint8(filtered_image));
title('Filtered Image');


subplot(2,2,3);
imshow(uint8(high_pass_image));
title('High-pass Filtered Image');


subplot(2,2,4);
imshow(uint8(Sharpened_image));
title('Sharpened Image');
%% 


