function focused_image = azimuth_focus(SAR_signal, r, lambda, dx)
%AZIMUTH_FOCUS Focus SAR image with exact range-dependent matched filter
%   Performs azimuth compression for stripmap SAR, zero squint.

[Nr, Na] = size(SAR_signal);
focused_image = zeros(Nr, Na);
fa = (-Na/2 : Na/2-1) / (Na * dx); % spatial frequencies [1/m]

for ir = 1:Nr
    Ka = 2 / (lambda * r(ir)); % chirp rate spatial domain
    Hm = exp(1j * pi * (fa.^2) / Ka);
    Sf = fftshift(fft(SAR_signal(ir,:)));
    focused_image(ir,:) = ifft(ifftshift(Sf .* conj(Hm)));
end
end
