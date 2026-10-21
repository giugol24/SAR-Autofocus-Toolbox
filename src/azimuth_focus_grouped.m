function focused_group = azimuth_focus_grouped(SAR_signal, r, lambda, dx, ngroups)
%AZIMUTH_FOCUS_GROUPED Focus using shared matched filter per range block
%   Used to demonstrate Depth of Focus / defocusing effect.

[Nr, Na] = size(SAR_signal);
group_idx = ceil((1:Nr) / (Nr/ngroups));
focused_group = zeros(Nr, Na);
fa = (-Na/2 : Na/2-1) / (Na * dx);

for g = 1:ngroups
    ind = group_idx == g;
    Rmean = mean(r(ind));
    Ka = 2 / (lambda * Rmean);
    Hm = exp(1j * pi * (fa.^2) / Ka);
    for ir = find(ind)
        Sf = fftshift(fft(SAR_signal(ir,:)));
        focused_group(ir,:) = ifft(ifftshift(Sf .* conj(Hm)));
    end
end
end
