function [focused_img, phase_error] = pga_autofocus(SAR_signal, r, lambda, Xa, varargin)
%PGA_AUTOFOCUS Phase Gradient Autofocus with LUMV estimator
%   Handles both clean and noisy data via adaptive windowing and thresholding.

p = inputParser; addParameter(p,'W',64); addParameter(p,'maxIter',20);
parse(p,varargin{:});

Nr = size(SAR_signal,1); Na = length(Xa);
phase_error = zeros(1, Na);
focused_img = azimuth_focus(SAR_signal, r, lambda, Xa(2)-Xa(1));

W = p.Results.W;

for iter = 1:p.Results.maxIter
    % Thresholding - keep only strong targets (crucial for noisy case)
    energy = max(abs(focused_img),[],2);
    idx_valid = energy >= prctile(energy, 85);
    if sum(idx_valid) < 3, idx_valid = true(Nr,1); end

    g = zeros(sum(idx_valid), Na);
    cnt = 1;
    for ir = find(idx_valid)'
        [~, pmax] = max(abs(focused_img(ir,:)));
        r_shifted = circshift(focused_img(ir,:), floor(Na/2)+1 - pmax);
        mask = zeros(1,Na); mask(max(1,floor(Na/2)-W/2):min(Na,floor(Na/2)+W/2)) = 1;
        g(cnt,:) = ifft(ifftshift(r_shifted .* mask));
        cnt = cnt+1;
    end

    g_dot = gradient(g,2);
    phi_dot = sum(imag(conj(g).*g_dot),1) ./ (sum(abs(g).^2,1)+eps);
    phi_e = detrend(cumsum(phi_dot),1);

    if std(phi_e) <= pi/4, break; end

    raw = ifft(ifftshift(focused_img,2),[],2);
    focused_img = fftshift(fft(raw .* exp(-1j*phi_e),[],2),2);
    phase_error = phase_error + phi_e;
    W = max(floor(W*0.8), 8);
end
end
