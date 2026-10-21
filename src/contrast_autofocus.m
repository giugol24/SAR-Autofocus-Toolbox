function [focused_opt, beta_opt] = contrast_autofocus(SAR_rcmc, r, lambda, Xa, H_az_freq)
%CONTRAST_AUTOFOCUS Two-stage (Coarse + Fine) contrast optimization
% Metric: contrast = std(intensity)/mean(intensity)

Na = length(Xa);
asse_az = (0:Na-1) - floor(Na/2);
s_rcmc_time = SAR_rcmc;

% Fase 1: Coarse
beta_coarse = linspace(-4e-5, 4e-5, 20);
metric_c = zeros(size(beta_coarse));

for i=1:length(beta_coarse)
    phase = exp(-1j*beta_coarse(i)*(asse_az.^2));
    img = apply_beta(s_rcmc_time, phase, H_az_freq);
    I = abs(img).^2;
    metric_c(i) = std(I(:))/mean(I(:));
end
[~, idx] = max(metric_c);
beta_best = beta_coarse(idx);

% Fase 2: Fine
delta = (8e-5)/20;
beta_fine = linspace(beta_best-delta, beta_best+delta, 20);
metric_f = zeros(size(beta_fine));

for i=1:length(beta_fine)
    phase = exp(-1j*beta_fine(i)*(asse_az.^2));
    img = apply_beta(s_rcmc_time, phase, H_az_freq);
    I = abs(img).^2;
    metric_f(i) = std(I(:))/mean(I(:));
end
[~, idx] = max(metric_f);
beta_opt = beta_fine(idx);

% Final focus
phase_opt = exp(-1j*beta_opt*(asse_az.^2));
focused_opt = apply_beta(s_rcmc_time, phase_opt, H_az_freq);
end

function img = apply_beta(data, phase_vec, H_az)
    data_corr = data .* phase_vec;
    S = fftshift(fft(data_corr,[],2),2);
    img = fftshift(ifft(ifftshift(S.*conj(H_az),2),[],2),2);
end
