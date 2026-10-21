function SAR_rcmc = rcmc_correction(SAR_signal, r, Xa)
%RCMC_CORRECTION Range Cell Migration Correction via interpolation
[Nr, Na] = size(SAR_signal);
SAR_rcmc = zeros(Nr, Na);

for ia = 1:Na
    delta_r = (Xa(ia)^2) ./ (2 * r);
    r_migrato = r + delta_r;
    SAR_rcmc(:, ia) = interp1(r, SAR_signal(:, ia), r_migrato, 'linear', 0);
end
end
