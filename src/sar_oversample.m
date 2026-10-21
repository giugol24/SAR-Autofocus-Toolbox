function [I_over, S_over] = sar_oversample(I_sub, overs)
% Zero-padding in frequenza per oversampling
[nr, nc] = size(I_sub);
S = fftshift(fft2(I_sub));
S_over = zeros(nr*overs, nc*overs);
r0 = floor((nr*overs-nr)/2)+1; c0 = floor((nc*overs-nc)/2)+1;
S_over(r0:r0+nr-1, c0:c0+nc-1) = S;
I_over = overs^2 * ifft2(ifftshift(S_over));
end
