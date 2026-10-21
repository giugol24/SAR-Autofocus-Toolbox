function I_wind = hamming_filter(S_over, r0, c0, nr, nc, overs)
w_r = hamming(nr); w_a = hamming(nc);
W_2D = w_r * w_a';
W_over = zeros(size(S_over));
W_over(r0:r0+nr-1, c0:c0+nc-1) = W_2D;
I_wind = overs^2 * ifft2(ifftshift(S_over .* W_over));
end
