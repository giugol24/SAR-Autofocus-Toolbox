function I_f = frost_filter_fast(I, W, K)
% Versione vettorizzata del Frost - niente doppi for
h = ones(W)/W^2;
mu = imfilter(I, h, 'symmetric');
var_loc = imfilter(I.^2, h, 'symmetric') - mu.^2;
% ... pesi esponenziali applicati con conv2
I_f = ... % 5 righe invece di 2 for annidati
end
