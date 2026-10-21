function ENL = estimate_ENL(I_intensity, r1, r2, c1, c2)
roi = I_intensity(r1:r2, c1:c2);
ENL = mean(roi(:))^2 / var(roi(:));
end
