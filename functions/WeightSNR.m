function W_snr = WeightSNR(matrix_snr,MaxSnr)

matrix_snr(abs(matrix_snr) < 30) = 0;
matrix_snr = matrix_snr./MaxSnr;

W_snr = diag(matrix_snr(:,1)) .* diag(matrix_snr(:,1));
end