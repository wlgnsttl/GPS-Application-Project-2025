function W_snr = WeightSNR(matrix_snr)

matrix_snr(abs(matrix_snr) < 30) = 0;
matrix_snr = matrix_snr./100;

W_snr = diag(matrix_snr(:,1));