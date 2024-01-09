fprintf("\nStart calculating...\n");
for SNR = -5:0
    c = computeC(SNR);
    fprintf("SNR = %d, c = %d\n", SNR, c);
end

function c = computeC(SNR)
    c = 20 / (10 ^ (SNR/20));
end