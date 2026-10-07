clc; clear; close all;
rng(1);
numTx = 2;
numRx = 2;
f = 2.4e9;
Pt_dBm = 0;
Pt = 10^((Pt_dBm-30)/10);
K = 10;
shadowingStd = 6;
TxCorr = toeplitz([1, 0.7]);
RxCorr = toeplitz([1, 0.5]);
RxCorrS = sqrtm(RxCorr);
TxCorrS = sqrtm(TxCorr);
BW = 10e6;
N0_dBmHz = -174;
NF = 7;
Noise_dBm = N0_dBmHz + 10*log10(BW) + NF;
N0 = 10^((Noise_dBm-30)/10);
numDistances = 50;
d = linspace(10, 500, numDistances);
numTrials = 1000;
R_target = 10;
antennaGain = @(theta) 10*log10(cos(theta).^2 + 0.01);
PL_3GPP = @(dist_m) 28.0 + 22*log10(dist_m) + 20*log10(f/1e9);
Pr_dBm     = zeros(1, numDistances);
Capacity   = zeros(1, numDistances);
OutageProb = zeros(1, numDistances);
theta = 0;
G_tx = 10^(antennaGain(0)/10);
G_rx = 10^(antennaGain(theta)/10);
useParfor = license('test','Distrib_Computing_Toolbox');
for idx = 1:numDistances
    di = d(idx);
    PL_dB  = PL_3GPP(di) + shadowingStd*randn(1, numTrials);
    PL_lin = 10.^(-PL_dB/10);
    Pr = Pt .* G_tx .* G_rx .* PL_lin;
    capacity_trials = zeros(1, numTrials);
    outage_flags    = false(1, numTrials);
    if useParfor
        parfor t = 1:numTrials
            H_LoS  = ones(numRx, numTx);
            H_NLoS = (randn(numRx, numTx) + 1j*randn(numRx, numTx))/sqrt(2);
            H = sqrt(K/(K+1))*H_LoS + sqrt(1/(K+1))*H_NLoS;
            estErrorVar = 0.05;
            H_est = H + sqrt(estErrorVar)*(randn(numRx,numTx)+1j*randn(numRx,numTx))/sqrt(2);
            [~, ~, V] = svd(H_est,'econ');
            f_vec = V(:,1);
            h_eff = H * f_vec;
            SNR = (norm(h_eff)^2 * Pr(t)) / N0;
            C = log2(1 + SNR);
            capacity_trials(t) = C;
            outage_flags(t)    = (C < R_target);
        end
    else
        for t = 1:numTrials
            H_LoS  = ones(numRx, numTx);
            H_NLoS = (randn(numRx, numTx) + 1j*randn(numRx, numTx))/sqrt(2);
            H_NLoS = RxCorrS * H_NLoS * TxCorrS;
            H = sqrt(K/(K+1))*H_LoS + sqrt(1/(K+1))*H_NLoS;
            estErrorVar = 0.05;
            H_est = H + sqrt(estErrorVar)*(randn(numRx,numTx)+1j*randn(numRx,numTx))/sqrt(2);
            [~, ~, V] = svd(H_est,'econ');
            f_vec = V(:,1);
            h_eff = H * f_vec;
            SNR = (norm(h_eff)^2 * Pr(t)) / N0;
            C = log2(1 + SNR);
            capacity_trials(t) = C;
            outage_flags(t)    = (C < R_target);
        end
    end
    Pr_dBm(idx)     = 10*log10(mean(Pr)*1000);
    Capacity(idx)   = mean(capacity_trials);
    OutageProb(idx) = mean(outage_flags);
end
figure;
plot(d, Pr_dBm, 'LineWidth', 2);
xlabel('Distance (m)'); ylabel('Received Power (dBm)');
title('Received Power vs Distance'); grid on;
figure;
plot(d, Capacity, 'LineWidth', 2);
xlabel('Distance (m)'); ylabel('Capacity (bps/Hz)');
title('MIMO Capacity vs Distance'); grid on;
figure;
semilogy(d, OutageProb, 'LineWidth', 2);
xlabel('Distance (m)'); ylabel('Outage Probability');
title(['Outage Probability (Target: ' num2str(R_target) ' bps/Hz)']);
grid on; ylim([1e-3 1]);
theta_deg = -90:1:90;
G = antennaGain(deg2rad(theta_deg));
figure;
polarplot(deg2rad(theta_deg), G, 'LineWidth', 2);
title('Antenna Radiation Pattern (dBi)');
