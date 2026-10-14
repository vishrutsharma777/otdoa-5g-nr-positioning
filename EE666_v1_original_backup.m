% EE 666 Final Project: Performance Analysis of OTDOA Positioning in 5G NR
%  Authors: Vishrut Sharma (25210157), Vikash Kumar Sinha (25210156)
%  IIT Gandhinagar
%
%  CHANNELS SIMULATED:
%   - AWGN only (baseline)
%   - 3GPP TDL-A (NLOS, indoor/urban, 100 ns delay spread)  <-- nrTDLChannel
%   - 3GPP TDL-C (NLOS, outdoor urban, 300 ns delay spread) <-- nrTDLChannel
%
%  BANDWIDTH CONFIGS:
%   - Config A : 52 RBs  ~18.7 MHz  (30 kHz SCS)
%   - Config B : 264 RBs ~95.0 MHz  (30 kHz SCS)
%
%  MONTE CARLO: 1000 trials, random UE positions
%
%  SECTIONS:
%   1. Scenario setup & constants
%   2. Corrected TDoA GDOP heatmap
%   3. CRLB vs bandwidth and SNR
%   4. 5G NR PRS generation (both bandwidth configs)
%   5. Single-point demo: AWGN vs TDL-A cross-correlation comparison
%   6. Single-point positioning demo across all 6 configurations
%   7. Monte Carlo: AWGN vs TDL-A vs TDL-C, 20 MHz vs 100 MHz
%   8. CRLB overlay validation plot


clear; clc; close all;
rng(42);

fprintf('  EE 666 -- 5G NR OTDOA Positioning Final Simulation\n');


%% SECTION 1: SCENARIO SETUP & PHYSICAL CONSTANTS
fprintf('[1/8] Setting up simulation scenario...\n');

c        = 3e8;
areaSize = 1000;

gNB_coords = [0,        0;
              0,        areaSize;
              areaSize, 0;
              areaSize, areaSize];

N_gNB  = size(gNB_coords, 1);
N_TDoA = N_gNB - 1;

ref_UE = [600, 380];

N_MC   = 1000;
margin = 50;
rng(42);
rand_UEs = margin + (areaSize - 2*margin) * rand(N_MC, 2);

SNR_dB_MC  = [0, 10, 20];
SNR_dB_vec = -10:2:30;
Delta_f    = 30e3;


%% SECTION 2: CORRECTED TDoA GDOP HEATMAP
fprintf('[2/8] Computing corrected TDoA GDOP heatmap...\n');

invC = eye(N_TDoA) - ones(N_TDoA) / (N_TDoA + 1);

gridRes  = 10;
[Xg, Yg] = meshgrid(0:gridRes:areaSize, 0:gridRes:areaSize);
GDOP_map  = zeros(size(Xg));

for ii = 1:size(Xg,1)
    for jj = 1:size(Xg,2)
        xu = Xg(ii,jj);
        yu = Yg(ii,jj);
        R  = max(sqrt((gNB_coords(:,1)-xu).^2 + (gNB_coords(:,2)-yu).^2), 1);
        H  = zeros(N_TDoA, 2);
        for k = 2:N_gNB
            H(k-1,1) = (xu-gNB_coords(k,1))/R(k) - (xu-gNB_coords(1,1))/R(1);
            H(k-1,2) = (yu-gNB_coords(k,2))/R(k) - (yu-gNB_coords(1,2))/R(1);
        end
        FIM = H' * invC * H + eye(2)*1e-9;
        GDOP_map(ii,jj) = sqrt(trace(inv(FIM)));
    end
end

R_ref = max(sqrt((gNB_coords(:,1)-ref_UE(1)).^2 + ...
                 (gNB_coords(:,2)-ref_UE(2)).^2), 1);
H_ref = zeros(N_TDoA, 2);
for k = 2:N_gNB
    H_ref(k-1,1) = (ref_UE(1)-gNB_coords(k,1))/R_ref(k) - ...
                   (ref_UE(1)-gNB_coords(1,1))/R_ref(1);
    H_ref(k-1,2) = (ref_UE(2)-gNB_coords(k,2))/R_ref(k) - ...
                   (ref_UE(2)-gNB_coords(1,2))/R_ref(1);
end
GDOP_ref = sqrt(trace(inv(H_ref' * invC * H_ref)));
fprintf('   GDOP at reference UE location = %.4f\n', GDOP_ref);

figure('Name','Fig1_GDOP','Position',[50 50 620 520]);
contourf(Xg, Yg, min(GDOP_map,4), 60, 'LineStyle','none');
colorbar; colormap(turbo); clim([0 4]);
hold on;
scatter(gNB_coords(:,1), gNB_coords(:,2), 180, '^', ...
        'MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1.5, ...
        'DisplayName','gNBs');
scatter(ref_UE(1), ref_UE(2), 140, 'o', ...
        'MarkerFaceColor','g','MarkerEdgeColor','k','LineWidth',1.5, ...
        'DisplayName','Ref UE');
title(sprintf('TDoA GDOP Heatmap (Corrected), central GDOP = %.2f', GDOP_ref));
xlabel('X (m)'); ylabel('Y (m)');
axis equal; axis([0 areaSize 0 areaSize]);
legend('Location','best'); hold off;


%% SECTION 3: CRLB vs BANDWIDTH AND SNR
fprintf('[3/8] Computing CRLB curves...\n');

BW_MHz_list = [5, 10, 20, 40, 100];
N_RB_list   = round(BW_MHz_list * 1e6 / (12 * Delta_f));
SNR_lin_vec = 10.^(SNR_dB_vec/10);
colors_bw   = lines(length(BW_MHz_list));

figure('Name','Fig2_CRLB','Position',[50 50 700 500]);
hold on;
for bIdx = 1:length(BW_MHz_list)
    N_SC  = N_RB_list(bIdx) * 12;
    k_idx = (-(N_SC/2):(N_SC/2)-1)';
    brms2 = (Delta_f^2) * mean(k_idx.^2);
    crlb_pos = GDOP_ref * c * sqrt(1 ./ (8*pi^2 * SNR_lin_vec * brms2));
    plot(SNR_dB_vec, crlb_pos, '-', 'Color', colors_bw(bIdx,:), ...
         'LineWidth', 2, 'DisplayName', sprintf('%d MHz', BW_MHz_list(bIdx)));
end
set(gca,'YScale','log'); grid on; grid minor;
xlabel('SNR (dB)'); ylabel('CRLB on 2D Position Error (m)');
title('Cramer-Rao Lower Bound vs. SNR for Different PRS Bandwidths');
legend('Location','southwest'); hold off;


%% SECTION 4: 5G NR PRS GENERATION (STANDARD-COMPLIANT)
fprintf('[4/8] Generating 5G NR PRS waveforms (both BW configs)...\n');

carrierA               = nrCarrierConfig;
carrierA.SubcarrierSpacing = 30;
carrierA.NSizeGrid     = 52;
Fs_A                   = nrOFDMInfo(carrierA).SampleRate;

carrierB               = nrCarrierConfig;
carrierB.SubcarrierSpacing = 30;
carrierB.NSizeGrid     = 264;
Fs_B                   = nrOFDMInfo(carrierB).SampleRate;

fprintf('   Config A: %3d RBs, Fs = %6.2f MHz (~%3.0f MHz PRS BW)\n', ...
        carrierA.NSizeGrid, Fs_A/1e6, carrierA.NSizeGrid*12*Delta_f/1e6);
fprintf('   Config B: %3d RBs, Fs = %6.2f MHz (~%3.0f MHz PRS BW)\n', ...
        carrierB.NSizeGrid, Fs_B/1e6, carrierB.NSizeGrid*12*Delta_f/1e6);

[txA, delaysA] = generate_prs(carrierA, gNB_coords, ref_UE, c);
[txB, delaysB] = generate_prs(carrierB, gNB_coords, ref_UE, c);

prs_demo              = make_prs_config(carrierA, 1);
prsGrid               = zeros(carrierA.NSizeGrid*12, 14);
prsGrid(nrPRSIndices(carrierA, prs_demo)) = nrPRS(carrierA, prs_demo);

figure('Name','Fig3_PRS_Grid','Position',[50 50 600 420]);
imagesc(abs(prsGrid)); set(gca,'YDir','normal');
colormap(parula); colorbar;
title('5G NR PRS Resource Grid - gNB 1 (CombSize=4, 52 RBs, 30 kHz SCS)');
xlabel('OFDM Symbol Index'); ylabel('Subcarrier Index');
occ = nnz(prsGrid) / numel(prsGrid);
text(0.5, 0.95, sprintf('PRS occupancy: %.1f%%', occ*100), ...
     'Units','normalized','Color','w','FontSize',10,'FontWeight','bold', ...
     'HorizontalAlignment','center');


%% SECTION 5: CROSS-CORRELATION — AWGN vs TDL-A, 20 MHz vs 100 MHz
fprintf('[5/8] Generating CCF comparison plots (AWGN vs TDL-A)...\n');

SNR_demo = 20;

[rxA_awgn, toaA_awgn] = rx_awgn(txA, delaysA, SNR_demo);
[rxA_tdla, toaA_tdla] = rx_tdl(txA, delaysA, SNR_demo, Fs_A, 'TDL-A', 100e-9);
[rxB_awgn, toaB_awgn] = rx_awgn(txB, delaysB, SNR_demo);
[rxB_tdla, toaB_tdla] = rx_tdl(txB, delaysB, SNR_demo, Fs_B, 'TDL-A', 100e-9);

figure('Name','Fig4_CCF_Comparison','Position',[50 50 920 560]);

combo_rx  = {rxA_awgn, rxA_tdla, rxB_awgn, rxB_tdla};
combo_tx  = {txA,      txA,      txB,      txB     };
combo_dl  = {delaysA,  delaysA,  delaysB,  delaysB };
combo_lbl = {'20 MHz + AWGN', '20 MHz + TDL-A (100 ns DS)', ...
             '100 MHz + AWGN', '100 MHz + TDL-A (100 ns DS)'};
win_half  = [80, 80, 350, 350];

for s = 1:4
    [co, lg] = xcorr(combo_rx{s}, combo_tx{s}{1});
    co       = abs(co) / max(abs(co));
    pk       = combo_dl{s}(1);
    ww       = max(1, find(lg >= pk-win_half(s), 1)) : ...
               min(length(lg), find(lg <= pk+win_half(s), 1, 'last'));

    subplot(2,2,s);
    plot(lg(ww) - pk, co(ww), 'b', 'LineWidth', 1.5);
    xline(0, '--r', 'True LoS', 'LineWidth', 1.2, 'LabelVerticalAlignment','bottom');
    grid on; ylim([0 1.15]);
    title(combo_lbl{s}, 'FontSize', 9);
    xlabel('Lag offset (samples)'); ylabel('Norm. Correlation');
end
sgtitle('Cross-Correlation Peak: BW and Multipath Channel Effect (gNB 1)', ...
        'FontWeight','bold');


%% SECTION 6: SINGLE-POINT POSITIONING DEMO (ALL 6 CONFIGURATIONS)
fprintf('[6/8] Single-point positioning demo (SNR = %d dB)...\n', SNR_demo);

[rxA_tdlc, toaA_tdlc] = rx_tdl(txA, delaysA, SNR_demo, Fs_A, 'TDL-C', 300e-9);
[rxB_tdlc, toaB_tdlc] = rx_tdl(txB, delaysB, SNR_demo, Fs_B, 'TDL-C', 300e-9);

configs_demo = {
    toaA_awgn, gNB_coords, ref_UE, c, Fs_A, '~20 MHz, AWGN';
    toaA_tdla, gNB_coords, ref_UE, c, Fs_A, '~20 MHz, TDL-A (100 ns)';
    toaA_tdlc, gNB_coords, ref_UE, c, Fs_A, '~20 MHz, TDL-C (300 ns)';
    toaB_awgn, gNB_coords, ref_UE, c, Fs_B, '~100 MHz, AWGN';
    toaB_tdla, gNB_coords, ref_UE, c, Fs_B, '~100 MHz, TDL-A (100 ns)';
    toaB_tdlc, gNB_coords, ref_UE, c, Fs_B, '~100 MHz, TDL-C (300 ns)';
};

est_pos_all = zeros(6, 2);
err_all     = zeros(6, 1);
for i = 1:6
    [est_pos_all(i,:), err_all(i)] = tdoa_nlls( ...
        configs_demo{i,1}, configs_demo{i,2}, configs_demo{i,3}, ...
        configs_demo{i,4}, configs_demo{i,5}, areaSize);
end

fprintf('\n   %-32s  Error (m)\n', 'Configuration');
fprintf('   %-32s  ---------\n', repmat('-',1,32));
for i = 1:6
    fprintf('   %-32s  %.2f\n', configs_demo{i,6}, err_all(i));
end

figure('Name','Fig5_Positioning_Map','Position',[50 50 640 540]);
hold on;
scatter(gNB_coords(:,1), gNB_coords(:,2), 200, '^r', 'filled', 'DisplayName','gNBs');
scatter(ref_UE(1), ref_UE(2), 180, 'og', 'filled', 'DisplayName','True UE');
mkr_styles = {'sb','db','sb','sd','dd','sd'};
mkr_colors = {[0.1 0.4 0.9],[0.1 0.4 0.9],[0.1 0.4 0.9], ...
              [0.9 0.2 0.1],[0.9 0.2 0.1],[0.9 0.2 0.1]};
for i = 1:6
    scatter(est_pos_all(i,1), est_pos_all(i,2), 120, ...
            mkr_styles{i}(1), 'LineWidth', 2, 'MarkerEdgeColor', mkr_colors{i}, ...
            'DisplayName', sprintf('%s (%.1fm)', configs_demo{i,6}, err_all(i)));
end
grid on; axis equal; axis([-100 areaSize+100 -100 areaSize+100]);
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Positioning Map - All Configurations (SNR = %d dB)', SNR_demo));
legend('Location','best','FontSize',7); hold off;


%% SECTION 7: MONTE CARLO — 1000 TRIALS
fprintf('\n[7/8] Running Monte Carlo (%d trials x 6 configs x 3 SNRs)...\n', N_MC);
fprintf('      [~18000 total runs -- expect several minutes]\n\n');

ch_names   = {'AWGN', 'TDL-A', 'TDL-C'};
ch_types   = {'awgn', 'tdl',   'tdl'  };
ch_prof    = {'',     'TDL-A', 'TDL-C'};
ch_ds      = {0,      100e-9,  300e-9 };

bw_labels  = {'20 MHz',  '100 MHz'};
bw_carrier = {carrierA,   carrierB};
bw_Fs      = {Fs_A,       Fs_B    };

errors = cell(2, 3, length(SNR_dB_MC));

total_combos = 2 * 3 * length(SNR_dB_MC);
combo_num    = 0;

for bIdx = 1:2
    carr_mc = bw_carrier{bIdx};
    Fs_mc   = bw_Fs{bIdx};

    for cIdx = 1:3
        for sIdx = 1:length(SNR_dB_MC)
            combo_num = combo_num + 1;
            snr_db    = SNR_dB_MC(sIdx);
            err_vec   = zeros(N_MC, 1);

            fprintf('   [%2d/%2d] %-8s + %-20s SNR=%2d dB ... ', ...
                    combo_num, total_combos, ...
                    bw_labels{bIdx}, ch_names{cIdx}, snr_db);

            for trial = 1:N_MC
                ue_pos = rand_UEs(trial, :);

                [tx_mc, dl_mc] = generate_prs(carr_mc, gNB_coords, ue_pos, c);

                if strcmp(ch_types{cIdx}, 'awgn')
                    [~, toa_mc] = rx_awgn(tx_mc, dl_mc, snr_db);
                else
                    [~, toa_mc] = rx_tdl(tx_mc, dl_mc, snr_db, Fs_mc, ...
                                         ch_prof{cIdx}, ch_ds{cIdx});
                end

                [~, err_vec(trial)] = tdoa_nlls(toa_mc, gNB_coords, ue_pos, ...
                                                c, Fs_mc, areaSize);
            end

            errors{bIdx, cIdx, sIdx} = err_vec;
            fprintf('median=%5.1f m  90th=%6.1f m\n', ...
                    median(err_vec), prctile(err_vec, 90));
        end
    end
end


%% SECTION 8: CDF PLOTS AND CRLB VALIDATION
fprintf('\n[8/8] Generating final CDF and validation plots...\n');

N_pts  = N_MC;
x_cdf  = (1:N_pts) / N_pts;

ch_ls   = {'-', '--', ':'};
col_20  = [0.12 0.47 0.71];
col_100 = [0.84 0.15 0.16];
bw_col  = {col_20, col_100};

snr20_idx = find(SNR_dB_MC == 20);
figure('Name','Fig6a_CDF_BW','Position',[50 50 850 520]);
hold on;
for bIdx = 1:2
    for cIdx = 1:3
        es = sort(errors{bIdx, cIdx, snr20_idx});
        plot(es, x_cdf, ...
             'Color',     bw_col{bIdx}, ...
             'LineStyle', ch_ls{cIdx}, ...
             'LineWidth', 2, ...
             'DisplayName', sprintf('%s + %s', bw_labels{bIdx}, ch_names{cIdx}));
    end
end
xline(1, ':k', 'LineWidth', 1, 'HandleVisibility','off');
grid on; xlim([0 200]);
xlabel('2D Positioning Error (m)'); ylabel('CDF');
title(sprintf('CDF of Positioning Error: Bandwidth and Channel Effect  (SNR=20 dB, N=%d)', N_MC));
legend('Location','southeast','FontSize',8,'NumColumns',2); hold off;

snr_colors = {[0.80 0.10 0.10], [0.10 0.55 0.10], [0.10 0.20 0.80]};
figure('Name','Fig6b_CDF_SNR','Position',[50 50 850 520]);
hold on;
for sIdx = 1:length(SNR_dB_MC)
    for cIdx = 1:3
        es = sort(errors{2, cIdx, sIdx});
        plot(es, x_cdf, ...
             'Color',     snr_colors{sIdx}, ...
             'LineStyle', ch_ls{cIdx}, ...
             'LineWidth', 2, ...
             'DisplayName', sprintf('SNR=%d dB + %s', SNR_dB_MC(sIdx), ch_names{cIdx}));
    end
end
grid on; xlim([0 400]);
xlabel('2D Positioning Error (m)'); ylabel('CDF');
title(sprintf('CDF of Positioning Error: SNR and Channel Effect  (~100 MHz PRS, N=%d)', N_MC));
legend('Location','southeast','FontSize',8); hold off;

figure('Name','Fig7_CRLB_Validation','Position',[50 50 720 500]);
hold on;
SNR_lin_sm = 10.^(SNR_dB_vec/10);

for bIdx = 1:2
    N_SC  = bw_carrier{bIdx}.NSizeGrid * 12;
    k_idx = (-(N_SC/2):(N_SC/2)-1)';
    brms2 = (Delta_f^2) * mean(k_idx.^2);

    crlb_sm = GDOP_ref * c * sqrt(1 ./ (8*pi^2 * SNR_lin_sm * brms2));
    plot(SNR_dB_vec, crlb_sm, '-', 'Color', bw_col{bIdx}, 'LineWidth', 2.5, ...
         'DisplayName', sprintf('CRLB (AWGN) %s', bw_labels{bIdx}));

    med_awgn = arrayfun(@(s) median(errors{bIdx,1,s}), 1:length(SNR_dB_MC));
    med_tdla = arrayfun(@(s) median(errors{bIdx,2,s}), 1:length(SNR_dB_MC));
    med_tdlc = arrayfun(@(s) median(errors{bIdx,3,s}), 1:length(SNR_dB_MC));

    plot(SNR_dB_MC, med_awgn, 'o--', 'Color', bw_col{bIdx}, 'LineWidth',2, ...
         'MarkerSize',9, 'DisplayName', sprintf('MC (AWGN) %s', bw_labels{bIdx}));
    plot(SNR_dB_MC, med_tdla, 's:', 'Color', bw_col{bIdx}*0.7, 'LineWidth',2, ...
         'MarkerSize',9, 'DisplayName', sprintf('MC (TDL-A) %s', bw_labels{bIdx}));
    plot(SNR_dB_MC, med_tdlc, '^-.', 'Color', bw_col{bIdx}*0.5, 'LineWidth',2, ...
         'MarkerSize',9, 'DisplayName', sprintf('MC (TDL-C) %s', bw_labels{bIdx}));
end
set(gca,'YScale','log'); grid on; grid minor;
xlabel('SNR (dB)'); ylabel('Median Positioning Error (m)');
title('Monte Carlo Median vs. AWGN CRLB: Multipath Degradation');
legend('Location','southwest','FontSize',7,'NumColumns',2); hold off;

fprintf('  Simulation Complete. %d Monte Carlo trials per config.\n', N_MC);
fprintf('  Save all figures as PNG/EPS for the LaTeX report.\n');


%% LOCAL HELPER FUNCTIONS
function [txWaveforms, true_delays] = generate_prs(carrier, gNB_coords, ue_pos, c)
    ofdmInfo        = nrOFDMInfo(carrier);
    Fs              = ofdmInfo.SampleRate;
    N_gNB           = size(gNB_coords, 1);
    txWaveforms     = cell(N_gNB, 1);
    true_delays     = zeros(N_gNB, 1);
    for i = 1:N_gNB
        prs              = make_prs_config(carrier, i);
        rg               = zeros(carrier.NSizeGrid*12, 14);
        rg(nrPRSIndices(carrier, prs)) = nrPRS(carrier, prs);
        txWaveforms{i}   = nrOFDMModulate(carrier, rg);
        dist             = norm(gNB_coords(i,:) - ue_pos);
        true_delays(i)   = round((dist/c) * Fs);
    end
end

function prs = make_prs_config(carrier, cell_id)
    prs                      = nrPRSConfig;
    prs.NumRB                = carrier.NSizeGrid;
    prs.PRSResourceSetPeriod = [10 0];
    prs.NumPRSSymbols        = 4;
    prs.SymbolStart          = 8;
    prs.NPRSID               = cell_id;
    prs.CombSize             = 4;
end

function [rxSignal, est_ToA] = rx_awgn(txWaveforms, true_delays, SNR_dB)
    N_gNB    = length(txWaveforms);
    sig_len  = length(txWaveforms{1}) + max(true_delays) + 200;
    rxSignal = zeros(sig_len, 1);
    for i = 1:N_gNB
        padded = [zeros(true_delays(i),1); txWaveforms{i}];
        rxSignal(1:length(padded)) = rxSignal(1:length(padded)) + padded;
    end
    rxSignal = awgn(rxSignal, SNR_dB, 'measured');
    est_ToA  = estimate_toa(rxSignal, txWaveforms);
end

function [rxSignal, est_ToA] = rx_tdl(txWaveforms, true_delays, SNR_dB, ...
                                       Fs, profile, delay_spread)
    N_gNB    = length(txWaveforms);
    tail     = ceil(10 * delay_spread * Fs) + 500;
    sig_len  = length(txWaveforms{1}) + max(true_delays) + tail;
    rxSignal = zeros(sig_len, 1);

    for i = 1:N_gNB
        ch                         = nrTDLChannel;
        ch.SampleRate              = Fs;
        ch.DelayProfile            = profile;
        ch.DelaySpread             = delay_spread;
        ch.MaximumDopplerShift     = 0;
        ch.NumTransmitAntennas     = 1;
        ch.NumReceiveAntennas      = 1;
        ch.NormalizePathGains      = true;
        ch.NormalizeChannelOutputs = false; 

        [faded, ~] = ch(txWaveforms{i});

        padded = [zeros(true_delays(i), 1); faded];
        rxSignal(1:length(padded)) = rxSignal(1:length(padded)) + padded;
    end

    rxSignal = awgn(rxSignal, SNR_dB, 'measured');
    est_ToA  = estimate_toa(rxSignal, txWaveforms);
end

function est_ToA = estimate_toa(rxSignal, txWaveforms)
    N_gNB   = length(txWaveforms);
    est_ToA = zeros(N_gNB, 1);
    for i = 1:N_gNB
        [corr_out, lags] = xcorr(rxSignal, txWaveforms{i});
        [~, mx]          = max(abs(corr_out));
        est_ToA(i)       = lags(mx);
    end
end

function [est_pos, err] = tdoa_nlls(est_ToA, gNB_coords, true_pos, ...
                                    c, Fs, areaSize)
    rstd_m = (est_ToA(2:end) - est_ToA(1)) * (c/Fs);

    cost_fn = @(p) sum(( ...
        (sqrt((gNB_coords(2:end,1)-p(1)).^2+(gNB_coords(2:end,2)-p(2)).^2) - ...
         sqrt((gNB_coords(1,1)-p(1))^2+(gNB_coords(1,2)-p(2))^2)) ...
        - rstd_m ).^2);

    opts    = optimset('Display','off','TolX',1e-6,'TolFun',1e-10,'MaxIter',5000);
    est_pos = fminsearch(cost_fn, [areaSize/2, areaSize/2], opts);
    err     = norm(true_pos - est_pos);
end


%% SECTION 9: DOPPLER SHIFT ANALYSIS (UE MOBILITY)
fprintf('\n[9/10] Running Doppler Monte Carlo (TDL-A, SNR=20 dB)...\n');

doppler_vals  = [0, 10, 100, 400];     % Hz
doppler_labels = {'Static (0 Hz)', 'Pedestrian (10 Hz)', ...
                  'Vehicular (100 Hz)', 'High-speed (400 Hz)'};
SNR_doppler   = 20;

errors_dop = cell(2, length(doppler_vals));

for bIdx = 1:2
    carr_d = bw_carrier{bIdx};
    Fs_d   = bw_Fs{bIdx};

    for dIdx = 1:length(doppler_vals)
        fd     = doppler_vals(dIdx);
        ev     = zeros(N_MC, 1);

        fprintf('   %s + %s, f_D=%3d Hz ... ', ...
                bw_labels{bIdx}, 'TDL-A', fd);

        for trial = 1:N_MC
            ue_pos = rand_UEs(trial,:);
            [tx_d, dl_d] = generate_prs(carr_d, gNB_coords, ue_pos, c);
            [~, toa_d]   = rx_tdl_doppler(tx_d, dl_d, SNR_doppler, ...
                                           Fs_d, 'TDL-A', 100e-9, fd);
            [~, ev(trial)] = tdoa_nlls(toa_d, gNB_coords, ue_pos, ...
                                       c, Fs_d, areaSize);
        end

        errors_dop{bIdx, dIdx} = ev;
        fprintf('median=%.2f m  90th=%.2f m\n', ...
                median(ev), prctile(ev,90));
    end
end

%% SECTION 9b: CCF SHAPE UNDER DOPPLER (gNB 1, single realisation)
fprintf('\n[10/10] Generating Doppler CCF and CDF plots...\n');

figure('Name','Fig8_CCF_Doppler','Position',[50 50 900 560]);
dop_colors = {[0.2 0.7 0.2], [0.2 0.2 0.9], [0.9 0.6 0.1], [0.9 0.1 0.1]};

for dIdx = 1:4
    fd = doppler_vals(dIdx);

    [rx_ccf, ~] = rx_tdl_doppler(txA, delaysA, SNR_demo, Fs_A, ...
                                  'TDL-A', 100e-9, fd);
    [co, lg] = xcorr(rx_ccf, txA{1});
    co = abs(co)/max(abs(co));
    pk = delaysA(1);
    ww = find(lg >= pk-80, 1) : find(lg <= pk+80, 1,'last');

    subplot(2,4,dIdx);
    plot(lg(ww)-pk, co(ww), 'Color', dop_colors{dIdx}, 'LineWidth',1.5);
    xline(0,'--r','LineWidth',1);
    grid on; ylim([0 1.15]);
    title(sprintf('20 MHz\n%s', doppler_labels{dIdx}), 'FontSize',8);
    xlabel('Lag offset (samples)'); ylabel('Norm. Corr.');

    [rx_ccf2, ~] = rx_tdl_doppler(txB, delaysB, SNR_demo, Fs_B, ...
                                   'TDL-A', 100e-9, fd);
    [co2, lg2] = xcorr(rx_ccf2, txB{1});
    co2 = abs(co2)/max(abs(co2));
    pk2 = delaysB(1);
    ww2 = find(lg2 >= pk2-350, 1) : find(lg2 <= pk2+350, 1,'last');

    subplot(2,4,dIdx+4);
    plot(lg2(ww2)-pk2, co2(ww2), 'Color', dop_colors{dIdx}, 'LineWidth',1.5);
    xline(0,'--r','LineWidth',1);
    grid on; ylim([0 1.15]);
    title(sprintf('100 MHz\n%s', doppler_labels{dIdx}), 'FontSize',8);
    xlabel('Lag offset (samples)'); ylabel('Norm. Corr.');
end
sgtitle('CCF Peak Shape: TDL-A, Doppler Effect (gNB 1, SNR=20 dB)', ...
        'FontWeight','bold');

%% SECTION 9c: CDF — DOPPLER EFFECT
figure('Name','Fig9_CDF_Doppler','Position',[50 50 860 520]);
dop_ls     = {'-', '--', ':', '-.'};
bw_col_dop = {[0.12 0.47 0.71], [0.84 0.15 0.16]};

hold on;
for bIdx = 1:2
    for dIdx = 1:length(doppler_vals)
        es = sort(errors_dop{bIdx, dIdx});
        plot(es, x_cdf, ...
             'Color',     bw_col_dop{bIdx}, ...
             'LineStyle', dop_ls{dIdx}, ...
             'LineWidth', 2, ...
             'DisplayName', sprintf('%s + %s', bw_labels{bIdx}, doppler_labels{dIdx}));
    end
end
grid on; xlim([0 30]);
xlabel('2D Positioning Error (m)'); ylabel('CDF');
title(sprintf('CDF of Positioning Error: Doppler Effect  (TDL-A, SNR=20 dB, N=%d)', N_MC));
legend('Location','southeast','FontSize',7,'NumColumns',2); hold off;

%% SECTION 9d: MEDIAN AND 90th PERCENTILE vs DOPPLER
figure('Name','Fig10_Doppler_vs_BW','Position',[50 50 700 480]);
hold on;

pct50_A = cellfun(@(e) median(e),    errors_dop(1,:));
pct90_A = cellfun(@(e) prctile(e,90),errors_dop(1,:));
pct50_B = cellfun(@(e) median(e),    errors_dop(2,:));
pct90_B = cellfun(@(e) prctile(e,90),errors_dop(2,:));

fd_plot = doppler_vals + 0.5;

plot(fd_plot, pct50_A, 'o-',  'Color', bw_col_dop{1}, 'LineWidth',2, ...
     'MarkerSize',9, 'DisplayName','20 MHz — Median');
plot(fd_plot, pct90_A, 's--', 'Color', bw_col_dop{1}, 'LineWidth',2, ...
     'MarkerSize',9, 'DisplayName','20 MHz — 90th pct');
plot(fd_plot, pct50_B, 'o-',  'Color', bw_col_dop{2}, 'LineWidth',2, ...
     'MarkerSize',9, 'DisplayName','100 MHz — Median');
plot(fd_plot, pct90_B, 's--', 'Color', bw_col_dop{2}, 'LineWidth',2, ...
     'MarkerSize',9, 'DisplayName','100 MHz — 90th pct');

set(gca,'XScale','log','YScale','log');
xticks([0.5 10 100 400]);
xticklabels({'0 (static)','10 (ped.)','100 (veh.)','400 (high-speed)'});
grid on; grid minor;
xlabel('Doppler Shift f_D (Hz)'); ylabel('Positioning Error (m)');
title('50th and 90th Percentile Error vs. Doppler Shift (TDL-A, SNR=20 dB)');
legend('Location','northwest','FontSize',9); hold off;

fprintf('\n=== All Doppler figures generated. ===\n');


%% ADDITIONAL LOCAL FUNCTION: rx_tdl_doppler


function [rxSignal, est_ToA] = rx_tdl_doppler(txWaveforms, true_delays, ...
                                               SNR_dB, Fs, profile, ...
                                               delay_spread, doppler_hz)
    N_gNB    = length(txWaveforms);
    tail     = ceil(10 * delay_spread * Fs) + 500;
    sig_len  = length(txWaveforms{1}) + max(true_delays) + tail;
    rxSignal = zeros(sig_len, 1);

    for i = 1:N_gNB
        ch                         = nrTDLChannel;
        ch.SampleRate              = Fs;
        ch.DelayProfile            = profile;
        ch.DelaySpread             = delay_spread;
        ch.MaximumDopplerShift     = doppler_hz;
        ch.NumTransmitAntennas     = 1;
        ch.NumReceiveAntennas      = 1;
        ch.NormalizePathGains      = true;
        ch.NormalizeChannelOutputs = false;

        [faded, ~] = ch(txWaveforms{i});

        padded = [zeros(true_delays(i), 1); faded];
        rxSignal(1:length(padded)) = rxSignal(1:length(padded)) + padded;
    end

    rxSignal = awgn(rxSignal, SNR_dB, 'measured');
    est_ToA  = estimate_toa(rxSignal, txWaveforms);
end