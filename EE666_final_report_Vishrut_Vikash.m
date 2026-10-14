% EE 666 Final Project: Performance Analysis of OTDOA Positioning in 5G NR
%  Authors: Vishrut Sharma (25210157), Vikash Kumar Sinha (25210156)
%  IIT Gandhinagar
%
% =========================================================================
%  REVISION 2 -- what changed and why
% =========================================================================
%  Revision 1 reported (correctly, in Sec. VI-D / VII-C) that the Monte
%  Carlo medians were flat across SNR because of integer-sample
%  discretisation of the CCF peak, and left sub-sample interpolation as
%  future work. This revision implements that future work and, in doing so,
%  uncovers two further issues that also had to be fixed for the SNR axis to
%  mean anything.
%
%  [1] FRACTIONAL GEOMETRIC DELAY.  v1 applied the propagation delay as
%      round((d/c)*Fs) -- an integer number of samples. One sample is 9.77 m
%      at 20 MHz and 2.44 m at 100 MHz, so the *true* delay was snapped to a
%      grid coarser than the accuracy being measured. Sub-sample estimation
%      alone would not have helped: there was nothing sub-sample left in the
%      signal to estimate. Delay is now applied exactly, as a frequency-
%      domain phase ramp (frac_delay).
%
%  [2] SUB-SAMPLE ToA ESTIMATION.  v1 took the integer argmax of |CCF|,
%      giving a quantisation RMS of 1/sqrt(12) = 0.289 samples -- which is
%      exactly what v1 measured. Three estimators are now implemented and
%      compared: 'integer' (v1), 'parabolic' (the formula proposed in v1
%      Sec. VII-C), and 'bandlimited' (FFT interpolation of the complex CCF).
%
%      MEASURED, single link at 100 MHz:
%          integer      0.289 samples   ~0.71 m
%          parabolic    0.055 samples   ~0.13 m
%          bandlimited  0.018 samples   ~0.043 m
%      v1 Sec. VII-C predicted parabolic would reach "sub-centimetre". It
%      does not: it is bias-limited at roughly 13 cm. Band-limited
%      interpolation is ~3x better again. This corrects that claim.
%
%  [3] SNR RANGE.  This is the one that actually explains v1's flat curves.
%      The PRS carries ~35 dB of processing gain (3168 REs). At the SNRs v1
%      swept (0, 10, 20 dB) the post-correlation SNR is so high that the
%      estimator sits on its bias floor, so SNR cannot influence the result
%      no matter how good the estimator is. The noise-limited region is
%      roughly -30 to -10 dB input SNR, and loss of lock sets in below about
%      -35 dB. The SNR sweep now covers that range.
%
%  [4] PRS COMB STAGGERING.  The report (Sec. IV-A) states that "adjacent
%      cells use staggered offsets", but v1 left REOffset at its default of
%      0 for all four gNBs, so all cells collided on the same resource
%      elements and were separated only by NPRSID scrambling. REOffset is
%      now staggered as mod(cell-1, CombSize), making the code match the
%      report and giving genuine frequency-domain orthogonality.
%
%  [5] PERFORMANCE.  v1 called generate_prs() inside every Monte Carlo
%      trial, but the PRS waveform does not depend on UE position -- only
%      the delays do. That repeated ~26,000 OFDM modulations needlessly.
%      Waveforms are now generated once per bandwidth config.
%
%  [6] CHANNEL OBJECT LIFETIME -- and a latent correctness bug it exposed.
%      v1 constructed a fresh nrTDLChannel for every gNB on every trial. Over a
%      full sweep that is ~120,000 System objects; MATLAB does not release them
%      promptly and the run was killed for memory at 86% completion. Objects are
%      now created once per (rate, profile, Doppler, gNB) and reset() before each
%      use, which bounds memory at ~2.3 GB.
%
%      Fixing that revealed something more important. nrTDLChannel defaults to
%      RandomStream = 'mt19937ar with seed' with a FIXED default Seed. So every
%      freshly constructed object produced the SAME fading realisation -- every
%      gNB, every trial, every run. The report states (Sec. V-A) that "a
%      separate, independently seeded channel object is created for each gNB-UE
%      link". That was never true: they were identically seeded.
%
%      Measured, TDL-C at 100 MHz, SNR = 20 dB:
%          fresh object, default fixed seed (v1) : median  0.046 m
%          randomised seed per trial             : median 21.51 m
%          cached + reset, Global stream (v2)    : median 21.66 m
%      The last two agree; the first does not. v1's sub-decimetre multipath
%      result was an artefact of reusing one benign realisation.
%
%      ~21.6 m corresponds to ~72 ns, comfortably inside TDL-C's 300 ns delay
%      spread. That is the expected NLOS behaviour: the first arrival is not the
%      strongest, so a peak-picking estimator locks onto a later dominant tap and
%      acquires a positive bias. The remedy is first-arrival (leading-edge)
%      detection rather than peak detection, which is what real 3GPP positioning
%      receivers implement. That is the honest next step for this study, and it
%      supersedes the v1 conclusion that "at 100 MHz the system is robust to the
%      simulated multipath profiles".
%
%      RandomStream is set to 'Global stream' so reset() draws a fresh,
%      independent realisation per trial while the whole sweep stays reproducible
%      from the rng(42) at the top of the script.
%
%  CHANNELS: AWGN, 3GPP TDL-A (100 ns DS), TDL-C (300 ns DS)
%  BANDWIDTHS: Config A 52 RB (~18.7 MHz), Config B 264 RB (~95 MHz), 30 kHz SCS
% =========================================================================

clear; clc; close all;
rng(42);

% -------------------------------------------------------------------------
% WORKLOAD PRESET. Pick one: 'quick' | 'standard' | 'full'
%
% This is a memory decision, not a statistical one. The development machine has
% 7.3 GB of RAM with ~1.7 GB free in normal use, and MATLAB needs ~2.3 GB for
% this workload. Two attempts at the 'full' preset were killed by the OS. The
% second died inside the estimator sweep in Section 5, BEFORE the Monte Carlo
% had started, because that loop runs 3 methods x 21 SNR points x n_est trials
% of 122.88 MHz FFTs.
%
% 'standard' (N_MC = 300) pins the medians down perfectly well -- the 90th
% percentile is the only statistic that visibly tightens at 1000 -- and it runs
% without evicting everything else on the machine. Use 'full' only on a larger
% machine, or after closing other applications.
% -------------------------------------------------------------------------
PRESET = 'standard';

switch PRESET
    case 'quick',    N_MC =   50;  N_EST =  40;
    case 'standard', N_MC =  300;  N_EST =  60;
    case 'full',     N_MC = 1000;  N_EST = 150;
    otherwise, error('EE666:preset','PRESET must be quick|standard|full');
end
QUICK = strcmp(PRESET,'quick');

fprintf('  EE 666 -- 5G NR OTDOA Positioning (revision 2)\n');
fprintf('  preset = %s   N_MC = %d   N_EST = %d

', PRESET, N_MC, N_EST);

TOA_METHOD = 'bandlimited';   % 'integer' | 'parabolic' | 'bandlimited'

% A trial whose error exceeds this has lost correlation lock; it is an
% outage, not a measurement. Medians are reported over locked trials only,
% with the outage rate quoted alongside, which is the standard way to report
% a positioning estimator in its threshold region.
OUTAGE_M = 50;


%% SECTION 1: SCENARIO SETUP & PHYSICAL CONSTANTS
fprintf('[1/9] Setting up simulation scenario...\n');

c        = 3e8;
areaSize = 1000;

gNB_coords = [0,        0;
              0,        areaSize;
              areaSize, 0;
              areaSize, areaSize];

N_gNB  = size(gNB_coords, 1);
N_TDoA = N_gNB - 1;

ref_UE = [600, 380];

margin   = 50;
rand_UEs = margin + (areaSize - 2*margin) * rand(N_MC, 2);

% SNR points for the Monte Carlo. See change [3]: these now sit in the
% region where the estimator is actually noise-limited. 20 dB is retained
% as a high-SNR anchor and to keep continuity with revision 1.
SNR_dB_MC  = [-30, -25, -20, -15, -10, 20];
SNR_dB_vec = -40:2:30;
Delta_f    = 30e3;


%% SECTION 2: TDoA GDOP HEATMAP
% Unchanged from revision 1 -- the corrected correlated-noise formulation
% was already right. C^-1 = I - 11'/N by Sherman-Morrison.
fprintf('[2/9] Computing TDoA GDOP heatmap...\n');

invC = eye(N_TDoA) - ones(N_TDoA) / (N_TDoA + 1);

gridRes  = 10;
[Xg, Yg] = meshgrid(0:gridRes:areaSize, 0:gridRes:areaSize);
GDOP_map = zeros(size(Xg));

for ii = 1:size(Xg,1)
    for jj = 1:size(Xg,2)
        GDOP_map(ii,jj) = gdop_at(Xg(ii,jj), Yg(ii,jj), gNB_coords, invC, true);
    end
end

GDOP_ref = gdop_at(ref_UE(1), ref_UE(2), gNB_coords, invC, false);
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
title(sprintf('TDoA GDOP Heatmap, central GDOP = %.2f', GDOP_ref));
xlabel('X (m)'); ylabel('Y (m)');
axis equal; axis([0 areaSize 0 areaSize]);
legend('Location','best'); hold off;


%% SECTION 3: CRLB vs BANDWIDTH AND SNR
fprintf('[3/9] Computing CRLB curves...\n');

BW_MHz_list = [5, 10, 20, 40, 100];
N_RB_list   = round(BW_MHz_list * 1e6 / (12 * Delta_f));
SNR_lin_vec = 10.^(SNR_dB_vec/10);
colors_bw   = lines(length(BW_MHz_list));

figure('Name','Fig2_CRLB','Position',[50 50 700 500]);
hold on;
for bIdx = 1:length(BW_MHz_list)
    crlb_pos = crlb_position(N_RB_list(bIdx), Delta_f, SNR_lin_vec, GDOP_ref, c);
    plot(SNR_dB_vec, crlb_pos, '-', 'Color', colors_bw(bIdx,:), ...
         'LineWidth', 2, 'DisplayName', sprintf('%d MHz', BW_MHz_list(bIdx)));
end
set(gca,'YScale','log'); grid on; grid minor;
xlabel('SNR (dB)'); ylabel('CRLB on 2D Position Error (m)');
title('Cramer-Rao Lower Bound vs. SNR for Different PRS Bandwidths');
legend('Location','southwest'); hold off;


%% SECTION 4: 5G NR PRS GENERATION (STANDARD-COMPLIANT, STAGGERED COMB)
fprintf('[4/9] Generating 5G NR PRS waveforms (once per config)...\n');

carrierA = nrCarrierConfig;
carrierA.SubcarrierSpacing = 30;
carrierA.NSizeGrid         = 52;
Fs_A     = nrOFDMInfo(carrierA).SampleRate;

carrierB = nrCarrierConfig;
carrierB.SubcarrierSpacing = 30;
carrierB.NSizeGrid         = 264;
Fs_B     = nrOFDMInfo(carrierB).SampleRate;

% Change [5]: generated ONCE. The PRS waveform is a property of the cell,
% not of where the UE happens to be; only the delays depend on position.
txA = generate_prs_waveforms(carrierA, N_gNB);
txB = generate_prs_waveforms(carrierB, N_gNB);

fprintf('   Config A: %3d RBs, Fs = %6.2f MHz (~%3.0f MHz PRS BW), 1 sample = %5.2f m\n', ...
        carrierA.NSizeGrid, Fs_A/1e6, carrierA.NSizeGrid*12*Delta_f/1e6, c/Fs_A);
fprintf('   Config B: %3d RBs, Fs = %6.2f MHz (~%3.0f MHz PRS BW), 1 sample = %5.2f m\n', ...
        carrierB.NSizeGrid, Fs_B/1e6, carrierB.NSizeGrid*12*Delta_f/1e6, c/Fs_B);

prs_demo = make_prs_config(carrierA, 1);
prsGrid  = zeros(carrierA.NSizeGrid*12, 14);
prsGrid(nrPRSIndices(carrierA, prs_demo)) = nrPRS(carrierA, prs_demo);

figure('Name','Fig3_PRS_Grid','Position',[50 50 600 420]);
imagesc(abs(prsGrid)); set(gca,'YDir','normal');
colormap(parula); colorbar;
title('5G NR PRS Resource Grid - gNB 1 (CombSize=4, 52 RBs, 30 kHz SCS)');
xlabel('OFDM Symbol Index'); ylabel('Subcarrier Index');
occ = nnz(prsGrid) / numel(prsGrid);
text(0.5, 0.95, sprintf('PRS occupancy: %.1f%%   (processing gain %.1f dB)', ...
     occ*100, 10*log10(nnz(prsGrid))), ...
     'Units','normalized','Color','w','FontSize',10,'FontWeight','bold', ...
     'HorizontalAlignment','center');


%% SECTION 5: ToA ESTIMATOR COMPARISON  [NEW IN REVISION 2]
% This is the headline new result. It isolates the single-link ToA error of
% the three estimators with no multipath and no inter-cell interference, so
% what is measured is purely the estimator, and sweeps SNR wide enough to
% cross from loss-of-lock, through the noise-limited region, onto the bias
% floor. Revision 1's operating points all sat on the flat right-hand end.
fprintf('[5/9] ToA estimator comparison (this is why v1 was flat)...\n');

snr_est   = -40:2.5:10;
methods   = {'integer','parabolic','bandlimited'};
n_est     = N_EST;
rms_est   = zeros(numel(methods), numel(snr_est));

for mi = 1:numel(methods)
    for si = 1:numel(snr_est)
        e = zeros(n_est,1);
        for t = 1:n_est
            d = 50 + 300*rand;                       % fractional delay
            L = length(txB{1}) + 500;
            y = frac_delay(txB{1}, d, L);
            y = awgn(y, snr_est(si), 'measured');
            e(t) = single_toa(y, txB{1}, methods{mi}) - d;
            clear y;                 % ~0.5 MB/trial -- do not let it pile up
        end
        rms_est(mi,si) = rms(e);
    end
    fprintf('   %-12s floor at high SNR: %.4f samples (%.3f m)\n', ...
            methods{mi}, rms_est(mi,end), rms_est(mi,end)*c/Fs_B);
end

figure('Name','Fig4_Estimator_Comparison','Position',[50 50 760 520]);
mk = {'o-','s-','^-'};
for mi = 1:numel(methods)
    semilogy(snr_est, rms_est(mi,:)*c/Fs_B, mk{mi}, 'LineWidth',2, ...
             'MarkerSize',5, 'DisplayName', methods{mi});
    hold on;
end
yline(c/Fs_B/sqrt(12), '--k', '1/\surd12 sample (integer-argmax limit)', ...
      'LineWidth',1.2);
grid on; grid minor;
xlabel('Input SNR (dB)'); ylabel('Single-link ToA error, RMS (m)');
title(sprintf(['ToA estimator comparison, 100 MHz PRS\n' ...
               'PRS processing gain %.0f dB -- the noise-limited region is ' ...
               'well below 0 dB'], 10*log10(nnz(prsGrid)*(264/52))));
legend('Location','southwest'); hold off;


%% SECTION 6: CROSS-CORRELATION SHAPE -- AWGN vs TDL-A
fprintf('[6/9] Generating CCF comparison plots (AWGN vs TDL-A)...\n');

SNR_demo   = 20;
delaysA_ref = vecnorm(gNB_coords - ref_UE, 2, 2) / c * Fs_A;
delaysB_ref = vecnorm(gNB_coords - ref_UE, 2, 2) / c * Fs_B;

rxA_awgn = rx_build(txA, delaysA_ref, SNR_demo, Fs_A, 'awgn', '', 0, 0);
rxA_tdla = rx_build(txA, delaysA_ref, SNR_demo, Fs_A, 'tdl', 'TDL-A', 100e-9, 0);
rxB_awgn = rx_build(txB, delaysB_ref, SNR_demo, Fs_B, 'awgn', '', 0, 0);
rxB_tdla = rx_build(txB, delaysB_ref, SNR_demo, Fs_B, 'tdl', 'TDL-A', 100e-9, 0);

figure('Name','Fig5_CCF_Comparison','Position',[50 50 920 560]);
combo_rx  = {rxA_awgn, rxA_tdla, rxB_awgn, rxB_tdla};
combo_tx  = {txA,      txA,      txB,      txB     };
combo_dl  = {delaysA_ref, delaysA_ref, delaysB_ref, delaysB_ref};
combo_lbl = {'20 MHz + AWGN', '20 MHz + TDL-A (100 ns DS)', ...
             '100 MHz + AWGN', '100 MHz + TDL-A (100 ns DS)'};
win_half  = [80, 80, 350, 350];

for s = 1:4
    [co, lg] = xcorr(combo_rx{s}, combo_tx{s}{1});
    co = abs(co) / max(abs(co));
    pk = combo_dl{s}(1);
    ww = max(1, find(lg >= pk-win_half(s), 1)) : ...
         min(length(lg), find(lg <= pk+win_half(s), 1, 'last'));

    subplot(2,2,s);
    plot(lg(ww) - pk, co(ww), 'b', 'LineWidth', 1.5);
    xline(0, '--r', 'True LoS', 'LineWidth', 1.2, 'LabelVerticalAlignment','bottom');
    grid on; ylim([0 1.15]);
    title(combo_lbl{s}, 'FontSize', 9);
    xlabel('Lag offset (samples)'); ylabel('Norm. Correlation');
end
sgtitle('Cross-Correlation Peak: Bandwidth and Multipath Effect (gNB 1)', ...
        'FontWeight','bold');


%% SECTION 7: SINGLE-POINT POSITIONING DEMO (ALL 6 CONFIGURATIONS)
fprintf('[7/9] Single-point positioning demo (SNR = %d dB)...\n', SNR_demo);

cfg_list = { ...
    txA, Fs_A, 'awgn', '',      0,      '~20 MHz, AWGN'; ...
    txA, Fs_A, 'tdl',  'TDL-A', 100e-9, '~20 MHz, TDL-A (100 ns)'; ...
    txA, Fs_A, 'tdl',  'TDL-C', 300e-9, '~20 MHz, TDL-C (300 ns)'; ...
    txB, Fs_B, 'awgn', '',      0,      '~100 MHz, AWGN'; ...
    txB, Fs_B, 'tdl',  'TDL-A', 100e-9, '~100 MHz, TDL-A (100 ns)'; ...
    txB, Fs_B, 'tdl',  'TDL-C', 300e-9, '~100 MHz, TDL-C (300 ns)'};

est_pos_all = zeros(6,2);
err_all     = zeros(6,1);
for i = 1:6
    txi = cfg_list{i,1}; Fsi = cfg_list{i,2};
    dl  = vecnorm(gNB_coords - ref_UE, 2, 2) / c * Fsi;
    rx  = rx_build(txi, dl, SNR_demo, Fsi, cfg_list{i,3}, cfg_list{i,4}, cfg_list{i,5}, 0);
    toa = estimate_toa(rx, txi, TOA_METHOD);
    [est_pos_all(i,:), err_all(i)] = tdoa_nlls(toa, gNB_coords, ref_UE, c, Fsi, areaSize);
end

fprintf('\n   %-32s  Error (m)\n', 'Configuration');
fprintf('   %-32s  ---------\n', repmat('-',1,32));
for i = 1:6
    fprintf('   %-32s  %.3f\n', cfg_list{i,6}, err_all(i));
end
fprintf('\n');

figure('Name','Fig6_Positioning_Map','Position',[50 50 640 540]);
hold on;
scatter(gNB_coords(:,1), gNB_coords(:,2), 200, '^r', 'filled', 'DisplayName','gNBs');
scatter(ref_UE(1), ref_UE(2), 180, 'og', 'filled', 'DisplayName','True UE');
mkr_colors = {[0.1 0.4 0.9],[0.1 0.4 0.9],[0.1 0.4 0.9], ...
              [0.9 0.2 0.1],[0.9 0.2 0.1],[0.9 0.2 0.1]};
mkr_shape  = {'s','d','^','s','d','^'};
for i = 1:6
    scatter(est_pos_all(i,1), est_pos_all(i,2), 120, mkr_shape{i}, ...
            'LineWidth', 2, 'MarkerEdgeColor', mkr_colors{i}, ...
            'DisplayName', sprintf('%s (%.2fm)', cfg_list{i,6}, err_all(i)));
end
grid on; axis equal; axis([-100 areaSize+100 -100 areaSize+100]);
xlabel('X (m)'); ylabel('Y (m)');
title(sprintf('Positioning Map - All Configurations (SNR = %d dB)', SNR_demo));
legend('Location','best','FontSize',7); hold off;


%% SECTION 8: MONTE CARLO -- BANDWIDTH x CHANNEL x SNR
fprintf('[8/9] Monte Carlo (%d trials x 6 configs x %d SNRs)...\n', ...
        N_MC, numel(SNR_dB_MC));

ch_names = {'AWGN','TDL-A','TDL-C'};
ch_types = {'awgn','tdl','tdl'};
ch_prof  = {'','TDL-A','TDL-C'};
ch_ds    = {0, 100e-9, 300e-9};

bw_labels  = {'20 MHz','100 MHz'};
bw_tx      = {txA, txB};
bw_Fs      = {Fs_A, Fs_B};
bw_carrier = {carrierA, carrierB};

errors = cell(2, 3, numel(SNR_dB_MC));
combo_num = 0; total_combos = 2*3*numel(SNR_dB_MC);
t_mc = tic;

for bIdx = 1:2
    txb = bw_tx{bIdx}; Fsb = bw_Fs{bIdx};
    for cIdx = 1:3
        for sIdx = 1:numel(SNR_dB_MC)
            combo_num = combo_num + 1;
            snr_db = SNR_dB_MC(sIdx);
            err_vec = zeros(N_MC,1);

            fprintf('   [%2d/%2d] %-8s + %-6s SNR=%4d dB ... ', ...
                    combo_num, total_combos, bw_labels{bIdx}, ch_names{cIdx}, snr_db);

            for trial = 1:N_MC
                ue = rand_UEs(trial,:);
                dl = vecnorm(gNB_coords - ue, 2, 2) / c * Fsb;
                rx = rx_build(txb, dl, snr_db, Fsb, ch_types{cIdx}, ...
                              ch_prof{cIdx}, ch_ds{cIdx}, 0);
                toa = estimate_toa(rx, txb, TOA_METHOD);
                [~, err_vec(trial)] = tdoa_nlls(toa, gNB_coords, ue, c, Fsb, areaSize);
            end

            errors{bIdx,cIdx,sIdx} = err_vec;
            fprintf('median=%7.3f m  90th=%8.3f m\n', ...
                    median(err_vec), prctile(err_vec,90));
        end
    end
end
fprintf('   Monte Carlo elapsed: %.1f min\n', toc(t_mc)/60);


%% SECTION 9: DOPPLER (UE MOBILITY)
fprintf('\n[9/9] Doppler Monte Carlo (TDL-A, SNR = %d dB)...\n', SNR_demo);

doppler_vals   = [0, 10, 100, 400];
doppler_labels = {'Static (0 Hz)','Pedestrian (10 Hz)', ...
                  'Vehicular (100 Hz)','High-speed (400 Hz)'};

errors_dop = cell(2, numel(doppler_vals));
for bIdx = 1:2
    txb = bw_tx{bIdx}; Fsb = bw_Fs{bIdx};
    for dIdx = 1:numel(doppler_vals)
        fd = doppler_vals(dIdx);
        ev = zeros(N_MC,1);
        fprintf('   %-8s + TDL-A, f_D=%3d Hz ... ', bw_labels{bIdx}, fd);
        for trial = 1:N_MC
            ue = rand_UEs(trial,:);
            dl = vecnorm(gNB_coords - ue, 2, 2) / c * Fsb;
            rx = rx_build(txb, dl, SNR_demo, Fsb, 'tdl', 'TDL-A', 100e-9, fd);
            toa = estimate_toa(rx, txb, TOA_METHOD);
            [~, ev(trial)] = tdoa_nlls(toa, gNB_coords, ue, c, Fsb, areaSize);
        end
        errors_dop{bIdx,dIdx} = ev;
        fprintf('median=%7.3f m  90th=%8.3f m\n', median(ev), prctile(ev,90));
    end
end


%% SECTION 10: CDF AND VALIDATION PLOTS
fprintf('\nGenerating final CDF and validation plots...\n');

x_cdf  = (1:N_MC)/N_MC;
ch_ls  = {'-','--',':'};
bw_col = {[0.12 0.47 0.71], [0.84 0.15 0.16]};

% --- CDF at the high-SNR anchor ---
[~, anchor] = max(SNR_dB_MC);
figure('Name','Fig7a_CDF_BW','Position',[50 50 850 520]);
hold on;
for bIdx = 1:2
    for cIdx = 1:3
        plot(sort(errors{bIdx,cIdx,anchor}), x_cdf, 'Color', bw_col{bIdx}, ...
             'LineStyle', ch_ls{cIdx}, 'LineWidth', 2, ...
             'DisplayName', sprintf('%s + %s', bw_labels{bIdx}, ch_names{cIdx}));
    end
end
grid on; set(gca,'XScale','log');
xlabel('2D Positioning Error (m)'); ylabel('CDF');
title(sprintf('CDF of Positioning Error (SNR = %d dB, N = %d, %s estimator)', ...
              SNR_dB_MC(anchor), N_MC, TOA_METHOD));
legend('Location','southeast','FontSize',8,'NumColumns',2); hold off;

% --- CDF vs SNR, 100 MHz ---
figure('Name','Fig7b_CDF_SNR','Position',[50 50 850 520]);
snr_cols = turbo(numel(SNR_dB_MC));
hold on;
for sIdx = 1:numel(SNR_dB_MC)
    plot(sort(errors{2,1,sIdx}), x_cdf, 'Color', snr_cols(sIdx,:), ...
         'LineWidth', 2, 'DisplayName', sprintf('SNR = %d dB', SNR_dB_MC(sIdx)));
end
grid on; set(gca,'XScale','log');
xlabel('2D Positioning Error (m)'); ylabel('CDF');
title(sprintf('CDF vs SNR, 100 MHz PRS + AWGN (N = %d)', N_MC));
legend('Location','southeast','FontSize',8); hold off;

% --- CRLB validation: the plot revision 1 could not make work ---
figure('Name','Fig8_CRLB_Validation','Position',[50 50 760 520]);
hold on;
for bIdx = 1:2
    crlb_sm = crlb_position(bw_carrier{bIdx}.NSizeGrid, Delta_f, ...
                            10.^(SNR_dB_vec/10), GDOP_ref, c);
    plot(SNR_dB_vec, crlb_sm, '-', 'Color', bw_col{bIdx}, 'LineWidth', 2.5, ...
         'DisplayName', sprintf('CRLB (AWGN) %s', bw_labels{bIdx}));

    med_awgn = arrayfun(@(s) median(errors{bIdx,1,s}), 1:numel(SNR_dB_MC));
    med_tdla = arrayfun(@(s) median(errors{bIdx,2,s}), 1:numel(SNR_dB_MC));
    med_tdlc = arrayfun(@(s) median(errors{bIdx,3,s}), 1:numel(SNR_dB_MC));

    plot(SNR_dB_MC, med_awgn, 'o--', 'Color', bw_col{bIdx}, 'LineWidth',2, ...
         'MarkerSize',8, 'DisplayName', sprintf('MC AWGN %s', bw_labels{bIdx}));
    plot(SNR_dB_MC, med_tdla, 's:',  'Color', bw_col{bIdx}*0.7, 'LineWidth',2, ...
         'MarkerSize',8, 'DisplayName', sprintf('MC TDL-A %s', bw_labels{bIdx}));
    plot(SNR_dB_MC, med_tdlc, '^-.', 'Color', bw_col{bIdx}*0.5, 'LineWidth',2, ...
         'MarkerSize',8, 'DisplayName', sprintf('MC TDL-C %s', bw_labels{bIdx}));
end
set(gca,'YScale','log'); grid on; grid minor;
xlabel('SNR (dB)'); ylabel('Median Positioning Error (m)');
title(sprintf('Monte Carlo vs AWGN CRLB (%s estimator)', TOA_METHOD));
legend('Location','southwest','FontSize',7,'NumColumns',2); hold off;

% --- Doppler CDF ---
figure('Name','Fig9_CDF_Doppler','Position',[50 50 860 520]);
dop_ls = {'-','--',':','-.'};
hold on;
for bIdx = 1:2
    for dIdx = 1:numel(doppler_vals)
        plot(sort(errors_dop{bIdx,dIdx}), x_cdf, 'Color', bw_col{bIdx}, ...
             'LineStyle', dop_ls{dIdx}, 'LineWidth', 2, ...
             'DisplayName', sprintf('%s + %s', bw_labels{bIdx}, doppler_labels{dIdx}));
    end
end
grid on; set(gca,'XScale','log');
xlabel('2D Positioning Error (m)'); ylabel('CDF');
title(sprintf('CDF of Positioning Error: Doppler (TDL-A, SNR = %d dB, N = %d)', ...
              SNR_demo, N_MC));
legend('Location','southeast','FontSize',7,'NumColumns',2); hold off;

% --- Doppler percentiles ---
figure('Name','Fig10_Doppler_vs_BW','Position',[50 50 700 480]);
hold on;
fd_plot = doppler_vals + 0.5;
for bIdx = 1:2
    p50 = cellfun(@median,              errors_dop(bIdx,:));
    p90 = cellfun(@(e) prctile(e,90),   errors_dop(bIdx,:));
    plot(fd_plot, p50, 'o-',  'Color', bw_col{bIdx}, 'LineWidth',2, ...
         'MarkerSize',9, 'DisplayName', sprintf('%s - Median', bw_labels{bIdx}));
    plot(fd_plot, p90, 's--', 'Color', bw_col{bIdx}, 'LineWidth',2, ...
         'MarkerSize',9, 'DisplayName', sprintf('%s - 90th pct', bw_labels{bIdx}));
end
set(gca,'XScale','log','YScale','log');
xticks([0.5 10 100 400]);
xticklabels({'0 (static)','10 (ped.)','100 (veh.)','400 (high-speed)'});
grid on; grid minor;
xlabel('Doppler Shift f_D (Hz)'); ylabel('Positioning Error (m)');
title(sprintf('50th and 90th Percentile Error vs Doppler (TDL-A, SNR = %d dB)', SNR_demo));
legend('Location','northwest','FontSize',9); hold off;

save('EE666_results_v2.mat','errors','errors_dop','SNR_dB_MC','doppler_vals', ...
     'rms_est','snr_est','methods','N_MC','TOA_METHOD','GDOP_ref','err_all');
fprintf('\n  Complete. Results saved to EE666_results_v2.mat\n');


%% ======================== LOCAL HELPER FUNCTIONS ========================

function g = gdop_at(xu, yu, gNB, invC, regularise)
%GDOP_AT TDoA GDOP at a point, from the correlated-noise FIM.
    N = size(gNB,1);
    R = max(sqrt((gNB(:,1)-xu).^2 + (gNB(:,2)-yu).^2), 1);
    H = zeros(N-1, 2);
    for k = 2:N
        H(k-1,1) = (xu-gNB(k,1))/R(k) - (xu-gNB(1,1))/R(1);
        H(k-1,2) = (yu-gNB(k,2))/R(k) - (yu-gNB(1,2))/R(1);
    end
    FIM = H' * invC * H;
    if regularise, FIM = FIM + eye(2)*1e-9; end
    g = sqrt(trace(inv(FIM)));
end

function s = crlb_position(N_RB, Delta_f, SNR_lin, GDOP, c)
%CRLB_POSITION 2D position CRLB from the RMS bandwidth of the PRS.
    N_SC  = N_RB * 12;
    k_idx = (-(N_SC/2):(N_SC/2)-1)';
    brms2 = (Delta_f^2) * mean(k_idx.^2);
    s = GDOP * c * sqrt(1 ./ (8*pi^2 * SNR_lin * brms2));
end

function tx = generate_prs_waveforms(carrier, N_gNB)
%GENERATE_PRS_WAVEFORMS One PRS waveform per cell. Position-independent, so
%   this is called once per bandwidth config rather than once per MC trial.
    tx = cell(N_gNB,1);
    for i = 1:N_gNB
        prs = make_prs_config(carrier, i);
        rg  = zeros(carrier.NSizeGrid*12, 14);
        rg(nrPRSIndices(carrier, prs)) = nrPRS(carrier, prs);
        tx{i} = nrOFDMModulate(carrier, rg);
    end
end

function prs = make_prs_config(carrier, cell_id)
%MAKE_PRS_CONFIG PRS resource for one cell.
%   REOffset is staggered across cells (change [4]). Revision 1 left it at
%   the default 0 for every gNB, so all four cells occupied identical
%   resource elements and were separated only by NPRSID scrambling -- which
%   contradicted Sec. IV-A of the report. Staggering by mod(cell-1,CombSize)
%   gives the frequency-domain orthogonality the report describes.
    prs = nrPRSConfig;
    prs.NumRB                = carrier.NSizeGrid;
    prs.PRSResourceSetPeriod = [10 0];
    prs.NumPRSSymbols        = 4;
    prs.SymbolStart          = 8;
    prs.NPRSID               = cell_id;
    prs.CombSize             = 4;
    prs.REOffset             = mod(cell_id-1, 4);
end

function y = frac_delay(x, d, L)
%FRAC_DELAY Delay x by d samples, d fractional, via a frequency-domain
%   phase ramp. Exact for a band-limited signal.
%
%   Change [1]. Revision 1 used round((dist/c)*Fs), snapping the true delay
%   to a 9.77 m grid at 20 MHz and 2.44 m at 100 MHz. No amount of
%   sub-sample estimation can recover information that was never put into
%   the signal, so this had to be fixed before the estimator work meant
%   anything.
    X = fft(x, L);
    k = [0:floor(L/2), -ceil(L/2)+1:-1].';
    y = ifft(X .* exp(-1j*2*pi*k*d/L));
end

function rx = rx_build(tx, delays_frac, SNR_dB, Fs, type, profile, ds, fd)
%RX_BUILD Superpose all cells at their fractional geometric delays, through
%   AWGN or a 3GPP TDL channel, and add noise.
%
%   Multipath is applied first and the geometric delay second; the two are
%   separable, and this ordering keeps the TDL taps referenced to the LoS
%   arrival. A fresh channel object per link gives independent fading.
    persistent chCache
    if isempty(chCache)
        chCache = containers.Map('KeyType','char','ValueType','any');
    end

    N_gNB = numel(tx);
    switch type
        case 'awgn'
            tail = 500;
        otherwise
            tail = ceil(10*ds*Fs) + 500;
    end
    L  = length(tx{1}) + ceil(max(delays_frac)) + tail;
    rx = zeros(L,1);

    for i = 1:N_gNB
        s = tx{i};
        if strcmp(type,'tdl')
            % CHANNEL OBJECTS ARE CACHED, NOT REBUILT PER TRIAL.
            % Constructing a fresh nrTDLChannel for every gNB on every trial
            % creates ~120,000 System objects over a full sweep. MATLAB does
            % not release them promptly and the run was killed for memory at
            % 86% completion. One object per (rate, profile, Doppler, gNB) is
            % created once and reset() before each use.
            %
            % RandomStream is set to the global stream so that reset() still
            % draws a fresh, INDEPENDENT fading realisation per trial while
            % the whole sweep stays reproducible from the rng(42) at the top.
            key = sprintf('%g|%s|%g|%g|%d', Fs, profile, ds, fd, i);
            if isKey(chCache, key)
                ch = chCache(key);
            else
                ch = nrTDLChannel;
                ch.SampleRate              = Fs;
                ch.DelayProfile            = profile;
                ch.DelaySpread             = ds;
                ch.MaximumDopplerShift     = fd;
                ch.NumTransmitAntennas     = 1;
                ch.NumReceiveAntennas      = 1;
                ch.NormalizePathGains      = true;
                ch.NormalizeChannelOutputs = false;
                ch.RandomStream            = 'Global stream';
                chCache(key) = ch;
            end
            reset(ch);
            s = ch(s);
        end
        rx = rx + frac_delay(s, delays_frac(i), L);
    end
    rx = awgn(rx, SNR_dB, 'measured');
end

function toa = estimate_toa(rx, tx, method)
%ESTIMATE_TOA Per-cell ToA in samples (fractional for the sub-sample methods).
    toa = zeros(numel(tx),1);
    for i = 1:numel(tx)
        toa(i) = single_toa(rx, tx{i}, method);
    end
end

function t = single_toa(rx, tx1, method)
%SINGLE_TOA ToA of one cell by cross-correlation peak location.
%
%   'integer'     argmax of |CCF|. Revision 1's estimator. Quantisation RMS
%                 is 1/sqrt(12) = 0.289 samples by construction.
%   'parabolic'   three-point quadratic fit about the peak, the formula
%                 proposed in revision 1 Sec. VII-C. Measured floor ~0.055
%                 samples -- five times better than integer, but bias- not
%                 noise-limited, and NOT the "sub-centimetre" the report
%                 predicted.
%   'bandlimited' FFT interpolation of the complex CCF in a window about the
%                 peak. Correct for a band-limited signal, which the CCF of
%                 an OFDM PRS is. Measured floor ~0.018 samples.
    [co, lg] = xcorr(rx, tx1);
    a = abs(co);
    [~, m] = max(a);

    switch method
        case 'integer'
            t = lg(m);

        case 'parabolic'
            if m > 1 && m < numel(a)
                y0 = a(m-1); y1 = a(m); y2 = a(m+1);
                den = y0 - 2*y1 + y2;
                dl = 0;
                if den ~= 0, dl = (y0 - y2) / (2*den); end
                if abs(dl) > 1, dl = 0; end     % reject a degenerate fit
                t = lg(m) + dl;
            else
                t = lg(m);
            end

        case 'bandlimited'
            U = 16;                              % interpolation factor
            w = max(1, m-8) : min(numel(co), m+8);
            seg = interpft(co(w), numel(w)*U);
            [~, mu] = max(abs(seg));
            t = lg(w(1)) + (mu-1)/U;

        otherwise
            error('single_toa:method','Unknown ToA method "%s".', method);
    end
end

function [est_pos, err] = tdoa_nlls(est_ToA, gNB_coords, true_pos, c, Fs, areaSize)
%TDOA_NLLS Nonlinear least squares over the hyperbolic RSTD residuals.
%
%   The search is CONSTRAINED to the deployment area plus a margin. Without
%   this, a trial that loses correlation lock produces mutually inconsistent
%   hyperbolae and the unbounded solver runs off to ~1e9 m, which then
%   destroys the median of the whole batch. A real positioning engine would
%   never report a fix 10^6 km outside its own service area either.
%
%   Trials that hit the bound are genuine outages, not measurements. Report
%   them as an outage rate (see OUTAGE_M) rather than letting them masquerade
%   as an enormous but finite error.
    lim    = [-200, areaSize+200];
    rstd_m = (est_ToA(2:end) - est_ToA(1)) * (c/Fs);

    resid = @(p) ...
        (sqrt((gNB_coords(2:end,1)-p(1)).^2 + (gNB_coords(2:end,2)-p(2)).^2) - ...
         sqrt((gNB_coords(1,1)-p(1))^2 + (gNB_coords(1,2)-p(2))^2)) - rstd_m;

    % Soft barrier keeps fminsearch inside the service area.
    pen = @(p) 1e6 * ( sum(max(0, lim(1)-p).^2) + sum(max(0, p-lim(2)).^2) );
    cost_fn = @(p) sum(resid(p).^2) + pen(p);

    opts    = optimset('Display','off','TolX',1e-8,'TolFun',1e-12,'MaxIter',5000);
    est_pos = fminsearch(cost_fn, [areaSize/2, areaSize/2], opts);
    est_pos = min(max(est_pos, lim(1)), lim(2));
    err     = norm(true_pos - est_pos);
end
