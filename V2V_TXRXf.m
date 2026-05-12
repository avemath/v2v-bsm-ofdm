% Avery Matherne and Jason Phan
%
% V2V_TXRX.m
% TX and RX both
% TX = usb:0, RX = usb:1
%
% Refs:
%   EE4003 labs:
%     Lab 3 (PSK_FULLCode_lab3.m, HammingCode_lab3.m): bit packing, FEC, RRC
%     Lab 9 (lab9_bpsk_pluto_loopback_pre_train_pll.m): pluto setup, MF,
%       full sync pipeline (PRE detection, tau-sweep, CFO est, TRAIN fit,
%       polarity, PLL, best-of-N decode)
%   Other:
%     OFDM via MATLAB Comm Toolbox
%     802.11p OFDM params: IEEE 802.11-2016 Clause 17
%     K=7 r=1/2 conv code: standard 802.11
%     BSM scaled-int encoding: SAE J2735 style
%     Channel est: SDR4Engineers Ch.10.7, MMSE EQ: standard form (book uses ZF)

clear; clc;

%% ----------------------- PARAMETERS -----------------------

USE_HARDWARE      = false;     % flip for pluto impl.
SNR_dB            = 30;        % sim only
FRAMES_TO_CAPTURE = 40;
CORR_THRESH       = 2;         % min peak/median ratio for PRE detection
NOISE_EST_LEN     = 200;
CFO_MAX_PLAUSIBLE = 200e3;     % Hz, skip if above this
TTC_WARN_THRESH   = 5.0;       % s, brake warning if closing TTC below this

% pluto config
% use 2.9 GHz instead of 5.9 GHz bc of pluto limits
TX_RADIO_ID       = 'usb:0';
RX_RADIO_ID       = 'usb:1';
CENTER_FREQUENCY  = 2.9e9;
TX_GAIN           = -10;
RX_GAIN           = 40;
RX_GAIN_SOURCE    = 'Manual';
RX_SAMPLES_PER_FRAME = 8192;
RX_OUTPUT_DTYPE   = 'double';

%% ----------------------- TX: BUILD PACKETS -----------------------

%% TX-1. BSM payloads (10 bytes = 80 bits each, two vehicles)
% Car 1: trailing, 30 m/s, no accel
% Car 2: lead, 25 m/s, decel -1 m/s^2, 20 m ahead
bsm1.vehicleID = 1; bsm1.velocity = 30; bsm1.accel = 0;  bsm1.posX = 100; bsm1.posY = 0;
bsm2.vehicleID = 2; bsm2.velocity = 25; bsm2.accel = -1; bsm2.posX = 120; bsm2.posY = 0;

PayloadBits1 = encodeBSMtoBits(bsm1);
PayloadBits2 = encodeBSMtoBits(bsm2);

fprintf('--- TRANSMITTING TWO BSMs ---\n');
fprintf('  Car 1: ID=%d  v=%.1f m/s  a=%+.1f m/s^2  pos=(%.0f,%.0f) m\n', ...
    bsm1.vehicleID, bsm1.velocity, bsm1.accel, bsm1.posX, bsm1.posY);
fprintf('  Car 2: ID=%d  v=%.1f m/s  a=%+.1f m/s^2  pos=(%.0f,%.0f) m\n\n', ...
    bsm2.vehicleID, bsm2.velocity, bsm2.accel, bsm2.posX, bsm2.posY);

%% TX-2. RRC Preamble (PRE) -- BPSK preamble, rng(42), shared by both packets
preLen   = 63;
sps      = 16;
rolloff  = 0.35;
span     = 6;
prePad   = 3*sps;
CFO_true = 5000;               % pluto XO is ~0.85 ppm

rng(42);                       % seed 42
preBits   = rand(preLen,1) > 0.5;
preSym    = 1 - 2*double(preBits);
rrc       = rcosdesign(rolloff, span, sps, 'sqrt');
preUp_raw = upfirdn(preSym, rrc, sps, 1);
preUp     = [zeros(prePad,1); preUp_raw; zeros(prePad,1)];
preUpLen  = length(preUp);

%% TX-2b. BPSK TRAIN sequence (256 known symbols), shared by both packets
trainLen   = 256;
trainBits  = rand(trainLen,1) > 0.5;
trainSym   = 1 - 2*double(trainBits);
trainUp    = upfirdn(trainSym, rrc, sps, 1);
trainUpLen = length(trainUp);

%% TX-3. OFDM + FEC parameters (802.11p PHY)
% OFDM via MATLAB Comm Toolbox, FEC concept from lab3
FFTLength            = 64;
NumGuardBandCarriers = [6; 5];
NumDataCarriers      = 48;
CyclicPrefixLength   = 16;            % 802.11p long GI = 1.6 us @ 10 MHz
PilotCarrierIndices  = [12;26;40;54]; % pilots @ -21,-7,+7,+21
Fs                   = 10e6;

trellis = poly2trellis(7,[171 133]);  % 802.11 K=7 r=1/2

EncodedBits1   = convenc(PayloadBits1, trellis);
EncodedBits2   = convenc(PayloadBits2, trellis);
NumCodedBits   = length(EncodedBits1);
NumOFDMSymbols = ceil(NumCodedBits / NumDataCarriers);
NumPadBits     = NumDataCarriers*NumOFDMSymbols - NumCodedBits;
PaddedBits1    = [EncodedBits1; zeros(NumPadBits,1)];
PaddedBits2    = [EncodedBits2; zeros(NumPadBits,1)];

mkMod = @(N) comm.OFDMModulator( ...
    'FFTLength',FFTLength, 'NumGuardBandCarriers',NumGuardBandCarriers, ...
    'InsertDCNull',true, 'PilotInputPort',true, ...
    'PilotCarrierIndices',PilotCarrierIndices, ...
    'CyclicPrefixLength',CyclicPrefixLength, 'NumSymbols',N);

%% TX-4. OFDM Preamble (seed 43, separate from PRE seed 42), shared by both packets
rng(43);
PreambleDataBits = randi([0 1], NumDataCarriers, 1);
preambleSym      = pskmod(PreambleDataBits, 2);
preMod           = mkMod(1);
Preamble         = preMod(preambleSym, ones(4,1));

%% TX-5. Data OFDM (one per vehicle)
DataSymbols1 = pskmod(PaddedBits1, 2);
DataSymbols2 = pskmod(PaddedBits2, 2);
dataMod1     = mkMod(NumOFDMSymbols);
dataMod2     = mkMod(NumOFDMSymbols);
DataOFDM1    = dataMod1( ...
    reshape(DataSymbols1, NumDataCarriers, NumOFDMSymbols), ones(4,NumOFDMSymbols));
DataOFDM2    = dataMod2( ...
    reshape(DataSymbols2, NumDataCarriers, NumOFDMSymbols), ones(4,NumOFDMSymbols));

%% TX-6. Assemble packets, concatenate Car 2 then Car 1, apply CFO, normalize
% [RRC PRE | TRAIN | OFDM Preamble | OFDM Data | guard]  per packet
guardSamps = round(Fs*50e-6);
Packet1    = [preUp; trainUp; Preamble; DataOFDM1; zeros(guardSamps,1)];
Packet2    = [preUp; trainUp; Preamble; DataOFDM2; zeros(guardSamps,1)];
TxWaveform = [Packet2; Packet1];
N_pkt      = length(Packet1);    % single-packet length
N_total    = length(TxWaveform);

cfoRot     = exp(1j*2*pi*(CFO_true/Fs)*(0:N_total-1).');
TxWaveform = TxWaveform .* cfoRot;
Packet_tx  = TxWaveform / max(abs(TxWaveform)) * 0.9;

fprintf('PRE        : %d samples (%.1f us)\n', length(preUp), length(preUp)/Fs*1e6);
fprintf('Each packet: %d bits -> coded %d bits -> %d OFDM syms (%d samples)\n', ...
    length(PayloadBits1), NumCodedBits, NumOFDMSymbols, N_pkt);
fprintf('TX waveform: %d samples total (%.1f us, two packets)\n\n', ...
    length(Packet_tx), length(Packet_tx)/Fs*1e6);

%% ----------------------- TX: TRANSMIT -----------------------

if USE_HARDWARE
    tx = sdrtx('Pluto','RadioID',TX_RADIO_ID,'CenterFrequency',CENTER_FREQUENCY, ...
               'BasebandSampleRate',Fs,'Gain',TX_GAIN);
    transmitRepeat(tx, Packet_tx);
    fprintf('TX running. Capturing RX now...\n\n');
end

%% ----------------------- RX: CAPTURE -----------------------

mkDemod = @(N) comm.OFDMDemodulator( ...
    'FFTLength',FFTLength, 'NumGuardBandCarriers',NumGuardBandCarriers, ...
    'RemoveDCCarrier',true, 'PilotOutputPort',true, ...
    'PilotCarrierIndices',PilotCarrierIndices, ...
    'CyclicPrefixLength',CyclicPrefixLength, 'NumSymbols',N);

%% RX-1. Rebuild reference preamble (must match TX-4)
rng(43);
preambleSym_ref = pskmod(randi([0 1], NumDataCarriers, 1), 2);
preMod_rx       = mkMod(1);
Preamble_ref    = preMod_rx(preambleSym_ref, ones(4,1));

%% RX-2. Capture or simulate
if USE_HARDWARE
    rx = sdrrx('Pluto','RadioID',RX_RADIO_ID,'CenterFrequency',CENTER_FREQUENCY, ...
        'BasebandSampleRate',Fs,'SamplesPerFrame',RX_SAMPLES_PER_FRAME, ...
        'GainSource',RX_GAIN_SOURCE,'Gain',RX_GAIN,'OutputDataType',RX_OUTPUT_DTYPE);

    % flush so AGC settles
    FLUSH_FRAMES = 5;
    for fi = 1:FLUSH_FRAMES
        rx();
    end

    fprintf('Capturing %d frames (%d samples, %.1f ms)...\n', ...
        FRAMES_TO_CAPTURE, FRAMES_TO_CAPTURE*RX_SAMPLES_PER_FRAME, ...
        FRAMES_TO_CAPTURE*RX_SAMPLES_PER_FRAME/Fs*1e3);
    rxBuf = [];
    for fi = 1:FRAMES_TO_CAPTURE
        rxBuf = [rxBuf; rx()];  %#ok<AGROW>
    end
    release(rx);
    rxSignal = rxBuf;
else
    % sim: replicate to mimic transmitRepeat, then channel + noise
    nReps   = 25;
    TxClean = repmat([Packet2; Packet1], nReps, 1);
    TxClean = TxClean / max(abs(TxClean)) * 0.9;
    N2      = length(TxClean);
    Tx_cfo  = TxClean .* exp(1j*2*pi*(CFO_true/Fs)*(0:N2-1).');

    % sim channel: moderate vehicular fading, less than highway V2V
    % (~1180 Hz Doppler) but more than lab env.
    chan = comm.RicianChannel('SampleRate',Fs,'KFactor',3, ...
        'MaximumDopplerShift',500,'PathDelays',[0 1e-7 2.5e-7], ...
        'AveragePathGains',[0 -5 -9], ...
        'RandomStream','mt19937ar with seed','Seed',99);
    rxFaded     = chan(Tx_cfo);
    rxWithNoise = awgn(rxFaded, SNR_dB, 'measured');

    noiseLen    = 500;
    sigPower    = mean(abs(Tx_cfo).^2);
    noiseSigma  = sqrt(sigPower / (2 * 10^(SNR_dB/10)));
    noisePrefix = noiseSigma * (randn(noiseLen,1) + 1j*randn(noiseLen,1));
    rxSignal    = [noisePrefix; rxWithNoise];
end

%% ----------------------- RX: SYNC + EQUALIZE + DECODE -----------------------
% receiver chain adapted from lab9

%% RX-3a. Pre-CFO matched filter
preambleLen   = length(Preamble_ref);
samplesNeeded = (FFTLength + CyclicPrefixLength) * NumOFDMSymbols;
rxMF = conv(rxSignal, rrc, 'same');

%% RX-3b. PRE cross-correlation for packet detection
M_pre   = length(preUp_raw);
rxSignal_raw = rxSignal;
corrPRE = conv(rxSignal_raw, flipud(conj(preUp_raw)), 'same');
noiseFloor = median(abs(corrPRE));
minTail    = preUpLen + trainUpLen + preambleLen + samplesNeeded + 100;

absCorrPRE = abs(corrPRE);
peakThresh = noiseFloor * CORR_THRESH;
[pks, locs] = findpeaks(absCorrPRE, 'MinPeakHeight', peakThresh, ...
    'MinPeakDistance', round(N_pkt * 0.5));

if isempty(locs)
    [~, locs] = max(absCorrPRE);
    pks = absCorrPRE(locs);
    warning('No peaks above threshold, using global max.');
end

% need tail room for full packet
hasRoom    = (locs + minTail) <= length(rxSignal_raw);
validLocs  = locs(hasRoom);
validPks   = pks(hasRoom);
if isempty(validLocs)
    [~, bi] = max(pks);
    validLocs = locs(bi);
    validPks  = pks(bi);
end

% try strongest peak first
[~, sortIdx] = sort(validPks, 'descend');
validLocs = validLocs(sortIdx);
validPks  = validPks(sortIdx);
isPredicted = false(size(validLocs));   % all findpeaks results are real

% augment w/ predictions: TX is [pkt2; pkt1] repeating, so peaks are N_pkt 
% apart, add predictions in case findpeaks misses any
extraLocs = [];
extraPks  = [];
for ip = 1:length(validLocs)
    p = validLocs(ip);
    for offset = [-N_pkt, N_pkt, -2*N_pkt, 2*N_pkt]
        c = p + offset;
        if c >= 1 && (c + minTail) <= length(rxSignal_raw)
            existing = [validLocs; extraLocs];
            if isempty(existing) || all(abs(existing - c) > round(N_pkt*0.25))
                extraLocs(end+1,1) = c;             %#ok<AGROW>
                extraPks(end+1,1)  = absCorrPRE(c); %#ok<AGROW>
            end
        end
    end
end
if ~isempty(extraLocs)
    fprintf('Added %d predicted peak positions.\n', length(extraLocs));
    validLocs   = [validLocs; extraLocs];
    validPks    = [validPks; extraPks];
    isPredicted = [isPredicted; true(length(extraLocs),1)];
end

fprintf('Found %d correlation peaks, trying %d total (%d real + %d predicted).\n', ...
    length(locs), length(validLocs), sum(~isPredicted), sum(isPredicted));

% reference bits per vehicle, for BER and assignment
refBits1 = PayloadBits1;
refBits2 = PayloadBits2;

%% Multi-packet decode loop, results keyed by best-match vehicle
results = struct('v1', [], 'v2', []);

for ci = 1:length(validLocs)
    ipk_try = validLocs(ci);
    corrSNR_try = validPks(ci) / noiseFloor;
    isPred_try = isPredicted(ci);
    if isPred_try
        fprintf('\n--- Trying peak %d/%d  (sample %d, corrSNR=%.1f, predicted) ---\n', ...
            ci, length(validLocs), ipk_try, corrSNR_try);
    else
        fprintf('\n--- Trying peak %d/%d  (sample %d, corrSNR=%.1f) ---\n', ...
            ci, length(validLocs), ipk_try, corrSNR_try);
    end

    if ~isPred_try && corrSNR_try < CORR_THRESH
        fprintf('  Skipping, below threshold.\n');
        continue;
    end

    % fresh copy each attempt
    rxSig = rxSignal_raw;
    rxMF_try = conv(rxSig, rrc, 'same');

    % tau-sweep over 0..sps-1 offsets
    tauList_try  = 0:(sps-1);
    scoreTau_try = zeros(1, sps);
    for ti = 1:length(tauList_try)
        tau = tauList_try(ti);
        idx = round((ipk_try - floor(M_pre/2) + span/2*sps + tau) + (0:preLen-1)*sps);
        idx = idx(idx >= 1 & idx <= length(rxMF_try));
        if isempty(idx), continue; end
        xs = rxMF_try(idx);
        scoreTau_try(ti) = abs(sum(xs .* conj(preSym(1:numel(xs)))));
    end
    [~, bestTauIdx_try] = max(scoreTau_try);
    tauBest_try = tauList_try(bestTauIdx_try);

    preStart_try    = ipk_try - floor(M_pre/2);
    trainStart_try  = round(preStart_try + M_pre + prePad + tauBest_try);
    packetStart_try = trainStart_try + trainUpLen;
    packetStart_try = max(1, min(packetStart_try, length(rxSig) - preambleLen - samplesNeeded));

    % CFO from PRE phase increments
    idxPRE_try  = round((ipk_try - floor(M_pre/2) + span/2*sps + tauBest_try) + (0:preLen-1)*sps);
    idxPRE_try  = idxPRE_try(idxPRE_try >= 1 & idxPRE_try <= length(rxMF_try));
    preTake_try = rxMF_try(idxPRE_try);

    preTakeDM = preTake_try .* conj(preSym(1:numel(preTake_try)));
    r_try     = preTakeDM(2:end) .* conj(preTakeDM(1:end-1));
    phi_try   = angle(sum(r_try));
    CFO_try   = (phi_try/(2*pi)) * (Fs/sps);

    if USE_HARDWARE && abs(CFO_try) > CFO_MAX_PLAUSIBLE
        fprintf('  CFO %.0f Hz implausible, skipping.\n', CFO_try);
        continue;
    end

    n_rx_try = (0:length(rxSig)-1).';
    rxSig    = rxSig .* exp(-1j*2*pi*(CFO_try/Fs)*n_rx_try);

    % TRAIN polyfit (linear phase + polarity)
    rxMF2_try      = conv(rxSig, rrc, 'same');
    trainSymI0_try = trainStart_try + round((span/2)*sps);
    idxTRAIN_try   = trainSymI0_try + (0:trainLen-1)*sps;
    idxTRAIN_try   = idxTRAIN_try(idxTRAIN_try >= 1 & idxTRAIN_try <= length(rxMF2_try));
    trainTake_try  = rxMF2_try(idxTRAIN_try);
    trainTake_try  = trainTake_try / rms(trainTake_try);
    refTrain_try   = trainSym(1:numel(trainTake_try));

    H_train_try = trainTake_try .* conj(refTrain_try);
    ang_t_try   = unwrap(angle(H_train_try));
    k_t_try     = (0:numel(ang_t_try)-1).';
    p_t_try     = polyfit(k_t_try, ang_t_try, 1);
    a_hat_try   = p_t_try(1);
    phi0_try    = p_t_try(2);

    trainBefore_try = trainTake_try;
    trainRot_try = trainTake_try .* exp(-1j*(phi0_try + a_hat_try*k_t_try));
    pol_try      = sign(real(sum(trainRot_try .* conj(refTrain_try))));
    if pol_try == 0, pol_try = 1; end
    trainAfter_try = trainRot_try * pol_try;

    a_samp_try = a_hat_try / sps;
    nTail_try  = length(rxSig) - trainStart_try + 1;
    mTail_try  = (0:nTail_try-1).';
    rotTail_try = pol_try * exp(-1j*(phi0_try + a_samp_try*mTail_try));
    rxSig(trainStart_try:end) = rxSig(trainStart_try:end) .* rotTail_try;

    % PLL on PRE (DD BPSK Costas)
    rxMF3_try    = conv(rxSig, rrc, 'same');
    preTake3_try = rxMF3_try(idxPRE_try);
    preTake3_try = preTake3_try / rms(preTake3_try);
    beta_pll = 0.001; burn_pll = 8; phi_hat_try = 0;
    zPLL_try = zeros(size(preTake3_try));
    for kk = 1:length(preTake3_try)
        zPLL_try(kk) = preTake3_try(kk) * exp(-1j*phi_hat_try);
        if kk > burn_pll
            err_pll = sign(real(zPLL_try(kk))) * imag(zPLL_try(kk));
            if abs(err_pll) < 5e-3, err_pll = 0; end
            phi_hat_try = phi_hat_try + beta_pll*err_pll;
        end
    end
    preRefAlign_try = preSym(1:numel(zPLL_try));
    phiAlign_try    = angle(sum(zPLL_try .* conj(preRefAlign_try)));
    phi_hat_try     = phi_hat_try + phiAlign_try;
    zPLL_try        = zPLL_try .* exp(-1j*phiAlign_try);

    rxSig = rxSig .* exp(-1j * phi_hat_try);

    % slice and demod
    if packetStart_try < 1 || packetStart_try + preambleLen - 1 > length(rxSig)
        fprintf('  Timing out of range, skipping.\n');
        continue;
    end
    rxPreamble_try = rxSig(packetStart_try : packetStart_try + preambleLen - 1);
    rxData_try     = rxSig(packetStart_try + preambleLen : end);
    if length(rxData_try) < samplesNeeded
        rxData_try = [rxData_try; zeros(samplesNeeded - length(rxData_try), 1)]; %#ok<AGROW>
    end
    rxData_try = rxData_try(1:samplesNeeded);

    preDemod_try      = mkDemod(1);
    [rxPreSub_try, ~] = preDemod_try(rxPreamble_try);
    % channel estimate (Ch.10.7, received/known)
    H_preamble_try    = rxPreSub_try ./ preambleSym_ref;

    dataDemod_try          = mkDemod(NumOFDMSymbols);
    [rxSubcarriers_try, ~] = dataDemod_try(rxData_try);
    H_est_try              = repmat(H_preamble_try, 1, NumOFDMSymbols);

    noiseWinEnd_try   = max(1, preStart_try - 50);
    noiseWinStart_try = max(1, noiseWinEnd_try - NOISE_EST_LEN);
    if noiseWinEnd_try > noiseWinStart_try
        N0_try = mean(abs(rxSig(noiseWinStart_try:noiseWinEnd_try)).^2);
    else
        N0_try = mean(abs(rxSubcarriers_try(:)).^2) * 0.01;
    end

    % MMSE EQ
    rxEq_try = (conj(H_est_try) .* rxSubcarriers_try) ./ (abs(H_est_try).^2 + N0_try);

    rxCoded_try = pskdemod(rxEq_try(:), 2);
    rxCoded_try = rxCoded_try(1 : end - NumPadBits);
    rxBits_try  = vitdec(rxCoded_try, poly2trellis(7,[171 133]), 35, 'term', 'hard');
    rxBits_try  = rxBits_try(1:80);

    % assign packet to vehicle by lower BER vs each reference
    BER1_try = sum(rxBits_try ~= refBits1) / 80;
    BER2_try = sum(rxBits_try ~= refBits2) / 80;
    if BER1_try <= BER2_try
        bestRef_try = 1;  bestBER_try = BER1_try;
    else
        bestRef_try = 2;  bestBER_try = BER2_try;
    end
    fprintf('  BER vs Car1=%.4f  vs Car2=%.4f  (CFO=%.0f Hz, tau=%d, pol=%+d) -> Car %d\n', ...
        BER1_try, BER2_try, CFO_try, tauBest_try, pol_try, bestRef_try);

    % save state for this peak
    state.BER         = bestBER_try;
    state.corrSNR     = corrSNR_try;
    state.ipk         = ipk_try;
    state.tauBest     = tauBest_try;
    state.bestTauIdx  = bestTauIdx_try;
    state.scoreTau    = scoreTau_try;
    state.tauList     = tauList_try;
    state.preStart    = preStart_try;
    state.trainStart  = trainStart_try;
    state.packetStart = packetStart_try;
    state.CFO         = CFO_try;
    state.idxPRE      = idxPRE_try;
    state.phi0        = phi0_try;
    state.a_hat       = a_hat_try;
    state.pol         = pol_try;
    state.phi_hat     = phi_hat_try;
    state.trainBefore = trainBefore_try;
    state.trainAfter  = trainAfter_try;
    state.preTakePreCFO = preTake_try;
    state.preTake3    = preTake3_try;
    state.zPLL        = zPLL_try;
    state.H_preamble  = H_preamble_try;
    state.rxSubcarriers = rxSubcarriers_try;
    state.rxEq        = rxEq_try;
    state.rxBits      = rxBits_try;
    state.N0          = N0_try;

    fld = sprintf('v%d', bestRef_try);
    if isempty(results.(fld)) || bestBER_try < results.(fld).BER
        results.(fld) = state;
    end

    % stop if both vehicles already at BER=0
    if ~isempty(results.v1) && ~isempty(results.v2) && ...
       results.v1.BER == 0 && results.v2.BER == 0
        fprintf('\nBoth vehicles decoded with BER=0, stopping.\n');
        break;
    end
end

if isempty(results.v1) || isempty(results.v2)
    error('Failed to decode both vehicles. v1=%s, v2=%s', ...
        mat2str(~isempty(results.v1)), mat2str(~isempty(results.v2)));
end

fprintf('\n--------------------------------------------\n');
fprintf('  Car 1: peak %d, BER = %.4f, corrSNR = %.1f, CFO = %.0f Hz\n', ...
    results.v1.ipk, results.v1.BER, results.v1.corrSNR, results.v1.CFO);
fprintf('  Car 2: peak %d, BER = %.4f, corrSNR = %.1f, CFO = %.0f Hz\n', ...
    results.v2.ipk, results.v2.BER, results.v2.corrSNR, results.v2.CFO);
fprintf('--------------------------------------------\n');

%% RX-8. Decode BSM fields per vehicle, compute TTC, brake warning
bsm1_rx = decodeBSMfromBits(results.v1.rxBits);
bsm2_rx = decodeBSMfromBits(results.v2.rxBits);

fprintf('\n--- BSM: GROUND TRUTH vs DECODED ---\n');
fprintf('  %-12s  %14s   %14s\n', 'Field', 'Truth', 'Decoded');
fprintf('  %-12s  %14d   %14d\n',     'Car 1 ID',   bsm1.vehicleID, bsm1_rx.vehicleID);
fprintf('  %-12s  %14.2f   %14.2f\n', 'Car 1 vel',  bsm1.velocity,  bsm1_rx.velocity);
fprintf('  %-12s  %14.2f   %14.2f\n', 'Car 1 acc',  bsm1.accel,     bsm1_rx.accel);
fprintf('  %-12s  %14.1f   %14.1f\n', 'Car 1 posX', bsm1.posX,      bsm1_rx.posX);
fprintf('  %-12s  %14.1f   %14.1f\n', 'Car 1 posY', bsm1.posY,      bsm1_rx.posY);
fprintf('  %-12s  %14d   %14d\n',     'Car 2 ID',   bsm2.vehicleID, bsm2_rx.vehicleID);
fprintf('  %-12s  %14.2f   %14.2f\n', 'Car 2 vel',  bsm2.velocity,  bsm2_rx.velocity);
fprintf('  %-12s  %14.2f   %14.2f\n', 'Car 2 acc',  bsm2.accel,     bsm2_rx.accel);
fprintf('  %-12s  %14.1f   %14.1f\n', 'Car 2 posX', bsm2.posX,      bsm2_rx.posX);
fprintf('  %-12s  %14.1f   %14.1f\n', 'Car 2 posY', bsm2.posY,      bsm2_rx.posY);

% TTC: car 1 trails car 2, closing if v1 > v2
gap         = bsm2_rx.posX - bsm1_rx.posX;
closingRate = bsm1_rx.velocity - bsm2_rx.velocity;
if closingRate > 0 && gap > 0
    TTC = gap / closingRate;
else
    TTC = Inf;
end

fprintf('\n--- COLLISION AVOIDANCE ---\n');
fprintf('  Gap         : %.2f m\n', gap);
fprintf('  Closing rate: %.2f m/s\n', closingRate);
if isfinite(TTC)
    fprintf('  TTC         : %.2f s\n', TTC);
else
    fprintf('  TTC         : Inf (not closing)\n');
end
brakeWarn = isfinite(TTC) && TTC < TTC_WARN_THRESH;
if brakeWarn
    fprintf('  ** WARNING: Brake Recommended **\n');
else
    fprintf('  No warning (TTC >= %.1f s threshold).\n', TTC_WARN_THRESH);
end

%% ----------------------- STOP TX -----------------------

% RX decode done, safe to stop
if USE_HARDWARE
    release(tx);
    fprintf('\nTX stopped.\n');
end

%% ----------------------- PLOTS -----------------------

figDir = fullfile(pwd, 'V2V_figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

% Stage 1: TX spectrum
figure(1); clf; set(gcf,'Name','V2V TX Spectrum');
pwelch([Preamble; DataOFDM1], 256, [], 1024, Fs, 'centered');
title('V2V OFDM Packet Spectrum');
saveas(gcf, fullfile(figDir, 'fig1_tx_spectrum_baseband.png'));

figure(2); clf; set(gcf,'Name','V2V TX Spectrum (RF, PlutoSDR)');
[Pxx, f] = pwelch([Preamble; DataOFDM1], 256, [], 1024, Fs, 'centered');
plot((f + CENTER_FREQUENCY)/1e9, 10*log10(Pxx));
grid on;
xlabel('Frequency (GHz)');
ylabel('Power/Frequency (dB/Hz)');
title(sprintf('V2V OFDM Packet Spectrum @ Fc = %.3f GHz', CENTER_FREQUENCY/1e9));
saveas(gcf, fullfile(figDir, 'fig2_tx_spectrum_rf.png'));

if USE_HARDWARE
    % Stage 2: RX spectrum
    figure(3); clf; set(gcf,'Name','V2V RX Spectrum (PlutoSDR Capture)');
    [Pxx_rx, f_rx] = pwelch(rxSignal, 256, [], 1024, Fs, 'centered');
    plot((f_rx + CENTER_FREQUENCY)/1e9, 10*log10(Pxx_rx));
    grid on;
    xlabel('Frequency (GHz)');
    ylabel('Power/Frequency (dB/Hz)');
    title(sprintf('V2V RX Captured Spectrum @ Fc = %.3f GHz', CENTER_FREQUENCY/1e9));
    saveas(gcf, fullfile(figDir, 'fig3_rx_spectrum_rf.png'));
end

% Stage 3: PRE corr (multi-packet loop)
figure(4); clf; set(gcf,'Name','PRE Correlation - Multi-Packet Detection');
tCorr = (0:length(corrPRE)-1) / Fs * 1e6;
plot(tCorr, abs(corrPRE), 'b'); hold on;
plot(results.v2.ipk/Fs*1e6, abs(corrPRE(results.v2.ipk)), 'rv', ...
    'MarkerSize',12, 'MarkerFaceColor','r', ...
    'DisplayName', sprintf('Car 2 (BER=%.4f)', results.v2.BER));
plot(results.v1.ipk/Fs*1e6, abs(corrPRE(results.v1.ipk)), 'gv', ...
    'MarkerSize',12, 'MarkerFaceColor','g', ...
    'DisplayName', sprintf('Car 1 (BER=%.4f)', results.v1.BER));
xlabel('Time (us)'); ylabel('|Correlation|');
title('PRE Correlation - Multi-Packet Detection');
legend('show','Location','best'); grid on;
saveas(gcf, fullfile(figDir, 'fig4_pre_correlation.png'));

% Stage 4: tau-sweep
figure(5); clf; set(gcf,'Name','tau-sweep');
subplot(1,2,1);
stem(results.v2.tauList, results.v2.scoreTau, 'b', 'filled'); hold on;
stem(results.v2.tauBest, results.v2.scoreTau(results.v2.bestTauIdx), 'r', 'filled', 'MarkerSize',8);
xlabel('tau (samples)'); ylabel('Score');
title(sprintf('CAR 2 - Best tau = %d', results.v2.tauBest));
grid on; xticks(results.v2.tauList);
subplot(1,2,2);
stem(results.v1.tauList, results.v1.scoreTau, 'b', 'filled'); hold on;
stem(results.v1.tauBest, results.v1.scoreTau(results.v1.bestTauIdx), 'r', 'filled', 'MarkerSize',8);
xlabel('tau (samples)'); ylabel('Score');
title(sprintf('CAR 1 - Best tau = %d', results.v1.tauBest));
grid on; xticks(results.v1.tauList);
saveas(gcf, fullfile(figDir, 'fig5_tau_sweep.png'));

% Stage 5: PRE Before/After CFO
figure(6); clf; set(gcf,'Name','PRE: Before vs After CFO');
subplot(2,2,1);
plot(real(results.v2.preTakePreCFO), imag(results.v2.preTakePreCFO), '.', 'MarkerSize',10);
axis equal; grid on; xline(0,'k:'); yline(0,'k:');
xlabel('I'); ylabel('Q'); title('CAR 2 PRE Before CFO');
subplot(2,2,2);
plot(real(results.v1.preTakePreCFO), imag(results.v1.preTakePreCFO), '.', 'MarkerSize',10);
axis equal; grid on; xline(0,'k:'); yline(0,'k:');
xlabel('I'); ylabel('Q'); title('CAR 1 PRE Before CFO');
subplot(2,2,3);
plot(real(results.v2.preTake3), imag(results.v2.preTake3), '.', 'MarkerSize',10);
axis equal; grid on; xline(0,'k:'); yline(0,'k:');
xlabel('I'); ylabel('Q'); title(sprintf('CAR 2 PRE After CFO (est=%.0f Hz)', results.v2.CFO));
subplot(2,2,4);
plot(real(results.v1.preTake3), imag(results.v1.preTake3), '.', 'MarkerSize',10);
axis equal; grid on; xline(0,'k:'); yline(0,'k:');
xlabel('I'); ylabel('Q'); title(sprintf('CAR 1 PRE After CFO (est=%.0f Hz)', results.v1.CFO));
saveas(gcf, fullfile(figDir, 'fig6_pre_cfo.png'));

% Stage 6: TRAIN linear-phase + pol fit
figure(7); clf; set(gcf,'Name','TRAIN: Before vs After Linear-Phase + Polarity Fit');
subplot(2,2,1);
scatter(real(results.v2.trainBefore), imag(results.v2.trainBefore), 20, ...
    1:numel(results.v2.trainBefore), 'filled');
axis equal; grid on; xline(0,'k:'); yline(0,'k:');
xlabel('I'); ylabel('Q');
title(sprintf('CAR 2 TRAIN Before  (\\phi_0=%+.2f rad)', results.v2.phi0));
colormap(gca, parula); cb = colorbar; cb.Label.String = 'symbol index';
subplot(2,2,2);
scatter(real(results.v1.trainBefore), imag(results.v1.trainBefore), 20, ...
    1:numel(results.v1.trainBefore), 'filled');
axis equal; grid on; xline(0,'k:'); yline(0,'k:');
xlabel('I'); ylabel('Q');
title(sprintf('CAR 1 TRAIN Before  (\\phi_0=%+.2f rad)', results.v1.phi0));
colormap(gca, parula); cb = colorbar; cb.Label.String = 'symbol index';
subplot(2,2,3);
scatter(real(results.v2.trainAfter), imag(results.v2.trainAfter), 20, ...
    1:numel(results.v2.trainAfter), 'filled');
axis equal; grid on; xline(0,'k:'); yline(0,'k:');
xlim([-2 2]); ylim([-2 2]);
xlabel('I'); ylabel('Q');
title(sprintf('CAR 2 TRAIN After  (polarity=%+d)', results.v2.pol));
colormap(gca, parula); cb = colorbar; cb.Label.String = 'symbol index';
subplot(2,2,4);
scatter(real(results.v1.trainAfter), imag(results.v1.trainAfter), 20, ...
    1:numel(results.v1.trainAfter), 'filled');
axis equal; grid on; xline(0,'k:'); yline(0,'k:');
xlim([-2 2]); ylim([-2 2]);
xlabel('I'); ylabel('Q');
title(sprintf('CAR 1 TRAIN After  (polarity=%+d)', results.v1.pol));
colormap(gca, parula); cb = colorbar; cb.Label.String = 'symbol index';
saveas(gcf, fullfile(figDir, 'fig7_train_fit.png'));

% Stage 7: PLL on PRE
figure(8); clf; set(gcf,'Name','PRE: Before vs After PLL');
subplot(2,2,1); plot(real(results.v2.preTake3), imag(results.v2.preTake3), 'o');
axis equal; grid on; xline(0,'k:'); yline(0,'k:'); title('CAR 2 PRE After CFO (Before PLL)');
subplot(2,2,2); plot(real(results.v1.preTake3), imag(results.v1.preTake3), 'o');
axis equal; grid on; xline(0,'k:'); yline(0,'k:'); title('CAR 1 PRE After CFO (Before PLL)');
subplot(2,2,3); plot(real(results.v2.zPLL), imag(results.v2.zPLL), 'o');
axis equal; grid on; xline(0,'k:'); yline(0,'k:'); title('CAR 2 PRE After PLL');
subplot(2,2,4); plot(real(results.v1.zPLL), imag(results.v1.zPLL), 'o');
axis equal; grid on; xline(0,'k:'); yline(0,'k:'); title('CAR 1 PRE After PLL');
saveas(gcf, fullfile(figDir, 'fig8_pre_pll.png'));

% Stage 8: Channel estimate
figure(9); clf; set(gcf,'Name','Channel |H|');
plot(1:NumDataCarriers, abs(results.v2.H_preamble), '-o', 'DisplayName','CAR 2');
hold on;
plot(1:NumDataCarriers, abs(results.v1.H_preamble), '-o', 'DisplayName','CAR 1');
grid on; legend('Location','best');
xlabel('Data subcarrier'); ylabel('|H|');
title('Channel Estimate Per Car');
saveas(gcf, fullfile(figDir, 'fig9_channel_estimate.png'));

% Stage 9: OFDM resource grid (pre-EQ and post-MMSE)
figure(10); clf; set(gcf,'Name','OFDM Resource Grid');
subplot(2,2,1);
imagesc(1:NumOFDMSymbols, 1:NumDataCarriers, abs(results.v2.rxSubcarriers));
axis xy; colorbar; xlabel('OFDM symbol'); ylabel('Data subcarrier');
title('CAR 2 |RX| pre-equalization');
subplot(2,2,2);
imagesc(1:NumOFDMSymbols, 1:NumDataCarriers, abs(results.v1.rxSubcarriers));
axis xy; colorbar; xlabel('OFDM symbol'); ylabel('Data subcarrier');
title('CAR 1 |RX| pre-equalization');
subplot(2,2,3);
imagesc(1:NumOFDMSymbols, 1:NumDataCarriers, real(results.v2.rxEq));
axis xy; colorbar; clim([-1.5 1.5]);
xlabel('OFDM symbol'); ylabel('Data subcarrier');
title('CAR 2 Re(RX) post-MMSE  (BPSK +/-1)');
subplot(2,2,4);
imagesc(1:NumOFDMSymbols, 1:NumDataCarriers, real(results.v1.rxEq));
axis xy; colorbar; clim([-1.5 1.5]);
xlabel('OFDM symbol'); ylabel('Data subcarrier');
title('CAR 1 Re(RX) post-MMSE  (BPSK +/-1)');
saveas(gcf, fullfile(figDir, 'fig10_ofdm_resource_grid.png'));

% Stage 10: Final BPSK constellation
figure(11); clf; set(gcf,'Name','V2V RX Constellation');
subplot(1,2,1);
scatter(real(results.v2.rxEq(:)), imag(results.v2.rxEq(:)), 25, 'filled'); hold on;
plot([-1 1],[0 0], 'r+', 'MarkerSize',14, 'LineWidth',2);
xlim([-2 2]); ylim([-2 2]); axis square; grid on;
title(sprintf('Car 2  BER=%.4f  corrSNR=%.1f', results.v2.BER, results.v2.corrSNR));
xlabel('I'); ylabel('Q');
subplot(1,2,2);
scatter(real(results.v1.rxEq(:)), imag(results.v1.rxEq(:)), 25, 'filled'); hold on;
plot([-1 1],[0 0], 'r+', 'MarkerSize',14, 'LineWidth',2);
xlim([-2 2]); ylim([-2 2]); axis square; grid on;
title(sprintf('Car 1  BER=%.4f  corrSNR=%.1f', results.v1.BER, results.v1.corrSNR));
xlabel('I'); ylabel('Q');
saveas(gcf, fullfile(figDir, 'fig11_bpsk_constellation.png'));

% Stage 11: Collision avoidance (TTC, gap, brake warning)
figure(12); clf; set(gcf,'Name','V2V Collision Avoidance');
fill([85 135 135 85],[-3.5 -3.5 3.5 3.5], [0.9 0.9 0.9], 'EdgeColor','none'); hold on;
plot([85 135], [0 0], 'w--', 'LineWidth',1.5);
fill(bsm1_rx.posX + [-3 3 3 -3], bsm1_rx.posY + [-1 -1 1 1], [1 0.6 0.2], 'EdgeColor','k');
text(bsm1_rx.posX, bsm1_rx.posY, 'Car 1', 'HorizontalAlignment','center', 'FontWeight','bold');
text(bsm1_rx.posX, bsm1_rx.posY-1.8, sprintf('%.0f m/s', bsm1_rx.velocity), ...
    'HorizontalAlignment','center');
fill(bsm2_rx.posX + [-3 3 3 -3], bsm2_rx.posY + [-1 -1 1 1], [0.3 0.6 1.0], 'EdgeColor','k');
text(bsm2_rx.posX, bsm2_rx.posY, 'Car 2', 'HorizontalAlignment','center', 'FontWeight','bold');
text(bsm2_rx.posX, bsm2_rx.posY-1.8, sprintf('%.0f m/s', bsm2_rx.velocity), ...
    'HorizontalAlignment','center');
midX = (bsm1_rx.posX + bsm2_rx.posX) / 2;
plot([bsm1_rx.posX+3, bsm2_rx.posX-3], [2.2 2.2], 'k-', 'LineWidth',1.5);
plot(bsm1_rx.posX+3, 2.2, 'k<', 'MarkerFaceColor','k', 'MarkerSize',6);
plot(bsm2_rx.posX-3, 2.2, 'k>', 'MarkerFaceColor','k', 'MarkerSize',6);
text(midX, 2.7, sprintf('Gap: %.1f m', gap), 'HorizontalAlignment','center');
if brakeWarn
    text(midX, -2.5, 'WARNING: Brake Recommended', ...
        'HorizontalAlignment','center', 'FontWeight','bold', ...
        'Color',[0.85 0.4 0], 'FontSize',12);
end
xlim([85 135]); ylim([-12 12]);
xlabel('Position (m)'); ylabel('Lane (m)');
if isfinite(TTC)
    title(sprintf('V2V Collision Avoidance  |  TTC = %.1f s', TTC));
else
    title('V2V Collision Avoidance  |  Not closing');
end
grid on;
saveas(gcf, fullfile(figDir, 'fig12_collision_avoidance.png'));

fprintf('\nAll figures saved to: %s\n', figDir);

%% ----------------------- LOCAL FUNCTIONS -----------------------

function bits = encodeBSMtoBits(bsm)
    % SAE J2735, scaled-int BSM, 10 bytes total
    velScaled  = uint16(round(bsm.velocity * 100));
    accScaled  = int8(round(bsm.accel * 10));
    posXScaled = uint16(round(bsm.posX * 10));
    posYScaled = int16(round(bsm.posY * 10));

    payload = uint8(zeros(1,10));
    payload(1)  = uint8(bsm.vehicleID);
    payload(2)  = uint8(bitshift(velScaled, -8));
    payload(3)  = uint8(bitand(velScaled, uint16(255)));
    payload(4)  = typecast(accScaled, 'uint8');
    payload(5)  = uint8(bitshift(posXScaled, -8));
    payload(6)  = uint8(bitand(posXScaled, uint16(255)));
    posYbytes   = typecast(posYScaled, 'uint8');
    payload(7)  = posYbytes(2);
    payload(8)  = posYbytes(1);
    payload(9)  = uint8(0);
    payload(10) = uint8(0);

    bits = double(reshape(de2bi(payload, 8, 'left-msb')', [], 1));
end

function bsm = decodeBSMfromBits(rxBits)
    % inverse of encodeBSMtoBits
    rxBytes = uint8(bi2de(reshape(rxBits, 8, 10)', 'left-msb')).';
    bsm.vehicleID = double(rxBytes(1));

    velScaled = uint16(rxBytes(2)) * 256 + uint16(rxBytes(3));
    bsm.velocity = double(velScaled) / 100;

    accScaled = typecast(uint8(rxBytes(4)), 'int8');
    bsm.accel = double(accScaled) / 10;

    posXScaled = uint16(rxBytes(5)) * 256 + uint16(rxBytes(6));
    bsm.posX = double(posXScaled) / 10;

    posYbytes = uint8([rxBytes(8), rxBytes(7)]);
    posYScaled = typecast(posYbytes, 'int16');
    bsm.posY = double(posYScaled) / 10;
end