# V2V BSM OFDM - ADALM-PLUTO SDR

**Authors:** Avery Matherne, Jason Phan

A vehicle-to-vehicle (V2V) basic safety message (BSM) transmission system implemented in MATLAB. The system encodes, modulates, transmits, and decodes BSMs from two simulated vehicles over a complete 802.11p-style OFDM physical layer with a full synchronization pipeline. It runs on ADALM-PLUTO SDR hardware (TX on `usb:0`, RX on `usb:1`) or in software simulation with a Rician fading channel.

---

## Overview

The transmitter constructs two BSM packets - one per vehicle - encodes each with a rate-1/2 convolutional code, maps them onto BPSK-OFDM, prepends an RRC-shaped preamble and a known training sequence for synchronization, and transmits both consecutively. The receiver detects each packet via cross-correlation, sweeps for sample-clock offset, estimates and corrects carrier frequency offset (CFO), fits and removes residual linear phase using the training sequence, applies a decision-directed PLL for fine phase lock, equalizes the OFDM subcarriers with MMSE, and decodes via Viterbi to recover the original BSM fields. Decoded kinematics from both vehicles are used to compute time-to-collision (TTC) and issue a brake warning.

---

## System Architecture

```
[BSM struct] --> encodeBSMtoBits --> ConvEnc(K=7,r=1/2) --> BPSK-OFDM modulator
                                                                         |
                        [RRC PRE (63 sym)] + [TRAIN (256 sym)] + [OFDM Preamble] + [OFDM Data]
                                                                         |
                                              CFO rotation applied, waveform normalized
                                                                         |
                                              ADALM-PLUTO TX (usb:0, 2.9 GHz) or AWGN sim
                                                                         |
                                              ADALM-PLUTO RX (usb:1, 2.9 GHz) or sim channel
                                                                         |
                        PRE cross-corr --> tau-sweep --> CFO est --> TRAIN fit --> PLL
                                                                         |
                                              OFDM preamble channel est --> MMSE EQ
                                                                         |
                                              Viterbi decode --> decodeBSMfromBits
                                                                         |
                                              TTC / gap / brake warning output
```

Each transmitted waveform is structured as:

```
| RRC PRE (63 sym × 16 sps) | TRAIN (256 sym × 16 sps) | OFDM Preamble (1 sym) | OFDM Data (N sym) | Guard |
```

Car 2's packet is transmitted first, immediately followed by Car 1's packet. The receiver tries each detected peak against both vehicle reference bit sequences and assigns by lower BER.

---

## Signal Processing Pipeline

### BSM Encoding (`encodeBSMtoBits`)

Messages follow a SAE J2735-style scaled-integer format packed into 10 bytes (80 bits):

| Byte(s) | Field      | Encoding                        |
|---------|------------|---------------------------------|
| 1       | Vehicle ID | raw `uint8`                     |
| 2–3     | Velocity   | `uint16`, scaled ×100 (m/s)     |
| 4       | Accel      | `int8`, scaled ×10 (m/s²)       |
| 5–6     | posX       | `uint16`, scaled ×10 (m)        |
| 7–8     | posY       | `int16`, scaled ×10 (m), big-endian |
| 9–10    | Reserved   | zeros                           |

Bits are packed MSB-first. The decoder (`decodeBSMfromBits`) is the exact inverse.

### Forward Error Correction

A K=7, rate-1/2 convolutional code with generator polynomials `[171, 133]` (octal) - the standard 802.11 FEC - is applied to the 80-bit BSM payload. This doubles the bit count to 160 coded bits. Decoding uses Viterbi with hard decisions and a traceback depth of 35.

### OFDM Modulation (802.11p PHY)

| Parameter              | Value                          |
|------------------------|--------------------------------|
| FFT length             | 64 points                      |
| Bandwidth              | 10 MHz                         |
| Data subcarriers       | 48                             |
| Guard band carriers    | 6 (lower) / 5 (upper)         |
| DC null                | yes                            |
| Pilot subcarriers      | 4 (indices 12, 26, 40, 54)    |
| Cyclic prefix length   | 16 samples (1.6 µs @ 10 MHz) |
| Modulation             | BPSK (one bit per subcarrier) |

A separate OFDM preamble symbol (seed 43, known BPSK data) precedes the data OFDM symbols for channel estimation.

### RRC Preamble (PRE) - Timing Detection

A 63-symbol BPSK sequence (seeded with `rng(42)`) is RRC pulse-shaped with rolloff=0.35, sps=16 (samples per symbol), and span=6. The receiver cross-correlates the raw received signal against the known shaped preamble to detect packet boundaries. Peak positions above a threshold of 2× the median correlation floor are accepted. The receiver also predicts additional candidate peaks at ±N_pkt and ±2·N_pkt offsets from each found peak to recover packets missed by `findpeaks`.

A tau-sweep over all 0–15 sample offsets resolves sub-sample timing ambiguity by maximizing correlation energy against the known preamble symbols at the matched-filter output.

### TRAIN Sequence - Linear Phase and Polarity Correction

A 256-symbol known BPSK sequence (same RRC shaping) follows the PRE. After CFO removal, the receiver extracts the training symbols, normalizes their RMS, and computes the channel response `H = rx .* conj(ref)`. The unwrapped phase of H is fit to a line via `polyfit` to estimate and remove residual linear phase (a second-order CFO residual). The sign of the projection onto the reference determines and corrects a ±1 polarity ambiguity. The fitted rotation is then extended over the remainder of the receive buffer.

### CFO Estimation and PLL

**Coarse CFO:** Differential phase is computed across consecutive PRE symbols after matched filtering and symbol-rate downsampling. The angle of the summed product gives a per-sps phase increment; CFO = φ/(2π) × (Fs/sps). In hardware mode, estimates above 200 kHz are rejected as implausible.

**Fine PLL:** A decision-directed BPSK Costas loop (loop gain β=0.001, 8-symbol burn-in) is run on the PRE symbols after coarse CFO and TRAIN correction. The residual phase offset is measured by projecting the PLL output onto the known preamble and subtracted globally.

### MMSE Channel Equalization

The known OFDM preamble symbol is demodulated and divided by the reference to obtain a per-subcarrier channel estimate H. This estimate is replicated across all data OFDM symbols and applied as:

```
rxEq = (conj(H) .* rx) ./ (|H|^2 + N0)
```

where N0 is estimated from a 200-sample noise window prior to the detected PRE. This is minimum mean-square error (MMSE) equalization; zero-forcing is the special case N0 → 0.

### Collision Avoidance Output

From the decoded BSMs, the receiver computes:

- **Gap:** posX(Car 2) − posX(Car 1)
- **Closing rate:** velocity(Car 1) − velocity(Car 2)
- **TTC:** gap / closing rate (infinite if not closing)
- **Brake warning:** issued when TTC < 5.0 s

---

## Hardware Requirements

- Two ADALM-PLUTO SDRs (one TX, one RX) connected via USB
- MATLAB Communications Toolbox
- MATLAB Support Package for ADALM-PLUTO Radio
- USB cables; short coax or antenna pair for the RF link

The system uses 2.9 GHz as the center frequency (instead of the DSRC band at 5.9 GHz) because the ADALM-PLUTO's AD9363 transceiver does not reach 5.9 GHz. TX gain is set to −10 dB and RX gain to 40 dB (manual).

---

## Running the Code

### Simulation Mode (default)

```matlab
% In V2V_TXRXf.m, line 24:
USE_HARDWARE = false;
SNR_dB       = 30;     % adjust to test link margin
```

Run `V2V_TXRXf` from the MATLAB command window. The simulation passes the waveform through a Rician fading channel (K=3, 500 Hz max Doppler, three paths at delays 0 / 100 ns / 250 ns with gains 0 / −5 / −9 dB) followed by AWGN.

### Hardware Mode

```matlab
USE_HARDWARE = true;
```

Connect TX Pluto to `usb:0` and RX Pluto to `usb:1`. Adjust `TX_GAIN` and `RX_GAIN` for the physical link budget. The transmitter calls `transmitRepeat` so the waveform loops continuously while the receiver captures `FRAMES_TO_CAPTURE` frames (default 40 × 8192 samples ≈ 32.8 ms).

---

## Output and Figures

Console output includes per-attempt BER against each vehicle reference, decoded BSM field comparison (truth vs. received), gap, closing rate, TTC, and brake warning status. All 12 diagnostic figures are saved as PNG to `./V2V_figures/`:

| Figure | Content |
|--------|---------|
| fig1   | TX OFDM baseband power spectrum |
| fig2   | TX OFDM spectrum shifted to RF center frequency |
| fig3   | RX captured spectrum (hardware mode only) |
| fig4   | PRE cross-correlation with detected packet peaks annotated |
| fig5   | tau-sweep score vs. offset for both vehicles |
| fig6   | PRE constellation before and after CFO correction |
| fig7   | TRAIN constellation before and after linear-phase + polarity fit |
| fig8   | PRE constellation before and after PLL |
| fig9   | Per-subcarrier channel magnitude estimate |H| for both vehicles |
| fig10  | OFDM resource grid (magnitude pre-EQ, real part post-MMSE) |
| fig11  | Final BPSK received constellation with BER annotation |
| fig12  | Top-down collision avoidance visualization with TTC and gap |
