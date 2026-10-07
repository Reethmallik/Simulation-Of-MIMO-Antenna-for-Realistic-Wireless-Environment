# MIMO Channel Simulation: Rician Fading, Spatial Correlation and SVD Beamforming

Simulation of a MIMO wireless link under realistic channel conditions, developed as part of my MEng thesis at Dublin City University: *Simulation of MIMO Antenna Systems for Realistic Wireless Environments*.

## What it models

- 2x2 MIMO link at 2.4 GHz, 10 MHz bandwidth
- Rician fading (K = 10) with Tx/Rx spatial correlation (0.7 / 0.5)
- Path loss (3GPP urban macro LOS style) with 6 dB log-normal shadowing
- Imperfect channel state information (CSI error variance 0.05)
- SVD-based beamforming using the estimated channel
- Thermal noise from -174 dBm/Hz with a 7 dB noise figure
- Antenna radiation pattern (boresight gain used in the link budget)

## Outputs

- Received power vs distance
- Capacity (bps/Hz) vs distance
- Outage probability vs distance (target rate 10 bps/Hz)
- Antenna radiation pattern (polar plot)

Each distance point is averaged over 1000 random channel trials, across 50 distances from 10 m to 500 m.

## Results

Received power falls log-linearly with distance. Capacity degrades as the effective SNR drops, and outage probability rises beyond a critical range. The antenna pattern shows maximum gain at boresight.

## Files

| File | Description |
|------|-------------|
| `mimo_simulation.m` | Original MATLAB simulation |
| `mimo_simulation.html` | Browser version of the same simulation (JavaScript). Open the file in any browser, no installation needed |

## How to run

**MATLAB:** open `mimo_simulation.m` and press Run. The Parallel Computing Toolbox is optional.

**Browser:** double-click `mimo_simulation.html`. The simulation runs automatically and draws the four plots.

## Notes

The HTML file is a port of the MATLAB code, so exact numbers differ slightly because of the different random number generator. The trends are the same.
