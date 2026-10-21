# SAR Autofocus Toolbox (MATLAB)

High-resolution Stripmap SAR processing chain: Range-Doppler focusing, Depth of Focus analysis, and robust autofocus.

**Features:**
- Exact range-dependent azimuth compression
- Depth of Focus analysis (5 vs 30 groups) with phase error < pi/4
- PGA with LUMV estimator and adaptive windowing
- Multi-Look Registration (MLR) via cross-correlation
- Tested on both clean and noisy (SNR) raw data
- RCMC via range interpolation
- Contrast-based autofocus (coarse-to-fine search on beta)

**Results:** 
<img width="700" height="525" alt="Figura_2" src="https://github.com/user-attachments/assets/fa29c561-3503-4a42-bc01-6ef048dde853" />
<img width="700" height="525" alt="Figura_3_1" src="https://github.com/user-attachments/assets/b7494b42-e62e-44d2-b3ad-8b2c127912f0" />

## Part 2: SAR Image Analysis
- Oversampling 4x via zero-padding
- Sidelobe control with 2D Hamming
- Multilook 3x3, ENL estimation
- Gamma PDF vs empirical histogram (natural vs anthropogenic)
- Frost despeckling with optimal window selection
