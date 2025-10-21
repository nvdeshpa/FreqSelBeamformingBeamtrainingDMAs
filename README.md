# Frequency-Selective Beamforming and Single-Shot Beam Training with Dynamic Metasurface Antennas

MATLAB simulation code for generating results from the paper "Frequency-selective beamforming and single-shot beam training with dynamic metasurface antennas."

**Authors:** Nitish Vikas Deshpande, Joseph Carlson, Miguel R. Castellanos, and Robert W. Heath Jr.

## Overview

This repository contains MATLAB code that simulates dynamic metasurface antenna (DMA) systems for wideband beamforming and beam training applications. The code compares the performance of frequency-selective beamforming with DMAs against traditional approaches including:
- True Time Delay (TTD) arrays
- Phase shifter arrays
- Binary PIN diode-based beamforming

## System Description

The simulation models a DMA system operating at:
- **Center frequency:** 15 GHz
- **Bandwidth:** 6 GHz (configurable)
- **Number of subcarriers:** 64 (configurable)
- **Array configurations:** Linear and planar arrays with configurable dimensions

### Key Features

1. **Frequency-Selective Beamforming:** Optimal target frequency selection for each beam angle
2. **Single-Shot Beam Training:** Estimation of target angles from frequency-domain responses
3. **Non-Ideal DMA Modeling:** Includes attenuation effects in metasurface wave propagation
4. **Rate Analysis:** Achievable rate calculations under various bandwidth and array configurations
5. **Comparative Benchmarking:** Performance comparison with TTD, phase shifters, and binary beamforming

## Code Structure

### Main Script: `DMA_main.m`

The script is organized into several sections:

#### 1. Universal Constants and Parameters
- Physical constants (speed of light, center frequency, wavelength)
- DMA-specific parameters (coupling factor, damping factor)
- Array dimensions (Ny, Nz) and configuration

#### 2. Linear Array Beamforming Analysis
- Computes optimal target frequencies for different beam angles
- Generates beamforming gain patterns
- Compares continuous vs. binary weight configurations
- Analyzes ideal and non-ideal DMA responses

#### 3. Planar Array Simulation
- Creates 2D beamforming patterns in angle-frequency domain
- Implements frequency-selective target beam allocation
- Visualizes radiation patterns as 3D surfaces

#### 4. Beam Training and Angle Estimation
- Identifies optimal frequency for each target angle
- Estimates angle from frequency response
- Analyzes estimation error across different angles

#### 5. Achievable Rate Calculations
- Computes data rates for various system configurations
- Sweeps bandwidth (50-500 MHz) and training bandwidth (0.25-6 GHz)
- Compares:
  - Proposed optimal frequency tuning
  - Beam-trained frequency tuning
  - Fixed frequency tuning (continuous and binary weights)
  - TTD benchmark
  - Phase shifter benchmark

### Helper Functions

#### `generate_alpha_weight_vector(f, f_res, F, Gamma, N, K)`
Generates DMA weight vectors based on resonant frequencies.

**Parameters:**
- `f`: Operating frequency
- `f_res`: Resonant frequencies for each antenna element
- `F`: Coupling factor
- `Gamma`: Damping factor
- `N`: Number of elements
- `K`: Number of frequency points

#### `generate_alpha_weight(f, f_0, Gamma, F)`
Computes individual DMA element weight for a given frequency and resonant frequency.

#### `compute_phi_target_arr(ng, Ny, theta_0, delta_dB)`
Calculates target beam angles for planar array configuration.

**Parameters:**
- `ng`: Refractive index
- `Ny`: Number of elements in y-dimension
- `theta_0`: Maximum steering angle
- `delta_dB`: Beam separation in dB

#### `compute_achievable_rate_LOS_channel(...)`
Calculates achievable data rate for line-of-sight channel.

**Includes:**
- Multi-carrier OFDM modeling
- SNR calculations per subcarrier
- Support for DMA, TTD, and phase shifter architectures
- Optional PIN diode binary weight configuration

#### `generate_intrinsic_phase_vector(f, ng, dy, N, c, atten_coeff, include_attenuation)`
Computes intrinsic phase shift due to wave propagation within the metasurface.

#### `generate_extrinsic_phase_vector(f, dy, N, c, phi)`
Computes extrinsic phase shift due to spatial propagation in free space.

## Usage

### Running the Simulation

1. Open MATLAB
2. Navigate to the repository directory
3. Run the main script:
   ```matlab
   DMA_main
   ```

### Configuration Options

#### Array Configuration
Modify `example_no` to switch between different array configurations:
- `example_no = 0`: 8×4 array (Ny=8, Nz=4)
- `example_no = 1`: 4×8 array (Ny=4, Nz=8)

#### Attenuation Modeling
Toggle non-ideal effects:
```matlab
include_attenuation = 1;  % 1: include attenuation, 0: ideal case
atten_coeff = 6;           % Attenuation coefficient
```

#### Bandwidth and Subcarriers
```matlab
B = 6*10^9;        % System bandwidth (Hz)
K = 64;            % Number of subcarriers
B_arr = [50:25:500]*10^6;  % Data bandwidth sweep (Hz)
Tr_arr = [0.25:0.25:6]*10^9;  % Training bandwidth sweep (Hz)
```

## Output Visualizations

The code generates multiple figures:

1. **Optimal target frequency vs. beam angle**
2. **Beamforming gain comparison** (optimal frequency, center frequency, binary weights)
3. **3D angle-frequency radiation pattern** (surface plot)
4. **Maximum beamforming gain vs. angle**
5. **Optimal frequency index vs. angle**
6. **Estimated vs. actual angle**
7. **Angle estimation error**
8. **Normalized beamforming gain** (estimated vs. exact angles)
9. **Achievable rate vs. training bandwidth** (for different data bandwidths)

## System Requirements

- MATLAB (tested on R2019b and later)
- No additional toolboxes required

## Key Results

The simulation demonstrates:
- **Frequency-selective beamforming** significantly improves wideband performance over fixed-frequency designs
- **Single-shot beam training** enables rapid angle estimation using frequency diversity
- **DMAs outperform phase shifters** for wideband applications
- **DMAs achieve comparable performance to TTD arrays** with simpler hardware

## License

See `LICENSE` file for details.

## Citation

If you use this code in your research, please cite:

```
N. V. Deshpande, J. Carlson, M. R. Castellanos, and R. W. Heath Jr., 
"Frequency-selective beamforming and single-shot beam training with dynamic metasurface antennas," 
[Journal/Conference details to be added]
```

## Contact

For questions or issues, please open an issue on the GitHub repository.
