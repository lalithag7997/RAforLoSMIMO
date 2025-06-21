
# Reflect-Arrays for Long-Range LoS MIMO Links

This repository contains a MATLAB-based simulation framework for evaluating Reflect-Array (RA) aided Line-of-Sight (LoS) MIMO systems designed for long-range, high-capacity wireless links. The approach leverages large virtual apertures formed by placing reflect-arrays near transceivers, enabling spatial degrees of freedom (DoF) even when physical form factors are constrained.

Reflect-arrays placed close to the transmitter and receiver require:

* Quadratic phase compensation for beam focusing, and
* Linear phase compensation for beamforming toward the corresponding RA/receiver.

This dual-phase compensation strategy allows the system to maintain compact transceiver size while delivering both beamforming gain and spatial multiplexing performance.

This technique enabled power-efficient transfer of 5.2 Gbps over a 1,500 m LoS link in the 28 GHz band.

Please refer to the full paper [here](https://wcsl.ece.ucsb.edu/sites/default/files/publications/creating_spatial_degrees_of_freedom_for_long_range_los_mimo_using_reflect_arrays_spawc_2024_final.pdf) for detailed system modeling and experimental insights.

---

## Simulation Overview

This simulation models:

* A 4×4 LoS MIMO system with 4 independent QPSK data streams
* RA-based virtual aperture creation with user-defined reflect-array size (`M`) and focal depth (`δ_z`)
* End-to-end channel construction that includes:

  * Direct LoS paths
  * RA-assisted paths with phase compensation (linear + quadratic)
  * Receiver-side displacement modeled as random ±2 m offsets

---

## How to Run

### Requirements

* MATLAB (R2021a or later recommended)
* Communications Toolbox

### Files

* `RA_LoS_MIMO_simulation.m`: Main simulation script

### Usage

1. Open the `.m` file in MATLAB.
2. Modify the following parameters as needed:

   * `M` → Reflect-array size
   * `delta_z` → Depth parameter
   * `Nrun` → Number of Monte Carlo runs
   * Channel geometry settings and offsets
3. Run the script.

---

## Inputs

| Parameter       | Description                                                   |
| --------------- | ------------------------------------------------------------- |
| `M`             | Number of elements in each reflect-array                      |
| `delta_z`       | Vertical depth placement of RAs (focus control)               |
| `Nrun`          | Number of Monte Carlo trials for averaging BER                |
| QPSK Streams    | 4 parallel symbol streams, mapped across 4×4 LoS MIMO link    |
| RX Displacement | Random ±2 m horizontal/vertical offsets applied to RX and RAs |

---

## Outputs

* BER vs SNR curve for the RA-aided 4×4 LoS MIMO system
* Comparison with a benchmark SISO system
* Insight into the performance gains from spatial DoF and beamforming via RA-assisted virtual aperture synthesis


