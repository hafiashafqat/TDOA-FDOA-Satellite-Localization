# Automated Satellite-Based RF Emitter Localization

## Overview

This project develops a MATLAB-based simulation framework for satellite-based Radio Frequency (RF) emitter localization using Time Difference of Arrival (TDOA) and Frequency Difference of Arrival (FDOA) techniques.

The framework integrates Systems Tool Kit (STK) with MATLAB to simulate satellite orbital motion, determine common satellite visibility, extract satellite states, generate synthetic TDOA/FDOA measurements, and estimate the geographical location of a ground-based RF emitter.

The project consists of 11 modules, covering the complete workflow from scenario initialization and satellite configuration to localization, accuracy assessment, and Monte Carlo robustness testing.

## Objectives

- Develop a satellite-based RF emitter localization simulation framework.
- Integrate STK 11 with MATLAB R2023b for satellite scenario modelling and data extraction.
- Simulate multiple satellites observing a common ground location.
- Generate synthetic TDOA and FDOA measurements from satellite-emitter geometry.
- Implement a cumulative localization algorithm to estimate emitter latitude and longitude.
- Evaluate localization accuracy, residuals, convergence, and uncertainty.
- Assess system robustness under different measurement noise levels, satellite configurations, observation durations, and model uncertainties.

## Methodology and Project Workflow

The project is organized into the following 11 modules:

| Module | Description |
|---|---|
| 01. STK–MATLAB Initialization | Establishes the connection between MATLAB and STK and initializes the simulation scenario. |
| 02. Satellite Creation and Propagation | Creates three satellites and configures their orbital parameters and propagation. |
| 03. Satellite Property Inspection | Retrieves and examines satellite properties and orbital information. |
| 04. Satellite Geometry and Visibility | Determines satellite visibility, access intervals, and common visibility periods over the ground emitter location. |
| 05. Satellite State Extraction | Extracts satellite position and velocity information at common observation epochs. |
| 06. Emitter Geometry and RIC Transformation | Defines the ground emitter geometry, performs coordinate transformations, and calculates relative satellite-emitter geometry. |
| 07. Synthetic TDOA/FDOA Generation | Generates synthetic TDOA and FDOA measurements with assumed Gaussian measurement noise. |
| 08. Measurement Validation | Performs consistency checks on satellite states, geometric calculations, and generated measurements. |
| 09. Cumulative Localization Solver | Estimates the emitter's latitude and longitude using cumulative TDOA and FDOA measurements. |
| 10. Accuracy and Residual Analysis | Evaluates localization error, residual statistics, uncertainty, convergence, and sensitivity to initial guesses. |
| 11. Monte Carlo Robustness Testing | Repeats localization simulations under different conditions to evaluate reliability and sensitivity to measurement and model uncertainties. |

## Simulation Configuration

The current simulation uses the following configuration:

| Parameter | Value |
|---|---|
| Software | MATLAB R2023b |
| Satellite Simulation | STK 11 |
| Number of Satellites | 3 |
| Satellite Altitude | 650 km |
| Orbital Inclination | 60° |
| Carrier Frequency | 2.4 GHz |
| TDOA Noise Standard Deviation | 1 microsecond |
| FDOA Noise Standard Deviation | 1 Hz |
| Ground Emitter Latitude | 31.00° |
| Ground Emitter Longitude | 73.50° |
| Ground Emitter Altitude | 0.2 km |
| Observation Epochs | 9 |

## Localization Approach

The localization process uses TDOA and FDOA measurements derived from the relative positions and velocities of the satellites and the simulated ground emitter.

The solver evaluates candidate emitter locations by predicting their corresponding TDOA and FDOA measurements and comparing them with the generated measurements.

An optimization process then refines the estimated latitude and longitude to minimize the measurement residuals.

Measurements from successive observation epochs are accumulated to improve the localization estimate.

The current implementation estimates two geographical parameters:

- Latitude
- Longitude

The emitter altitude is assumed to be known and fixed at 0.2 km.

## Simulation Results

The current simulation produced the following baseline results:

- **Single-run localization error:** Approximately 8.3 metres.
- **Monte Carlo trials:** 100 per test case.
- **Baseline Monte Carlo RMS error:** Approximately 13.9 metres.
- **Baseline failure rate:** 0%.
- **95% uncertainty coverage:** 95%.
- **Formal covariance calibration:** Consistent with the empirical Monte Carlo results.

Additional tests examined the effects of measurement noise, reference satellite selection, observation duration, satellite count, emitter position, assumed altitude error, and carrier-frequency mismatch.
##Results and Analysis
Module 10 – Accuracy and Residual Analysis
The localization solver achieved a horizontal positioning error of approximately 8.3 metres in the simulated scenario using nine measurement epochs.

The residual analysis produced a chi-square per degree of freedom of 0.99, indicating consistency with the assumed measurement noise model. The solver also demonstrated improved convergence when additional measurement epochs were included, with the success rate increasing from 52% for one epoch to 93% for nine epochs across 69 initial guesses.

These results demonstrate the benefit of accumulating measurements over time, although some sensitivity to initial guesses remains.
Module 11 – Monte Carlo and Robustness Testing
Monte Carlo simulations were conducted using 100 trials per test configuration to evaluate the system under different measurement noise levels, satellite geometries, reference satellite selections, and emitter locations.

The analysis examined four main aspects:
- Measurement noise: Localization error generally increased as TDOA and FDOA noise levels increased.
- Covariance calibration: Compared empirical localization errors with the uncertainty predicted by the solver to evaluate uncertainty estimation.
- Altitude sensitivity: Incorrect assumptions about emitter altitude affected horizontal localization accuracy.
- Carrier frequency sensitivity: Carrier frequency mismatch influenced FDOA predictions and consequently affected localization performance.

These results demonstrate the performance of the implemented localization framework under the simulated conditions and assumptions.

## Repository Structure

The repository is organized around the 11 individual modules:

```text
Automated-Satellite-Emitter-Localization/
│
├── README.md
│
├── Module01_Initialization/
├── Module02_Satellite_Creation/
├── Module03_Satellite_Inspection/
├── Module04_Satellite_Visibility/
├── Module05_State_Extraction/
├── Module06_Emitter_Geometry/
├── Module07_TDOA_FDOA_Generation/
├── Module08_Measurement_Validation/
├── Module09_Localization_Solver/
├── Module10_Accuracy_Analysis/
├── Module11_Monte_Carlo/
│
├── Results/
│
└── Figures/
```

*The folder names above are a suggested organization; adjust them to match the actual repository files.*

## Requirements

- MATLAB R2023b
- STK 11
- MATLAB–STK COM connectivity on Windows
- Optimization Toolbox for the `lsqnonlin` solver option, if used
- Statistics and Machine Learning Toolbox is not required for the manual residual-correlation calculation.

## How to Run

1. Install MATLAB R2023b and STK 11.
2. Open STK and establish the simulation environment.
3. Open the project folder in MATLAB.
4. Run the modules in numerical order, from Module 01 to Module 11.
5. Ensure that the output files generated by each module are available to the subsequent modules.
6. Review the generated CSV files, localization results, statistical summaries, and figures.

**Note:** Some modules depend on outputs generated by earlier modules. They should therefore be executed in the correct sequence.

## Future Research Directions

Beyond its current application to satellite-based RF emitter localization, this framework provides a foundation for exploring other problems involving multi-observer measurements, signal analysis, and source localization.

One potential research direction is investigating the mathematical and computational methods developed in this project in the context of **gravitational-wave astronomy**.

Gravitational-wave observatories, such as LIGO, Virgo, and KAGRA, use measurements from geographically separated detectors to study gravitational-wave signals and constrain the locations of their astrophysical sources.

Some concepts explored in this project, particularly multi-observer geometry, measurement differences, uncertainty estimation, residual analysis, and localization optimization, may provide a useful starting point for studying related computational approaches in gravitational-wave source localization.

Potential extensions include:

- **Multi-detector source localization:** Exploring optimization and geometric methods for estimating source locations using measurements from multiple gravitational-wave detectors.
- **Time-delay-based localization:** Studying the use of differences in gravitational-wave arrival times between detectors to constrain source directions.
- **Signal analysis and parameter estimation:** Investigating computational techniques for extracting source information from noisy measurements.
- **Uncertainty and robustness analysis:** Applying Monte Carlo simulations and statistical validation to assess the reliability of source-location estimates under different noise conditions.
- **Cross-domain computational methods:** Exploring how optimization, numerical modelling, and multi-observer measurement techniques can be adapted across satellite-based localization and gravitational-wave data analysis.

---

**Technologies:** MATLAB | STK | TDOA | FDOA | Satellite Orbital Mechanics | Numerical Optimization | Monte Carlo Simulation
