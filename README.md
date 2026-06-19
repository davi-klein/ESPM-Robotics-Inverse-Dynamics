# Advanced Inverse Dynamics for Serial and Parallel Manipulators

![MATLAB](https://img.shields.io/badge/MATLAB-R2023a%2B-blue.svg)
![Robotics](https://img.shields.io/badge/Domain-Robotics%20%26%20Control-orange.svg)
![Architecture](https://img.shields.io/badge/Architecture-Modular%20%26%20Vectorized-success.svg)

This repository contains a production-ready MATLAB physics engine for computing the **Inverse Dynamics** of increasingly complex robotic topologies: a **2R planar serial robot**, a **3-RRR planar parallel robot**, and a spatial **6-SPS Stewart-Gough platform**. 

The physics engine is strictly built upon **D'Alembert's Principle** combined with **Equimomental Systems of Point Masses** and scaled using **Davies' Method (Graph Theory)** for robust closed-chain constraint solving.

## Mathematical and Engineering Highlights

Unlike standard academic scripts, this repository has been refactored to meet strict software engineering standards, focusing on modularity, memory optimization, and separation of concerns:

* **Topology-Agnostic Plotting:** The visualization layer is decoupled from the physics engine.
* **Vectorized 3D Coriolis Kinematics:** Complex spatial tensors and pseudo-inertia matrices are handled efficiently without massive for-loops.
* **Singularity-Free Kinematics:** The 6-SPS platform uses a custom 18-DOF Newton-Raphson solver based on **Natural Coordinates** to completely avoid trigonometric singularities commonly found in Euler-angle approaches.
* **Graph Theory Cutset Matrices:** Davies' method handles up to 13 simultaneous rigid bodies and 98 exact cutsets (in the 6-SPS model) with strict algebraic parity.

## Repository Architecture

The project is structured to separate the orchestration flow from the heavy mathematical modules:

```text
robotics-inverse-dynamics/
├── data/                       # Pre-calculated trajectory tensors and ground-truth validation (GIM/Tsai)
├── src/                        # The Physics API (Modular toolset)
│   ├── init_robot_*_params.m   # Data models (Geometry, Mass, Pseudo-Inertia)
│   ├── compute_kinematics_*.m  # Jacobians and Newton-Raphson solvers
│   ├── compute_point_masses_*.m# Non-rigid spatial transformations
│   ├── compute_inverse_*.m     # Graph-theory based constraint solver (Davies)
│   └── plot_simulation.m       # Decoupled visualization module
├── main_2R_simulation.m        # Orchestrator for the 2-Link Serial Robot
├── main_3RRR_simulation.m      # Orchestrator for the 3-RRR Parallel Robot
└── main_6SPS_simulation.m      # Orchestrator for the 6-SPS Spatial Platform
```

## Getting Started

No external heavy toolboxes (like Simscape) are required. The engine relies purely on base MATLAB linear algebra operations for maximum performance.

1. Clone the repository:
   ```bash
   git clone [https://github.com/davi-klein/ESPM-Robotics-Inverse-Dynamics.git](https://github.com/davi-klein/ESPM-Robotics-Inverse-Dynamics.git)
   cd ESPM-Robotics-Inverse-Dynamics
   ```
2. Open MATLAB and ensure you are in the root directory of the project.
3. Run any of the main orchestrator scripts directly. For example:
   ```matlab
   run('main_6SPS_simulation.m')
   ```
4. Standardized plots detailing the required actuation forces/torques over time will be generated automatically.

## Author

**Davi Klein**
M.Sc. Student in Computer Science | Robotics Researcher
Federal University of Santa Maria (UFSM)
