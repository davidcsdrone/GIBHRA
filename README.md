# GIBHRA: Geometric Inversion-Based Hole Restoration Algorithm

[![Conference](https://img.shields.io/badge/Published-IEEE%20WMNC%202025-blue.svg)](https://doi.org/10.23919/WMNC67099.2025.11299275)
https://ieeexplore.ieee.org/document/11299275/ 
## Overview
This repository contains the MATLAB simulation code and algorithmic implementation for the **Geometric Inversion-Based Hole Restoration Algorithm (GIBHRA)**. 

In Wireless Sensor Networks (WSNs), random node deployment often results in "coverage holes" or unmonitored areas, which severely degrade network performance.
Existing healing algorithms frequently suffer from high computational overhead or yield sub-optimal placements for restorative nodes. 

GIBHRA addresses these limitations by formulating the hole restoration problem as the geometric challenge of finding a maximally inscribed circle perfectly tangent 
to the three boundary nodes defining a coverage void. By leveraging the principles of geometric inversion, GIBHRA calculates the mathematically optimal location for 
a new restorative node, maximizing coverage while minimizing deployment cost.

## Key Contributions & Results
Simulations conducted in MATLAB on a 1000m x 1000m network area demonstrate that GIBHRA is decisively superior to benchmark methods (HCHA, HPA, NLCHR).
* **Coverage Efficiency:** Achieves **99.2%** $k$-coverage efficiency in sparse networks (compared to 50.8% for HPA).
* **Resource Efficiency:** Yields a **48% reduction** in physical sensor deployment costs to achieve 98% network coverage.
* **Restoration Capacity:** Consistently restores an average of approximately **2800 m²** of newly covered area per deployed node.
* **Computational Optimization:** Reduces the combinatorial search process time complexity from $O(n^3)$ to $O(m \log m)$ via a Delaunay triangulation-based approach.

## Algorithm Architecture
The GIBHRA framework operates in three distinct phases:

### Phase 1: Boundary Node Identification
The algorithm filters the entire set of network nodes to identify the specific subset of nodes that actively border a coverage hole. This pre-processing step significantly prunes the search space, preventing wasted computational resources on interior nodes.

### Phase 2: Adaptive Triplet Selection via Delaunay Triangulation
GIBHRA utilizes a Delaunay triangulation approach to partition boundary points into a mesh of non-overlapping triangles. It employs an adaptive distance threshold ($T_d$) to ensure all identified triplets consist of mutually non-overlapping and non-intersecting nodes, providing a well-defined geometric basis for restoration.

### Phase 3: The Geometric Inversion Framework
To solve for a restorative circle tangent to three non-tangent boundary nodes, the algorithm applies a spatial transformation:
1. **Temporary Tangency:** A uniform radial expansion ($r_{new}$) forces the two closest nodes into a tangent configuration.
2. **Inversion:** A geometric inversion is performed centered at this new point of tangency, mapping the tangent circles into two parallel lines and the third circle into a new distinct circle.
3. **Inverted Solution:** The center and radius of the restorative circle are easily calculated in this simplified linear space.
4. **Re-inversion:** The inverse transformation is applied to map the solution back to the original Euclidean plane, yielding the precise optimal coordinates $(x_o, y_o)$ for the new sensor deployment.

## Repository Structure
* `/algorithms`: Contains the core logic for GIBHRA and the benchmark comparative models (e.g., DDPCH, HORA).
* `/simulations`: Execution scripts for running the coverage, efficiency, and cost comparison pipelines.
* `/src`: Standalone helper functions, mathematical operations, and spatial plotting scripts.
* `/results`: Output figures, visualizations, and comparative data charts.

## Execution
To initiate the primary simulation pipeline:
1. Ensure all folders (`/algorithms`, `/simulations`, `/src`) are added to your MATLAB path.
2. Run `main_GIBHRA_simulation.m` or `master.m` from the root directory.

## Citation
If you utilize this code or algorithm in your research, please cite the following publication:
> D. Cruz Santiago and H. M. Ammari, "An Efficient Geometric Inversion Based Approach to Hole Detection and Restoration in Wireless Sensor Networks," *2025 16th IFIP Wireless and Mobile Networking Conference (WMNC)*. DOI: 10.23919/WMNC67099.2025.11299275
