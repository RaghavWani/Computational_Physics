# Computational Physics Repository

This repository contains Fortran codes for various short projects. These projects explore various physics problems that can be solved using computational techniques. Here is a short description of each project. More specific details can be found in the README file of each project folder. 

## 1. Molecular Dynamics:

Molecular Dynamics (MD) is a powerful computational simulation technique used to predict the time evolution of interacting particles by integrating their equations of motion. It's a specific N-body simulation focusing on the dynamics and interactions of atoms or molecules within a system. In this MD simulation, trajectories of particles are determined by numerically solving Newton's equations of motion, where forces and potential energies are calculated using **Lennard-Jones potential**.

##### **Verlet Algorithm**

This project features MD simulations utilizing the Verlet algorithm, a stable and computationally efficient numerical integration method. Widely employed in condensed matter physics, materials science, and biophysics, the Verlet algorithm accurately predicts the behavior of molecular systems over time. By implementing the Verlet algorithm, this repository facilitates the exploration of complex molecular dynamics, enabling researchers to study fundamental physical and chemical phenomena at the atomic scale.

## 2. Ising Model:

The Ising model is a mathematical model in statistical mechanics used to study the behavior of magnetic systems. It consists of a lattice of spins, where each spin can be in one of two states, typically representing the magnetic orientation of atoms. The Ising Hamiltonian describes the interactions between neighboring spins, which determines the system's total energy based on the spin configurations. 

##### **Monte Carlo simulations**
Here, I have employed Monte Carlo simulations to explore the behavior of the Ising model by randomly sampling configurations and evaluating their energies. This stochastic approach allows us to simulate the system's evolution over time and study various thermodynamic properties. By leveraging Monte Carlo methods, one can investigate phase transitions, critical phenomena, and other thermodynamic properties of magnetic systems.

## 3. Turing Patterns:

Turing Patterns arise from the instability of homogeneous steady states in reaction-diffusion systems. The formation of spatial patterns is a phenomenon observed in various biological, chemical, and physical systems. This project focuses on simulating Turing Patterns using the finite difference method to solve partial differential equations, particularly the **Fitzhugh-Nagumo equation**. By incorporating diffusion terms into the model, one can observe the emergence of intricate patterns and study the underlying mechanisms driving pattern formation in dynamic systems. 
