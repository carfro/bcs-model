## Background
One of the most fundamental questions of nuclear structure theory is to determine the spectra of atomic nuclei. A common approach to solve this many-body problem is to treat the nuclei as a superfluid; catching both the particle-particle as well as particle-hole correlations. This is often done using the Hartree-Fock-Bogoliubov (HFB) transformation from nucleons to quasi-particles, which gives a better approximate the ground state energy. Albeit being a powerful method, it yields a wave function that breaks the known symmetries of the system. Thus a way to increase the accuracy of the approximation would be to restore said symmetries. One way of doing this is to filter out the components of the HFB state that obeys the symmetries through the use of projectional techniques. By applying these projections one should in principle also be able to deduce the spectra of excited, rotational bands.

## Project goal
The thesis project objective is to write a program/code that takes HFB-states as input and through projection outputs the components that obeys the symmetries of angular momentum, parity and particle number.

This specific repo is the beginning of this project and so far contains MO + BCS states together with overlapp pfaffian calculations.

## How to run
Run the project using '$make', to clean the dir use '$make clean'.
The makefile is, in it's current state, set to work for os x 10.14.1. 

