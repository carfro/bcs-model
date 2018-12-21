## Background
One the most fundamental questions of nuclear structure theory is to determine the spectra of atomic nuclei. A common approach to solve this many-body problem is the Hartree-Fock-Bogoliubov (HFB) method to approximate the ground state energy. Albeit being a powerful method, it yields a wave function that breaks the known symmetries of the system. Thus a way to increase the accuracy of the approximation would be to restore said symmetries. One way of doing this is to find the components of the HFB state that obeys the symmetries through the use of projectional techniques. By applying these projections one should in principle also be able to deduce the spectra of excited states.

## Project goal
The thesis project objective is to write a program/code that takes HFB-states as input and through projection outputs the components that obeys the symmetries of angular momentum, parity and particle number.

This specific repo is the beginning of this project and so far contains MO + BCS states together with overlapp pfaffian calculations.

## How to run
The makefile is, in it's current state, set to work for os x 10.14.1. 
