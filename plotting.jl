#=
# Plotting for the paper

In this file I am going to load the file for the fidelity and all the different sensitivities and plot them.
The data is stored into two different kind of files: one containing only the sensitivities against particle loss, and the other one containing all the rest.

The name of the first type of file is `./data/pl_<n>.csv`, where `<n>` is the number of particles.
The name of the second type of file is `./data/fid_<n>.csv`, where `<n>` is the number of particles.

The columns for the first type of file are:

1. pl_esta: particle loss of the esta protocol for 5 corrections applied to the original Hamiltonian
2. pl_esta1: partilcle loss of the ESTA protocol with only one correction applied to the original Hamiltonian
3. pl_esta_cont: particle loss of the ESTA protocol with 5 corrections applied to the continous approximation of the original Hamiltonian
4. pl_sta: particle loss of the STA protocol

The columns for the second type of file are:

1. fid_esta: fidelity of the ESTA protocol
2. fid_esta1: fidelity of the ESTA protocol with only one correction
3. fid_esta_cont: fidelity of the ESTA protocol with continuous corrections
4. fid_sta: fidelity of the STA protocol
5. fid_ad: fidelity of the adiabatic protocol
6. tn_esta: sensitivity of the ESTA protocol
7. tn_sta: sensitivity of the STA protocol
8. tn_esta1: sensitivity of the ESTA protocol with only one correction
9. tn_esta_cont: sensitivity of the ESTA protocol with continuous corrections
10. mn_esta: sensitivity of the ESTA protocol
11. mn_sta: sensitivity of the STA protocol
12. mn_esta1: sensitivity of the ESTA protocol with only one correction
13. mn_esta_cont: sensitivity of the ESTA protocol with continuous corrections
14. final_times: the same vector passed as argument

I will load the files through the `DataFrames` and `CSV` packages and then plot using `PGFPlotsX`.
=#

