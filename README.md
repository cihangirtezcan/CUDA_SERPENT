# CUDA_SERPENT

This CUDA Optimization of **SERPENT** is used in the ToSC publication _Cryptanalysis: Theory versus Practice - Correcting Cryptanalysis Results on ASCON, ChaCha, and SERPENT using GPUs_ by Cihangir Tezcan, Gregor Leander, and Hosein Hadipour.

You can use benchmark to see how many SERPENT encryptions your GPU can perform in a second. These codes allowed us to perform 2^{34.70} SERPENT encryptions (including the key schedule) per second on an RTX 4090. 

We used these codes to experimentally verify the theoretically obtained differential-linear distinguishers for SERPENT and to find better distinguishers.

Random 100 keys are stored in "keys.txt". 

Number of data of an experiment can be increased after an experiment is completed by storing the bias results in the "startingpoint.txt" file.
