# Infinite-Nuclear-Matter
Simple Monte Carlo simulation to model momentum distributions in an infinite nucleus
With this code, we show that we need to enforce Q = 0, in order to avoid contaminations from non-SRC pairs in the 2-body momentum
distribution of nucleons in infinite nuclear matter. We create a simulation where we describe infinite nuclear matter as a correlated Fermi gas. From this simulation, we generate 1-body momentum distributions. Then, we use these 1-body momentum distributions to construct 2-body momentum distributions. The 2-body momentum distributions constitute a handy tool to study 2-body interactions. Thus, we draw some fundamental conclusions about the nature of SRC pairs from these distributions. We begin by modeling the nucleus as a correlated Fermi gas. This means that the momentum occupation distribution of the nucleons is divided into two regions, as shown here:

* Fermi gas region, 0 < k < kf
* Correlated high-momentum tail, kf < k < 5 fm−1

where kf is the Fermi momentum. The Fermi gas region, which is a good representation for the nuclear mean-field potential, ranges between 0 < k < kf . The correlated high-momentum tail ranges between kf < k < 5 fm−1 and quickly decays to zero as k^(−4). The mathematical expression of this distribution is as follows:

               | A0    , 0 < k < kf
     nCFG(k) = | C/k^4 , kf < k < 5 fm−1
               | 0     , otherwise
          
In this expression, A0 and C are normalization constants. C is called the contact and is related to the probability to find two nucleons close to each other. The value of kf = 250 MeV/c = 1.27 fm^(−1) was experimentally determined by E. J. Moniz et al.
In our study, we construct the 1-body momentum distribution by generating a large number of values that follow the distribution above. Experiments have shown that SRC pairs account for approximately 20% of all the nucleons in the nucleus. Thus, in our simulation, 80% of all the nucleons are created in the Fermi gas region, and 10% in the correlated high momentum tail. For each nucleon belonging to this 10%, a correlated partner is created. This accounts for the remaining 10%. From these conditions we obtain the normalization.

The process of creating the 1-body momentum distribution consists of integrating out all the angular dependence from the momentum occupation distribution.

As stated before, a good way to study 2-body systems is by analyzing their 2-body momentum distribution. We create these distributions by classifying every pair according to one of the following categories:

• mean field-mean field (MF-MF): a pair falls into this category if both nucleons have a 1-body momentum that belongs to the mean-field region (i.e. k < kf ).

• mean field-SRC (MF-SRC): a pair falls into this category if one nucleon has a 1-body momentum that belongs to the mean-field region (i.e. k < kf ) and the other nucleon has a 1-body momentum that belongs to the high-momentum tail (i.e. k > kf ).

• SRC-SRC: a pair falls into this category if both nucleons have a 1-body momentum that belongs to the high-momentum tail (i.e. k > kf ).

• True SRC: pairs in this category correspond the the pairs that were originally defined as SRC-SRC in the simulation.

Therefore, we loop over every possible nucleon pair and classify it according
to one of these categories.
