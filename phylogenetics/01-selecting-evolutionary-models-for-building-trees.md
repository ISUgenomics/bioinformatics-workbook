# Selecting Models of Evolution

For reconstructing the phylogenetic tree, we will need to make some assumptions on how DNA or protein sequences change over time. Mathematically, we would expect the 4 nucleotides or the 20 amino acids have equal chances of changing over a fixed course of time (substitution), biologically this is not the case. Many processes prevents this from happening. As a result there are many substitution models that have been developed, called evolutionary models, that make assumptions on how DNA/protein sequences change over time.

Evolutionary models try to explain the complex biological phenomenon as a simple mathematically explainable formula. Although phylogentic tree building programs can estimate all the different parameters of the substitution matrix, it is advantageous to provide a simple, relatively good fitting model for the data in hand to begin with. This ensures characterization of the evolutionary process accurate.

Process of selecting the model: we test our data for fit with a battery of models, perform likelihood ratio tests (LRT), and select the best fitting model by doing statistical tests. It is only appropriate to do this LRT testing when the 2 models compared differs by one or couple of parameters. The process is as follows (source: **Posada, D.** (2003). Selecting models of evolution. _The phylogenetic handbook: a practical approach to DNA and protein phylogeny_, 256-282.)


1. Select the competing models: one for the null hypothesis H0 and one for the
alternative hypothesis H1.
2. Estimate the tree and the parameters of the model under the null hypothesis.
3. Use the tree and the estimated parameters to simulate 200–1000 replicate data
sets of the same size as the original.
4. For each simulated data set, estimate a tree and calculate its likelihood under
the models representing H0 and H1 (L0 and L1, respectively). Calculate the LRT
statistic &#916; = 2 (loge L1 – loge L0). These simulated &#916;s form the distribution of
the LRT statistic if the null hypothesis was true (i.e., they constitute the null
distribution of the LRT statistic).
5. The probability of observing the LRT statistic from the original data set if the
null hypothesis is true is the number of simulated &#916;s bigger than the original &#916;,
divided by the total number of simulated data sets. If this probability is smaller
than a predefined value (usually 0.05), H0 is rejected
