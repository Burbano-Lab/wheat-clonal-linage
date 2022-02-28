# A pandemic clonal lineage of the wheat blast fungus
# 7. AVR Rmg8 analyses

Program              | Location
-------------------- | -----------------------------------
*MUSCLE v.3.5.1551*  | (https://github.com/rcedgar/muscle)
*RAxML-NG v.1.0.3*   | (https://github.com/amkozlov/raxml-ng)


We created the multi-fasta file [M_oryzae_AVR-Rmg8_blastresults_curated.fna](/data/07_AVR_Rmg8/M_oryzae_AVR-Rmg8_blastresults_curated.fna) containing animoacid sequences of the AVR Rmg8 protein in various isolates. Then, we used *MUSCLE* to perform a multiple sequence alignment.

```bash
muscle -in M_oryzae_AVR-Rmg8_blastresults_curated.fna -out M_oryzae_AVR-Rmg8_blastresults_curated.aligned.fasta
```

---
[Main README](/README.md) | [Previous - 06. Mating Type Analyses](/06_Mating_Type.md)
