# Protein Complexes

The protein complexes are taken from Complex Portal for _Homo sapiens_ and expression data (FPKM) for basemine hepaic cell models from Hecatos Project. The expression data is from multiple timepoints.

1. The **complex expression** is computed by first computing the expression of each subunit and then the lowest expressed subunit defines the expression of the subunit.
2. The prot_complex_exp_order.r employs dynamic Bayesian networks to expression data to generate networks. And from these networks **order of assembly** is predicted

