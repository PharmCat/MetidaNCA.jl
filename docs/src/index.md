# MetidaNCA

```@meta
CurrentModule = MetidaNCA
```

Non-compartment pharmacokinetic analysis (NCA).

*This program comes with absolutely no warranty. No liability is accepted for any loss and risk to public health resulting from use of this software.

## NCA

Pharmacokinetics, sometimes abbreviated as PK, is a branch of pharmacology dedicated to determine the fate of substances administered to a living organism.

When analyzing pharmacokinetic data, one generally employs either model fitting using nonlinear regression analysis or non-compartmental analysis techniques (NCA). The method one actually employs depends on what is required from the analysis. If the primary requirement is to determine the degree of exposure following administration of a drug (such as AUC), and perhaps the drug's associated pharmacokinetic parameters, such as clearance, elimination half-life, T (max), C (max), etc., then NCA is generally the preferred methodology to use in that it requires fewer assumptions than model-based approaches.

## Validation

Validation report: [validation_report.pdf](./validation_report.pdf)
Appendix 2: [Appendix2.1.pdf](./pdf/Appendix2.1.pdf), [Appendix2.2.pdf](./pdf/Appendix2.2.pdf), [Appendix2.3.pdf](./pdf/Appendix2.3.pdf)

## Contents

```@contents
Pages = [
        "index.md",
        "examples.md",
        "parameters.md",
        "api.md",
      ]
Depth = 3
```

## Reference

* Makoid C, Vuchetich J, Banakar V (1996-1999), Basic Pharmacokinetics;
* Gabrielsson and Weiner (1997), Pharmacokinetic and Pharmacodynamic Data Analysis: Concepts and Applications;
* Gibaldi and Perrier (1982), Pharmacokinetics;
* Wagner (1975), Fundamentals of Clinical Pharmacokinetics.
* Gabrielsson J, Weiner D. Non-compartmental analysis. Methods Mol Biol. 2012;929:377-89. doi: 10.1007/978-1-62703-050-2_16. PMID: 23007438.
