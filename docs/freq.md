# Allele frequencies of SNPs significantly associated with COVID-19

## Summary
Here, I look at the allele frequencies of the genetic variants significantly
associated (\\( p < 5 \times 10^{-8}\\) ) with COVID-19 in diverse populations,
to assess the potential of replicability of the association in different
populations. Please note that additional allele frequency information is
also available in the published GWAS summary statistics data files.

## Results

Here, the `A1` allele is the risk increasing allele.

### Analysis of ANA_A2_V2
(very severe respiratory confirmed covid vs. population)

### Analysis of ANA_B1_V2
(hospitalized covid vs. not hospitalized covid)

### Analysis of ANA_B2_V2
(hospitalized covid vs. population)

<div class="table-wrapper" markdown="block">



</div>


### Analysis of ANA_C1_V2
(covid vs. lab/self-reported negative)

### Analysis of ANA_C2_V2
(covid vs. population)

### Analysis of ANA_D1_V2
(predicted covid from self-reported symptoms vs. predicted or self-reported non-covid)

## Methods
I downloaded 6 versions of COVID-19 GWAS summary statistics data from
[COVID-19 Host Genetics Initiative](https://www.covid19hg.org/).

* `ANA_A2_V2`: very severe respiratory confirmed covid vs. population
* `ANA_B1_V2`: hospitalized covid vs. not hospitalized covid
* `ANA_B2_V2`: hospitalized covid vs. population
* `ANA_C1_V2`: covid vs. lab/self-reported negative
* `ANA_C2_V2`: covid vs. population
* `ANA_D1_V2`: predicted covid from self-reported symptoms vs. predicted or self-reported non-covid

I used 1000 Genomes Reference panel to obtain allele frequencies in 5 different
populations.

* `AFR`: African
* `AMR`: Admixed American
* `EAS`: East Asian
* `EUR`: European
* `SAS`: South Asian
