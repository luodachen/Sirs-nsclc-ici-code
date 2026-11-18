# deconstructSigs 1.9.0 -- hg38 patch note

**Source**: https://github.com/raerose01/deconstructSigs (master, Version 1.9.0; archived on CRAN)
**Obtained via**: the mirror https://gh-proxy.com/ (direct GitHub access was unreliable from the
build machine)

## Why the patch is needed

1. CRAN only keeps releases up to 1.8.0, and **1.8.0 hard-codes hg19 in the source** (there is no
   `bsg` argument), so it cannot be used with hg38 coordinates.
2. 1.9.0 adds the `bsg` argument for other genomes, but its DESCRIPTION still **hard-depends on
   `BSgenome.Hsapiens.UCSC.hg19`** (another ~700 MB package). A plain install therefore fails with:
   `ERROR: dependency 'BSgenome.Hsapiens.UCSC.hg19' is not available`

## Patch contents (2 items, plain text substitutions)

| # | File | Change |
|---|---|---|
| 1 | `DESCRIPTION` | remove `BSgenome.Hsapiens.UCSC.hg19` from Imports (reshape2 / BSgenome / GenomeInfoDb / grDevices / graphics / utils are kept) |
| 2 | `R/*.R` (7 places, mainly `mut.to.sigs.input.R` and `vcf.to.sigs.input.R`) | `BSgenome.Hsapiens.UCSC.hg19` -> `BSgenome.Hsapiens.UCSC.hg38`, so that the **default branch also uses hg38** (this prevents a silent fall-back to hg19 coordinates when `bsg` is omitted) |

## Installation

```r
install.packages("vendor/deconstructSigs-1.9.0-hg38patched", repos = NULL, type = "source")
# result: deconstructSigs 1.9.0 installed successfully (pure-R package, no Rtools required)
```

## Suggested wording for the manuscript Methods

> Somatic mutational signatures were derived with **deconstructSigs v1.9.0** (Rosenthal et al., 2016),
> using the 96-trinucleotide context computed from the **UCSC hg38** reference genome
> (`BSgenome.Hsapiens.UCSC.hg38`) via the `bsg` argument, and the COSMIC v2 30-signature set
> (`signatures.cosmic`) with exome-to-genome trinucleotide-bias correction.

**Note**: if a reviewer asks for the COSMIC v3.2 SBS nomenclature (SBS4 = tobacco smoking), the
analysis can be repeated with the Bioconductor package `MutationalPatterns` (native hg38 support,
no patch needed); the two agree on the smoking-related signature call.
