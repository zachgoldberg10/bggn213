Class 11: Structural Bioinformatics 1
================

### Structural Bioinformatics

## Section 1: Introduction to the RCSB Protein Data Bank (PDB)

Q1: Download a CSV file from the PDB site (accessible from “Analyze” -\>
“PDB Statistics” \> “by Experimental Method and Molecular Type”. Move
this CSV file into your RStudio project and determine the percentage of
structures solved by X-Ray and Electron Microscopy. Also can you
determine what proportion of structures are protein?

``` r
data <- read.csv("Data Export Summary.csv")
data
```

    ##   Experimental.Method Proteins Nucleic.Acids Protein.NA.Complex Other  Total
    ## 1               X-Ray   131463          2060               6768     8 140299
    ## 2                 NMR    11241          1304                262     8  12815
    ## 3 Electron Microscopy     2925            32               1004     0   3961
    ## 4               Other      280             4                  6    13    303
    ## 5        Multi Method      144             5                  2     1    152

``` r
ans <- data$Total/sum(data$Total) *100
names(ans) <- data$Experimental.Method
round(ans,2)
```

    ##               X-Ray                 NMR Electron Microscopy               Other 
    ##               89.06                8.13                2.51                0.19 
    ##        Multi Method 
    ##                0.10

Also can you determine what proportion of structures are protein?

``` r
areProt <- sum(data$Proteins)/sum(data$Total)
round(areProt, 2)
```

    ## [1] 0.93

## Section 3: Bio3d

``` r
library(bio3d)
pdb <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
str(pdb)
```

    ## List of 8
    ##  $ atom  :'data.frame':  1686 obs. of  16 variables:
    ##   ..$ type  : chr [1:1686] "ATOM" "ATOM" "ATOM" "ATOM" ...
    ##   ..$ eleno : int [1:1686] 1 2 3 4 5 6 7 8 9 10 ...
    ##   ..$ elety : chr [1:1686] "N" "CA" "C" "O" ...
    ##   ..$ alt   : chr [1:1686] NA NA NA NA ...
    ##   ..$ resid : chr [1:1686] "PRO" "PRO" "PRO" "PRO" ...
    ##   ..$ chain : chr [1:1686] "A" "A" "A" "A" ...
    ##   ..$ resno : int [1:1686] 1 1 1 1 1 1 1 2 2 2 ...
    ##   ..$ insert: chr [1:1686] NA NA NA NA ...
    ##   ..$ x     : num [1:1686] 29.4 30.3 29.8 28.6 30.5 ...
    ##   ..$ y     : num [1:1686] 39.7 38.7 38.1 38.3 37.5 ...
    ##   ..$ z     : num [1:1686] 5.86 5.32 4.02 3.68 6.34 ...
    ##   ..$ o     : num [1:1686] 1 1 1 1 1 1 1 1 1 1 ...
    ##   ..$ b     : num [1:1686] 38.1 40.6 42.6 43.4 37.9 ...
    ##   ..$ segid : chr [1:1686] NA NA NA NA ...
    ##   ..$ elesy : chr [1:1686] "N" "C" "C" "O" ...
    ##   ..$ charge: chr [1:1686] NA NA NA NA ...
    ##  $ xyz   : 'xyz' num [1, 1:5058] 29.36 39.69 5.86 30.31 38.66 ...
    ##  $ seqres: Named chr [1:198] "PRO" "GLN" "ILE" "THR" ...
    ##   ..- attr(*, "names")= chr [1:198] "A" "A" "A" "A" ...
    ##  $ helix :List of 4
    ##   ..$ start: Named num [1:2] 87 87
    ##   .. ..- attr(*, "names")= chr [1:2] "" ""
    ##   ..$ end  : Named num [1:2] 90 90
    ##   .. ..- attr(*, "names")= chr [1:2] "" ""
    ##   ..$ chain: chr [1:2] "A" "B"
    ##   ..$ type : chr [1:2] "1" "1"
    ##  $ sheet :List of 4
    ##   ..$ start: Named num [1:17] 10 18 32 75 52 43 62 69 96 96 ...
    ##   .. ..- attr(*, "names")= chr [1:17] "" "" "" "" ...
    ##   ..$ end  : Named num [1:17] 15 23 34 78 59 49 66 73 98 98 ...
    ##   .. ..- attr(*, "names")= chr [1:17] "" "" "" "" ...
    ##   ..$ chain: chr [1:17] "A" "A" "A" "A" ...
    ##   ..$ sense: chr [1:17] "0" "-1" "0" "1" ...
    ##  $ calpha: logi [1:1686] FALSE TRUE FALSE FALSE FALSE FALSE ...
    ##  $ remark:List of 1
    ##   ..$ biomat:List of 4
    ##   .. ..$ num   : int 1
    ##   .. ..$ chain :List of 1
    ##   .. .. ..$ : chr [1:2] "A" "B"
    ##   .. ..$ mat   :List of 1
    ##   .. .. ..$ :List of 1
    ##   .. .. .. ..$ A B: num [1:3, 1:4] 1 0 0 0 1 0 0 0 1 0 ...
    ##   .. ..$ method: chr "AUTHOR"
    ##  $ call  : language read.pdb(file = "1hsg")
    ##  - attr(*, "class")= chr [1:2] "pdb" "sse"

Atom selection

``` r
ca.inds <- atom.select(pdb, "calpha")
```

Write new PDB of only protein

``` r
prot <- atom.select(pdb, "protein", value =  TRUE)
write.pdb(prot, file = "1hsg_protein.pdb")
```

Write new PDB of only ligand

``` r
ligand <- atom.select(pdb, "ligand", value = TRUE)
write.pdb(ligand, file = "1hsg_ligand.pdb")
```
