# Factor Analysis on Ancient Damage Data 


Load the Lindo 2016 Ancient Damage data

```r
signature_counts <- get(load("../summary_data/signature-counts-clubbed-Lindo2016.rda"))
```

We briefly explore the data summary.

```r
print(dim(signature_counts))
signature_counts[1:5,1:5]
```

We convert the counts data to log CPM data by the `limma::voom` transformation.

```r

### voom transformation

voom_signature_counts <- t(limma::voom(t(signature_counts))$E);

### track the voom weights 

voom_weights <- t(limma::voom(t(signature_counts))$weights);

```

