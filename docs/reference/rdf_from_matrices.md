# Fast approximate residual degrees of freedom

Defining \\H = A^TA + B^TB\\ where \\A\\ and \\B\\ are low rank, compute
\\n - 2tr(H) + tr(HH)\\ in \\O(np^2)\\ instead of \\O(n^2p^2)\\.

## Usage

``` r
rdf_from_matrices(A, B)
```

## Arguments

- A:

  a `matrix` or `sparseMatrix`

- B:

  a `matrix` or `sparseMatrix`

## See also

rdf.merMod
