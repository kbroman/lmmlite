## Comparisons to [pylmm](http://github.com/nickFurlotte/pylmm)

The scripts in this directory are for comparing
[lmmlite](https://github.com/kbroman/lmmlite) to
[pylmm](http://github.com/nickFurlotte/pylmm) (with the example data
included in lmmlite), looking both at REML and maximum likelihood.

Basic findings:
  - There's one phenotype where pylmm gave an estimated heritability
    of 0.99 (firm), while lmmlite gave an estimate closer to 1.
  - The beta estimates are essentially the same except that one case.
    with pylmm giving heritability of 0.99.
  - The lmmlite formula for the log likelihood ignores some terms that depend
    only on the sample size, so the log likelihoods for the two
    programs are not the same. (But they're basically the same up to a
    constant that depends on sample size.)

So it looks like lmmlite is working fine.
