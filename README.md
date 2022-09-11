# BCMF

## Notes on Code:

- Need to think about fitting of model used for extracting clever covariates
  for the ordinal model. Right now we fit a Gaussian model to the binary
  expanded dataset, which is probably not _wrong_ (we have lots of flexibility
  in how we do this), but is at least unnatural.

  - Two reasonable options: just fit a model to the non-expanded data with a
    Gaussian model, or fit a binary probit to the expanded.