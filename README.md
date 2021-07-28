Continuous LWE Research
=======================

This is a repository housing scripts for research on how to use Continuous LWE as a basis for security of lattice-based cryptographic primitives.

* `mlwr` contains code for a MLWRSign implementation, which upon working will be adapted to have continuous-valued secrets.
* `saber` contains an implementation of Saber and uSaber, as well as a (not working) adjustment of Saber to naively have continuous values, and a modification of it which is correct but insecure (removes modulus on reply from party sending encrypted message).
* `sliced_associativity` contains tests to see if "slicing" a public key into a sum of binary matrices to avoid multiplicative rollover mod q could allow for associativity. This partially resolves some issues with using continuous values in Saber, but unfortunately the main issue that the part of the encrypted message removed by modding by q is multiplied by a continuous-valued secret message, adding (equidistributed mod q) error to the resulting decrypted message. 
