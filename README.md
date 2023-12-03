# QC-MDPC Markovian model

This repository contains implementation code for a model of a simple variant of bit-flipping decoding, named step-by-step decoding. Despite its higher Decoding Failure Rate (DFR), the behavior of this variant can be analyzed as a Markov chain, in line with the theoretical framework from Julia Chaulet's PhD thesis[^1].

More details can be found in the corresponding paper
> Nicolas Sendrier & Valentin Vasseur: On the Decoding Failure Rate of QC-MDPC Bit-Flipping Decoders. In: PQCrypto 2019. LNCS, vol 11505. Springer, Cham (July 2019). <https://doi.org/10.1007/978-3-030-25510-7_22>

[^1]: <https://theses.hal.science/tel-01599347>


## Usage

```sh
./model r d t_pass t_fail [noX|XSt|Xt] [uniform|single|pair|triple]

Arguments:
  r                     - Block length of the QC-MDPC code
  d                     - Block weight of the QC-MDPC code
  t_pass                - Threshold value; decoder is assumed to always succeed below this
  t_fail                - Threshold value; decoder is assumed to always fail above this
  [noX|XSt|Xt]          - Counters modelling choices:
                           noX: assume X is always 0
                           XSt: approximate X depending on S and t
                           Xt: only consider t
  [uniform|single|pair|triple]
                        - Position choosing strategies:
                           uniform: choose a position uniformly at random
                           single: choose a random position involved in a random unverified equation
                           pair: choose a random position involved in two random unverified equations
                           triple: choose a random position involved in three random unverified equations
```


## Example

```sh
$ ./model 12323 71 5 150 XSt single
[...]
134 8686 1.5584092502450608068002708549156985159297998811537e-38
134 8688 1.5758903267089172240645171161458487415017988782202e-38
134 8690 1.673091742453008705786241127347703254834452558173e-38
134 8692 5.9145902060027073480858907898187387746980615426662e-39
134 8694 4.7131401367262487874693617850630851440251814087452e-39
134 8696 2.6785897848013041846947680358438709992965071236336e-39
134 8698 2.1395391582221008492292814688905799096378425267951e-39
134 8700 2.147109004254025258393609940596833102914633622043e-39
134 8702 1.9410696871843215182015955838967805211324456241996e-39
134 8704 4.7886006688963440857226427460464777837970945068226e-40
[...]
```

In this example, the model predicts that with an error weight of 134 and a syndrome weight of 8686, the probability of decoding failure using the step-by-step decoder is approximately 1.56e-38.
