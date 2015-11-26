# qassoc.rkt
An implementation of quantum associative memory based on the quantum computer simulator quetzal.rkt.
## Usage
```
> (define P '((0 0 0 0) (0 0 1 1) (0 1 1 0) (1 0 0 1) (1 1 0 0) (1 1 1 1)))
> (learn P)
> (Grover-part '(0 0 1 ?) P)
> (measure-register)
The most likely result is |0011> with a probability of 0.9401041666666672
```
