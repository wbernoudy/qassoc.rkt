# qassoc.rkt
An implementation of quantum associative memory based [quetzal](https://github.com/rhyzomatic/quetzal), a quantum computer simulator.
## Usage
```
> (define P '((0 0 0 0) (0 0 1 1) (0 1 1 0) (1 0 0 1) (1 1 0 0) (1 1 1 1)))
> (learn P)
> (Grover-part '(0 0 1 ?) P)
> (measure-register)
The most likely result is |0011> with a probability of 0.9401041666666672
```
