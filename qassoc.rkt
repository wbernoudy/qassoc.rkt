#lang racket
(require math)
(require racket/vector)

(require 2htdp/batch-io)
(define-namespace-anchor a)
(define ns (namespace-anchor->namespace a))

(require "quetzal.rkt") ; note that quetzal.rkt must be present in current directory

(define (sq n) (* n n))

(define measure-register2 (λ () ; Gives more info than (measure-register)
    (let ([Ψ (matrix->list register)] [q-index 0] [max 0] [probabilities null])
        (set! probabilities (map (λ (qubit) (magnitude qubit)) Ψ))
        (for ([qubit (length Ψ)])
            (if (> (list-ref probabilities qubit) 0) (printf "~a |~a>, " (sq (list-ref probabilities qubit)) 
                (~r qubit #:base 2 #:min-width (exact-round (/ (log (length Ψ)) (log 2))) #:pad-string "0")) (printf ""))
            (when (< max (list-ref probabilities qubit))
                (set! max (list-ref probabilities qubit))
                (set! q-index qubit)))
        (displayln "")
        (display "The most likely result is |")
        (display (~r q-index #:base 2 #:min-width (exact-round (/ (log (length Ψ)) (log 2))) #:pad-string "0"))
        (display "> with a probability of ") (displayln (* max max))
        )))

;-----------------Start qassoc.rkt--------------------;

(define (S-con p)
    (matrix [
       [1 0 0 0]
       [0 1 0 0]
       [0 0 (sqrt (/ (sub1 p) p)) (/ 1 (sqrt p))]
       [0 0 (/ 1 (sqrt p)) (sqrt (/ (sub1 p) p))]
    ]))

(define (learn patterns)
    (letrec ([bits (length (car patterns))]
        [qubits (+ (length (car patterns)) 3)] ; n qubits for storage, 1 intermediate qubit, 2 control qubits
        [C0NOT (matrix* (G-nqubit-constructor 4 '(0) Pauli-X-gate) ; C0NOT is the same as CNOT surrounded by nots on the control qubit
            CNOT-gate
            (G-nqubit-constructor 4 '(0) Pauli-X-gate))]
        [big-control-X (matrix-stack (list (submatrix (identity-matrix (expt 2 (+ 1 bits))) (- (expt 2 (+ 1 bits)) 2) (expt 2 (+ 1 bits)))
                                         (matrix-row (identity-matrix (expt 2 (+ 1 bits))) (- (expt 2 (+ 1 bits)) 1))
                                         (matrix-row (identity-matrix (expt 2 (+ 1 bits))) (- (expt 2 (+ 1 bits)) 2))))])

        (initialize-register (build-list qubits (λ (x) 0))) ; Initialize all qubits to |0>

        (for ([pattern patterns] [p-index (length patterns)])

            (for ([bit pattern] [b-index bits])
                (if (= bit 1)
                    (apply-gate register (list (sub1 qubits) b-index) C0NOT) ; flip the corresponding qubit if the bit in the pattern is a 1
                    (void)))

            (apply-gate register (list (- qubits 1) (- qubits 2)) C0NOT) ; flip the c1 control qubit
            (apply-gate register (list (- qubits 2) (- qubits 1)) (S-con (- (length patterns) p-index))) ; apply the "save" gate to control qubits

            (for ([bit pattern] [b-index bits])
                (if (= bit 0)
                    (apply-gate register (list b-index) Pauli-X-gate) ; apply not gates on all qubits whose corresponding bit is 0
                    (void)))

            (apply-gate register (build-list (- qubits 2) values) big-control-X) ; apply not gate with controls on all storage qubits on intermediate bit
            (apply-gate register (list (- qubits 3) (- qubits 2)) CNOT-gate) ; CNOT with control on intermediate qubit, target on c1
            (apply-gate register (build-list (- qubits 2) values) big-control-X) ; reverse the controlled not

            (for ([bit pattern] [b-index bits]) ; reverse all the not gates to the storage qubits
                (if (= bit 0)
                    (apply-gate register (list b-index) Pauli-X-gate)
                    (apply-gate register (list (sub1 qubits) b-index) C0NOT))))

        (apply-gate register (list (- qubits 1)) Pauli-X-gate) ; the algorithm will leave all states with the c2 bit as |1>, so not all these

        (remove-unentangled-qubit register) ; We can now ignore the last three qubits as they are not entangled
        (remove-unentangled-qubit register)
        (remove-unentangled-qubit register)))

(define (remove-unentangled-qubit reg) ; removes the last qubit from the system (assumes it is unentangled)
    (let [(new-reg (submatrix reg 1 (build-list (/ (matrix-num-cols reg) 2) (λ (x) (* x 2)))))]
        (set-register new-reg)
        new-reg))


(define (constr-big-control-Z N)
    (diagonal-matrix (build-list N (lambda (x) (if (= x (- N 1)) -1 1)))))

(define (get-not-qubits qubits mtau)   ; gets all the qubits which we know should be flipped
            (cond [(null? mtau) '()]
                [(eq? (car mtau) 0) (cons (- qubits (length mtau)) (get-not-qubits qubits (cdr mtau)))]
                [else (get-not-qubits qubits (cdr mtau))]))

(define (constr-search-phaser tau)
    (letrec ([qubits (length tau)]
        [N (expt 2 qubits)]
        [find-unknown (lambda (mtau)     ; get an unknown to apply a controlled Z gate to
            (if (eq? (car mtau) '?) (- qubits (length mtau))
                (find-unknown (cdr mtau))))]
        [known-bits (filter (lambda (x) (not (eq? (list-ref tau x) '?))) (build-list qubits values))]
        [big-control-Z (constr-big-control-Z (expt 2 (+ (length known-bits) 1)))]
        [to-not-gate (lambda (q)   ; applies an X gate
            (G-nqubit-constructor N (list q) Pauli-X-gate))]
        [to-Z-gate (lambda (q)    ; applies the Z gate to the two possible states of the unknown qubit
            (matrix* (G-nqubit-constructor N (append known-bits (list q)) big-control-Z) ; first a Z on the unknown
                (to-not-gate q) ; flip the unknown
                (G-nqubit-constructor N (append known-bits (list q)) big-control-Z) ; another Z
                (to-not-gate q) ; flip it back to preserve the state
                ))])

        (eval (cons matrix* (append
            (map to-not-gate (get-not-qubits qubits tau)) ; first apply X to all the qubits that should be flipped (so that if we have the right pattern, all qubits are 1)
            (list (to-Z-gate (find-unknown tau))) ; controlled Z gate on on of the unknowns
            (map to-not-gate (reverse (get-not-qubits qubits tau))))) ns) ; apply the same Xs in reverse order to preserve state
        ))

(define (constr-patterns-phaser patterns)
    (letrec ([qubits (length (car patterns))]
        [N (expt 2 qubits)]
        [big-control-Z (constr-big-control-Z N)]
        [to-not-gate (lambda (q)   ; applies an X gate
            (G-nqubit-constructor N (list q) Pauli-X-gate))]
        [constr-helper (lambda (pats)
            (cond [(null? pats) (list (identity-matrix N))]
                [else (append (map to-not-gate (get-not-qubits qubits (car pats)))
                    (list big-control-Z)
                    (map to-not-gate (get-not-qubits qubits (car pats)))
                    (constr-helper (cdr pats)))]))])

        (eval (cons matrix* (constr-helper patterns)) ns)
        ))

(define make-Hadamard (λ (N)
    (let ([constant (exact->inexact (/ 1 (expt 2 (/ (/ (log N) (log 2)) 2))))])
        (build-matrix N N (lambda (i j)
            (* constant (expt -1 (matrix-dot (bits->row-matrix (bits i N)) (bits->row-matrix (bits j N))))))))))

(define (calc-T p r0 r1 N) ; p: total patterns, r0: # of states that don't correspond to patterns, r1: # that do
    (letrec ([a (/ (* 2 (- p (* 2 r1))) N)]
        [b (/ (* 4 (+ p r0)) N)]
        [k (+ (- (* 4 a) (* a b)) (/ r1 (+ r0 r1)))]
        [l (- (/ (* 2 a (+ N (- p r0 (* 2 r1)))) (- N r0 r1)) (* a b) (/ (- p r1) (- N r0 r1)))])

        (exact-round (/ (- (/ pi 2) (atan (/ (* k (sqrt (/ (+ r0 r1) (- N r0 r1)))) l))) (acos (- 1 (/ (* 2 (+ r1 r0)) N)))))))

(define (Grover-part tau patterns)
    (letrec ([qubits (length tau)]
        [N (expt 2 qubits)]
        [all-qubits (build-list qubits values)]
        [full-Hadamard (make-Hadamard N)]
        [Hadamard (lambda (reg)
            (apply-gate register all-qubits full-Hadamard))]
        [search-phase-flip-gates (constr-search-phaser tau)]
        [patterns-phase-flip-gates (constr-patterns-phaser patterns)]
        [T (calc-T (length patterns) (- (sq (count (lambda (x) (eq? x '?)) tau)) 1) 1 N)] ; Assume that there is only one match in the learned patterns
        [big-control-Z (constr-big-control-Z N)])

        (displayln T)

        (apply-gate register all-qubits search-phase-flip-gates) ; Apply the phase flip for the one we're searching for

        (Hadamard (apply-gate (Hadamard register) (reverse all-qubits) big-control-Z)) ; One Grover diffusion

        (apply-gate register all-qubits patterns-phase-flip-gates) ; Apply gate that flips the phase of all patterns

        (Hadamard (apply-gate (Hadamard register) (reverse all-qubits) big-control-Z)) ; A second Grover diffusion

        (for ([i T]) ; We apply the search phase flipper and the Gover diffusion T times
            (apply-gate register all-qubits search-phase-flip-gates) ; Apply the phase flip for the one we're searching for
            (Hadamard (apply-gate (Hadamard register) (reverse all-qubits) big-control-Z)) ; Apply a Grover diffusion
        )))
