# parameters
q = 2**13
p = 2**10
n = 16  # should be 256

# rings
Zq = Integers(q)
Zp = Integers(p)
Zqx.<x> = Zq[]
Zpx.<x> = Zp[]
Rq = QuotientRing(Zqx, x**256 + 1)
Rp = QuotientRing(Zpx, x**256 + 1)


# function to slice matrices
def slice_mat(a, k=1):
    a_copy = matrix(RR, a.nrows(), a.ncols(), a.list())
    slices = []
    while True:
        current_list = []
        for v in a_copy.list():
            if v >= k:
                current_list.append(k)
            elif v <= -k:
                current_list.append(-k)
            else:
                current_list.append(0.0)
        current = matrix(RR, a.nrows(), a.ncols(), current_list)
        a_copy -= current
        if list(set(current.list())) == [0]:
            slices.append(a_copy)
            break
        slices.append(current)

    return slices


def test_associativity():
    # public keys
    A = matrix(Zq, 2 * n, 2 * n, [randint(0, 8) for _ in range(4 * n**2)])
    # A_slices = slice_mat(A)
    A_slices = [matrix(RR, 2 * n, 2 * n, A.list())]

    # private keys
    s_real = matrix(RR, 2 * n, n, [uniform(-2, 1) for _ in range(2 * n**2)])
    s_prime = matrix(Zq, n, 2 * n, [Zq(randint(-2, 1)) for _ in range(2 * n**2)])
    s_prime_real = matrix(RR, n, 2 * n, [RR(v) % q for v in s_prime.list()])

    # compute left multiply
    left_real = matrix(RR, n, 2 * n, [RR(v) % q for v in (s_prime * A).list()])

    # slice left multiply
    # left_slices = slice_mat(left_real)
    left_slices = [left_real]

    # compute left and right rounds
    left_rounds = [matrix(Zp, n, 2 * n, [Zp(round(p / q * round(RR(v) % q)))
                                         for v in left_slice.list()])
                   for left_slice in left_slices]

    right_rounds = [matrix(Zp, 2 * n, n, [Zp(round(p / q * round(RR(v) % q)))
                                          for v in (A_slice * s_real % q).list()])
                    for A_slice in A_slices]

    # compute s' mod p
    s_prime_p = matrix(Zp, n, 2 * n, [Zp(v) for v in s_prime.list()])

    # convert left round to real
    left_rounds_real = [matrix(RR, n, 2 * n, [RR(v) % p for v in left_round.list()])
                        for left_round in left_rounds]

    # compute left and right associativity
    left_assoc = s_prime_p * sum(right_rounds)
    right_assoc = sum([matrix(Zp, n, n, [Zp(round(v)) for v in (left_round_real * s_real % p).list()])
                       for left_round_real in left_rounds_real])

    # compute difference, see if valid
    diffs = [int(v) for v in (left_assoc - right_assoc).list()]
    diffs_valid = all([v < p / 4 or v > 3 * p / 4 for v in diffs])

    # output answer
    print('passed' if diffs_valid else 'failed')


def show_old_steps():
    # public keys
    A = matrix(Zq, 2 * n, 2 * n, [randint(0, 256) for _ in range(4 * n**2)])  # 256 -> q for final
    A_real = matrix(RR, 2 * n, 2 * n, A.list())
    A_slices = slice_mat(A)

    # private keys
    s_real = matrix(RR, 2 * n, n, [uniform(-2, 1) for _ in range(2 * n**2)])
    s_prime = matrix(Zq, n, 2 * n, [Zq(randint(-2, 1)) for _ in range(2 * n**2)])
    s_prime_real = matrix(RR, n, 2 * n, [RR(v) % q for v in s_prime.list()])
    s_prime_p = matrix(Zp, n, 2 * n, [Zp(v) for v in s_prime.list()])

    # right assoc
    r0 = s_prime_real * (p/q * sum([(A_slice * s_real).apply_map(round)
                                 for A_slice in A_slices])).apply_map(round)

    r1 = s_prime_real * (p/q * sum([(A_slice * s_real)
                                 for A_slice in A_slices])).apply_map(round)

    r2 = s_prime_real * (p/q * (A_real * s_real)).apply_map(round)

    r3 = (p/q * (s_prime_real * A_real * s_real))

    # left assoc
    print('computing slices')
    left_slices = slice_mat(s_prime * A)
    l0 = (p/q * sum([(left_slice * s_real).apply_map(round)
                     for left_slice in left_slices])).apply_map(round)

    l1 = (p/q * sum([left_slice * s_real
                     for left_slice in left_slices])).apply_map(round)

    l2 = (p/q * (s_prime_real * A_real * s_real)).apply_map(round)


    l3 = (p/q * (s_prime_real * A_real * s_real))

    diffs = (l0 - r0).list()
    diffs_okay = all([int(d) < p/4 or int(d) > 3 * p/4 for d in diffs])

    print(r0)
    print('---')
    print(r1)
    print('---')
    print(r2)
    print('---')
    print(r3)
    print('---')
    print('---')
    print(l0)
    print('---')
    print(l1)
    print('---')
    print(l2)
    print('---')
    print(l3)
    print('---')
    print('---')
    print(diffs)
    print('---')
    print('success' if diffs_okay else 'failure')


def show_steps():
    # public keys
    A = matrix(RR, 2 * n, 2 * n, [randint(0, 1024) for _ in range(4 * n**2)])  # 16 -> q

    # private keys
    s = matrix(RR, 2 * n, n, [uniform(-2, 1) for _ in range(2 * n**2)])
    s_prime = matrix(RR, n, 2 * n, [randint(-2, 1) for _ in range(2 * n**2)])

    # left side
    A_slices = slice_mat(A)
    l0 = ((s_prime * (p/q * (sum([(A_slice * s).apply_map(round) for A_slice in A_slices]) % q)).apply_map(round)) % p).apply_map(round)

    # left errors
    e0l = (p/q * (sum([(A_slice * s).apply_map(round) for A_slice in A_slices]) % q)).apply_map(round) - p/q * (sum([(A_slice * s).apply_map(round) for A_slice in A_slices]) % q)
    e1l = [(A_slices[j] * s).apply_map(round) - A_slices[j] * s for j in range(len(A_slices))]

    # left final
    l1 = ((p/q * s_prime * A * s + (p/q * s_prime * sum(e1l) + s_prime * e0l)) % p).apply_map(round)

    # right side
    spA_slices = slice_mat((p/q * (s_prime * A % q)).apply_map(round))
    r0 = (matrix(RR, sum([(spA_slice * s).apply_map(round) for spA_slice in spA_slices])) % p).apply_map(round)

    # right errors
    e1r = [(spA_slice * s).apply_map(round) - spA_slice * s for spA_slice in spA_slices]
    e2r = (p/q * (s_prime * A % q)).apply_map(round) - p/q * (s_prime * A % q)
    e3r = (s_prime * A % q) - s_prime * A

    # right final
    r1 = ((p/q * s_prime * A * s + (p/q * e3r * s + e2r * s + sum(e1r))) % p).apply_map(round)

    print('Left real:')
    print(l0)
    print()
    print('Left ideal + errors:')
    print(l1)
    print()
    print('Right real:')
    print(r0)
    print()
    print('Right ideal + errors:')
    print(r1)
    print()

    '''
    print('e0l:')
    print(e0l)
    print()
    print('e1l:')
    print(e1l)
    print()
    print('e1r:')
    print(e1r)
    print()
    print('e2r:')
    print(e2r)
    print()
    print('e3r:')
    print(e3r)
    print()
    '''

    for k in range(100 * q, 200 * q):
        print(max(((0.01 * k) * s % q).list()))


show_steps()
