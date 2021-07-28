# generates Rq from paper
def gen_ring(q, n):
    Zqx.<x> = Integers(q)[]
    Rq = QuotientRing(Zqx, x**n + 1)
    return Rq

# generates random matrix in a ring
def random_matrix(k, l, R, n, q=None):
    if q is None:
        Zq = R.polynomial_ring().base_ring()
        q = Zq.order()
    entries = [[R([randint(0, q) for _ in range(n)])
                for _ in range(l)]
               for _ in range(k)]
    return matrix(entries)

# generates random centered binomial secret
def gen_secret(u, k, l, R, n):
    entries = [[R([randint(-2**(u-1), 2**(u-1) - 1) for _ in range(n)])
                for _ in range(l)]
               for _ in range(k)]
    return matrix(entries)

# generates random centered binomial secret
def gen_binary_secret(k, l, R, n):
    entries = [[R([randint(0, 1)] + [0 for _ in range(n - 1)])
                for _ in range(l)]
               for _ in range(k)]
    return matrix(entries)

# generates random centered binomial secret
def gen_cont_secret(u, k, l, R, n):
    entries = [[R([uniform(-2**(u-1), 2**(u-1) - 1) for _ in range(n)])
                for _ in range(l)]
               for _ in range(k)]
    return matrix(entries)

# bit shift each entry of a matirx
def shift(b, e, R):
    entries = b.list()
    try:
        Z = R.polynomial_ring().base_ring()
    except AttributeError:  # no quotient, i.e. continuous
        Z = R.base_ring()
    entries_shift = [R([Z(floor(c) >> e) for c in p.list()]) for p in entries]
    b_shift = matrix(R, b.nrows(), b.ncols(), entries_shift)
    return b_shift

# bit shift a single polynomial
def shift_poly(b, e, R):
    return shift(matrix([[b]]), e, R)[0, 0]

# coerce each entry of a matrix
def coerce(b, R):
    return shift(b, 0, R)

# coerce a single polynomial
def coerce_poly(b, R):
    return shift_poly(b, 0, R)

# ulightsaber params
l = 2
n = 256
eq = 13
ep = 12
eT = 11
q = 2**eq
p = 2**ep
T = 2**eT
u = 2

# rings in the paper
Rq = gen_ring(q, n)
Rp = gen_ring(p, n)
RT = gen_ring(T, n)
R2 = gen_ring(2, n)
Rc.<x> = RR[]

# random message
m = R2([randint(0, 1) for _ in range(n)])

# key generation
A = random_matrix(l, l, Rq, n)
s = gen_cont_secret(u, l, 1, Rc, n)
h1 = Rq([2**(eq - ep - 1) for _ in range(n)])
h2 = Rq([2**(ep - 2) - 2**(ep - eT - 1) + 2**(eq - ep - 1) for _ in range(n)])
h = matrix([[h1] for _ in range(l)])
b = shift(coerce(coerce(A.T, Rc) * s, Rq) + h, eq - ep, Rp)

# encryption
s_prime = gen_binary_secret(l, 1, Rq, n)

# b_prime = shift(A * s_prime + h, eq - ep, Rp)
# NB: this change makes it closer to Regev
b_prime = coerce(A, Rc) * coerce(s_prime, Rc)

v_prime = (b.T * coerce(s_prime, Rp))[0, 0]
cm = shift_poly(v_prime + coerce_poly(h1, Rp) - 2**(ep - 1) * coerce_poly(m, Rp), ep - eT, RT)

# decryption
v = (coerce(coerce(b_prime, Rc).T * s, Rp))[0, 0]

print(v_prime)
print(v)

m_prime = shift_poly(v - 2**(ep - eT) * coerce_poly(cm, Rp) + coerce_poly(h2, Rp), ep - 1, R2)


# tests
def cont_mod(b, q):
    entries = b.list()
    entries_mod = [Rc([RR(c % q if c % q >= 0 else c % q + q)
                       for c in p.list()]) for p in entries]
    b_mod = matrix(Rc, b.nrows(), b.ncols(), entries_mod)
    return b_mod


def cont_mod_poly(b, q):
    return cont_mod(matrix([[b]]), q)[0, 0]


# check 1: associativity of real multiplication (sanity)
# val1 = coerce(s_prime, Rc).T * (coerce(A, Rc).T * s)
# val2 = (coerce(s_prime, Rc).T * coerce(A, Rc).T) * s
# print('Check 1:')
# print(val1 == val2)

# check 2: associativity of real multiplication mod q
# val1 = s_prime.T * coerce(coerce(A, Rc).T * s, Rq)
# val2 = coerce(coerce(s_prime.T * A.T, Rc) * s, Rq)
# print('Check 2:')
# print(val1 == val2)

# check 3: show the values of coercing
# print('Check 3:')
# print(A[0, 0])
# print('--------')
# print(coerce(A, Rc).T[0, 0])

# print('Check 4:')
# print(A[0, 0].list()[0])
# print('--------')
# print(coerce(A, Rc).T[0, 0].list()[0])

# print(A[0, 0].list()[0])
# print(RR(floor(A[0, 0].list()[0])))


# print(cont_mod_poly(cont_mod_poly(Rc([1, 0, 1, 1]) * Rc([2, 3]), 5) * Rc([4, 0.5]), 5))
# print(cont_mod_poly(Rc([1, 0, 1, 1]) * cont_mod_poly(Rc([2, 3]) * Rc([4, 0.5]), 5), 5))

# check correctness for m
print(m == m_prime)
