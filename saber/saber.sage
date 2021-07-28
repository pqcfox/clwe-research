# lightsaber params
l = 2
n = 256
eq = 13
ep = 10
eT = 3
q = 2**eq
p = 2**ep
T = 2**eT
mu = 10

# generates Rq from paper
def gen_ring(q, n):
    Zqx.<x> = Integers(q)[]
    Rq = QuotientRing(Zqx, x**n + 1)
    return Rq

Rq = gen_ring(q, n)
Rp = gen_ring(p, n)
RT = gen_ring(T, n)
R2 = gen_ring(2, n)
m = R2([randint(0, 1) for _ in range(n)])

# generates random matrix in a ring
def random_matrix(k, l, R, n, high=q-1):
    entries = [[R([randint(0, high) for _ in range(n)])
                for _ in range(l)]
               for _ in range(k)]
    return matrix(entries)

# generates random centered binomial secret
def gen_secret(mu, k, l, R, n):
    s = matrix(R, k, l)
    for _ in range(mu // 2):
        s += random_matrix(k, l, R, n, high=1)
        s -= random_matrix(k, l, R, n, high=1)
    return s

# key generation
A = random_matrix(l, l, Rq, n)
s = gen_secret(mu, l, 1, Rq, n)
h1 = Rq([2**(eq - ep - 1) for _ in range(n)])
h2 = Rq([2**(ep - 2) - 2**(ep - eT - 1) + 2**(eq - ep - 1) for _ in range(n)])
h = matrix([[h1] for _ in range(l)])

def shift(b, e, R):
    entries = b.list()
    Z = R.polynomial_ring().base_ring()
    entries_shift = [R([Z(c >> e) for c in p.list()]) for p in entries]
    b_shift = matrix(R, b.nrows(), b.ncols(), entries_shift)
    return b_shift

def shift_poly(b, e, R):
    return shift(matrix([[b]]), e, R)[0, 0]

def coerce(b, R):
    return shift(b, 0, R)

def coerce_poly(b, R):
    return shift_poly(b, 0, R)

b = shift(A.T * s + h, eq - ep, Rp)

# encryption
s_prime = gen_secret(mu, l, 1, Rq, n)
b_prime = shift(A * s_prime + h, eq - ep, Rp)
v_prime = (b.T * coerce(s_prime, Rp))[0, 0]
cm = shift_poly(v_prime + coerce_poly(h1, Rp) - 2**(ep - 1) * coerce_poly(m, Rp), ep - eT, RT)

# decryption
v = (b_prime.T * coerce(s, Rp))[0, 0]
m_prime = shift_poly(v - 2**(ep - eT) * coerce_poly(cm, Rp) + coerce_poly(h2, Rp), ep - 1, R2)

print(m == m_prime)
