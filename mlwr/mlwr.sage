# parameters
q, p = 2**23, 2**19
k, l = 3, 2
gamma1 = 2**19
n = 256
d = 10
eta = 8

# rings
Zq = Integers(q)
Zp = Integers(p)
Zqx.<x> = Zq[]
Zpx.<x> = Zp[]
Rc.<x> = RR[]
Rq = QuotientRing(Zqx, x**256 + 1)
Rp = QuotientRing(Zpx, x**256 + 1)

# main setup
A = matrix(Rq, k, l, [Rq([randint(0, q - 1) for _ in range(256)]) for _ in range(k * l)])
s1 = vector(Rq, l, [Rq([randint(-eta, eta) for _ in range(256)]) for _ in range(l)])
t = vector(Rp, k, [Rp([round(p/q * RR(v)) for v in poly.list()]) for poly in (A * s1).list()])
y = vector(Rq, l, [Rq([randint(-(gamma1 - 1), gamma1 - 1) for _  in range(256)]) for _ in range(l)])
w = vector(Rp, k, [Rp([round(p/q * RR(v)) for v in poly.list()])
                   for poly in (A * y).list()])

c_vec = [0 for _ in range(n)]
for i in sample(range(n), 60):
    c_vec[i] = choice([-1, 1])
c = Rc(c_vec)

z = y + Rq(c.list()) * s1
xi1 = vector(Rc, k, [[(round(p/q * RR(v)) - p/q * RR(v)) for v in poly.list()]
                     for poly in (A * y).list()])
s2 = vector(Rc, k, [Rc([(round(p/q * RR(v)) - p/q * RR(v)) for v in poly.list()])
                    for poly in (A * s1).list()])
xi2 = vector(Rc, k, [[Rc([round(v) - v]) for v in (c * poly).list()]
                     for poly in s2.list()])


nu = vector(Rp, [[round(v - w) for v, w in zip(poly1.list(), poly2.list())]
                 for poly1, poly2 in zip(xi1.list(), xi2.list())])

# associativity test
rndAz = vector(Rp, k, [Rp([round(p/q * RR(v)) for v in poly.list()])
                       for poly in (A * z).list()])
ct = vector(Rp, k, [Rp(c.list()) * poly for poly in t.list()])
cs2r = vector(Rp, k, [Rp([round(v) for v in (c * poly).list()]) for poly in s2.list()])

left = rndAz - ct
right = w - cs2r - nu

print(left - right)


