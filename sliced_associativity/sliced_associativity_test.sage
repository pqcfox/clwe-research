import os.path
from collections import namedtuple
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt


LWEConfig = namedtuple('LWEConfig', 'name k l s_max A_max s_prime_max q')


def get_sample(config):
    s = matrix(RR, config.l, 1, [uniform(0, config.s_max) for _ in range(config.l)])
    s_prime = matrix(RR, 1, config.k, [randint(0, config.s_prime_max) for _ in range(config.k)])
    A = matrix(RR, config.k, config.l, [randint(0, config.A_max) for _ in range(config.k * config.l)])
    return s, A, s_prime


def run_test(test_label, config, left_func, right_func, num_tests=10):
    diffs = []
    print(f'Testing {test_label}...')
    for _ in tqdm(range(num_tests)):
        s, A, s_prime = get_sample(config)
        diff = left_func(s, A, s_prime, config.q) - right_func(s, A, s_prime, config.q)
        try:
            diffs.append(float(diff))
        except TypeError:
            diffs.append(float(diff.list()[0]))

    plt.hist(diffs, bins=100)
    plt.savefig(os.path.join(config.name, test_label))
    plt.close()


def slice_mat(a):
    a_copy = matrix(RR, a.nrows(), a.ncols(), a.list())
    slices = []
    while True:
        current = matrix(RR, a.nrows(), a.ncols(), [0 if v < 1 else 1 for v in a_copy.list()])
        a_copy -= current
        if list(set(current.list())) == [0]:
            break
        slices.append(current)

    if len(slices) == 0:
        slices = [a_copy]
    return slices


def run_all_tests(config):
    print(f'Processing {config.name}...')

    try:
        os.mkdir(config.name)
    except FileExistsError:
        pass

    run_test('associativity', config,
             lambda s, A, s_prime, _: (s_prime * A) * s,
             lambda s, A, s_prime, _: s_prime * (A * s))

    run_test('rounded_associativity', config,
             lambda s, A, s_prime, _: (s_prime * A).apply_map(round) * s,
             lambda s, A, s_prime, _: s_prime * (A * s).apply_map(round))

    run_test('floored_associativity', config,
             lambda s, A, s_prime, _: (s_prime * A).apply_map(floor) * s,
             lambda s, A, s_prime, _: s_prime * (A * s).apply_map(floor))

    run_test('mod_associativity', config,
             lambda s, A, s_prime, q: ((s_prime * A) % q * s) % q,
             lambda s, A, s_prime, q: (s_prime * (A * s) % q) % q)

    run_test('rounded_mod_associativity', config,
             lambda s, A, s_prime, q: ((s_prime * A).apply_map(round) % q * s).apply_map(round) % q,
             lambda s, A, s_prime, q: (s_prime * (A * s).apply_map(round) % q).apply_map(round) % q)

    run_test('floored_mod_associativity', config,
             lambda s, A, s_prime, q: ((s_prime * A).apply_map(floor) % q * s).apply_map(floor) % q,
             lambda s, A, s_prime, q: (s_prime * (A * s).apply_map(floor) % q).apply_map(floor) % q)

    run_test('rounded_sliced_tweaked_associativity', config,
             lambda s, A, s_prime, _: sum([s_prime * (sl * s).apply_map(round) for sl in slice_mat(A)]),
             lambda s, A, s_prime, _: sum([sl * s for sl in slice_mat(s_prime * A)]).apply_map(round))

    run_test('sliced_associativity', config,
             lambda s, A, s_prime, _: sum([s_prime * (sl * s) for sl in slice_mat(A)]),
             lambda s, A, s_prime, _: sum([sl * s for sl in slice_mat(s_prime * A)]))

    run_test('rounded_sliced_associativity', config,
             lambda s, A, s_prime, _: sum([s_prime * (sl * s).apply_map(round) for sl in slice_mat(A)]),
             lambda s, A, s_prime, _: sum([(sl * s).apply_map(round) for sl in slice_mat(s_prime * A)]))

    run_test('rounded_sliced_mod_associativity', config,
             lambda s, A, s_prime, q: sum(
                 [s_prime * ((sl * s).apply_map(round) % q) for sl in slice_mat(A)]) % q,
             lambda s, A, s_prime, q: sum(
                 [sl * s for sl in slice_mat((s_prime * A) % q)]).apply_map(round) % q)


run_all_tests(LWEConfig(name='basic_values', k=10, l=10, s_max=10, A_max=10, s_prime_max=10, q=10))
run_all_tests(LWEConfig(name='no_mult_wrap', k=10, l=10, s_max=2, A_max=2**13, s_prime_max=2, q=2**13))
