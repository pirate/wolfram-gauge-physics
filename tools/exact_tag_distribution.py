"""Exact Bernoulli-ensemble tracer CDF from free occupancy motion and conserved rank."""
import math
from functools import lru_cache


@lru_cache(maxsize=8192)
def binomial_pmf(count, density):
    if count == 0: return (1.,)
    if density == 0: return (1.,)+(0.,)*count
    if density == 1: return (0.,)*count+(1.,)
    logs = [math.lgamma(count+1)-math.lgamma(k+1)-math.lgamma(count-k+1)
            +k*math.log(density)+(count-k)*math.log1p(-density) for k in range(count+1)]
    largest = max(logs)
    values = [math.exp(value-largest) for value in logs]
    total = sum(values)
    return tuple(value/total for value in values)


def phase_cdf(layers, position, density, phase):
    # Production checkpoints use complete pairs of layers, so each occupancy bit
    # has moved exactly +/-layers in the interior. The tag's color retains its rank.
    if layers % 2:
        raise ValueError('closed-form routing requires an even number of layers')
    if position < -layers: return 0.
    if position >= layers: return 1.
    plus = minus = fixed = 0
    for original in range(-2*layers-1, 2*layers+2):
        velocity = 1 if (original+phase) % 2 == 0 else -1
        coefficient = int(original+velocity*layers <= position)-int(original <= 0)
        if original == 0: fixed = coefficient
        elif coefficient == 1: plus += 1
        elif coefficient == -1: minus += 1
    positive, negative = binomial_pmf(plus, density), binomial_pmf(minus, density)
    cumulative, value = [], 0.
    for probability in negative:
        value += probability
        cumulative.append(value)
    # P(U - V + fixed >= 0), with U,V independent binomials.
    return sum(p*(0 if k+fixed < 0 else 1 if k+fixed >= minus else cumulative[k+fixed])
               for k, p in enumerate(positive))


def tag_distribution(layers, density):
    cdf = [(phase_cdf(layers, x, density, 0)+phase_cdf(layers, x, density, 1))/2
           for x in range(-layers, layers+1)]
    probabilities, previous = [], 0.
    for value in cdf:
        probabilities.append(max(0., value-previous))
        previous = value
    mean = sum((i-layers)*p for i, p in enumerate(probabilities))
    variance = sum(((i-layers)-mean)**2*p for i, p in enumerate(probabilities))
    return {'cdf': cdf, 'pmf': probabilities, 'mean': mean, 'variance': variance}
