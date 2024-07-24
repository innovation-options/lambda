import json
from math import exp, log, pow, sqrt
from operator import mul, sub
from operator import truediv as div


class Option:
    def __init__(self, rate=.05, strike=20, term=1, iterations=4, sigma=.75):
        self.rate = rate
        self.strike = strike
        self.spot = strike
        self.term = term
        self.iterations = iterations
        self.sigma = sigma

    @property
    def delta_t(self):
        return self.term / self.iterations

    @property
    def up_factor(self):
        return exp(
            mul(
                self.sigma,
                sqrt(
                    mul(
                        2,
                        self.delta_t
                    )
                )
            )
        )

    @property
    def down_factor(self):
        return div(
            1,
            self.up_factor
        )

    @property
    def flat_factor(self):
        return 1

    @property
    def up_prob(self):
        return pow(
            div(
                sub(
                    exp(
                        div(
                            mul(
                                self.rate,
                                self.delta_t,
                            ),
                            2
                        )
                    ),
                    exp(
                        mul(
                            -self.sigma,
                            sqrt(
                                div(
                                    self.delta_t,
                                    2
                                )
                            )
                        )
                    )
                ),
                sub(
                    exp(
                        mul(
                            self.sigma,
                            sqrt(
                                div(
                                    self.delta_t,
                                    2
                                )
                            )
                        )
                    ),
                    exp(
                        mul(
                            -self.sigma,
                            sqrt(
                                div(
                                    self.delta_t,
                                    2
                                )
                            )
                        )
                    )
                )
            ),
            2
        )

    @property
    def down_prob(self):
        return pow(
            div(
                sub(
                    exp(
                        mul(
                            self.sigma,
                            sqrt(
                                div(
                                    self.delta_t,
                                    2
                                )
                            )
                        )
                    ),
                    exp(
                        div(
                            mul(
                                self.rate,
                                self.delta_t,
                            ),
                            2
                        )
                    ),
                ),
                sub(
                    exp(
                        mul(
                            self.sigma,
                            sqrt(
                                div(
                                    self.delta_t,
                                    2
                                )
                            )
                        )
                    ),
                    exp(
                        mul(
                            -self.sigma,
                            sqrt(
                                div(
                                    self.delta_t,
                                    2
                                )
                            )
                        )
                    )
                )
            ),
            2
        )

    @property
    def flat_prob(self):
        return 1 - (self.up_prob + self.down_prob)

    @property
    def spot_tree(self):
        """ Step One: Find Spot Value at Terminus."""
        v = {}
        i = -self.iterations
        while i <= self.iterations:
            v[i] = mul(
                self.strike,
                pow(
                    self.up_factor,
                    i
                )
            )
            i += 1
        return v

    @property
    def term_value(self):
        """ Step Two: Find Option Value at Terminus."""
        p = self.spot_tree
        i = -self.iterations
        while i <= self.iterations:
            p[i] = max(
                p[i] - self.strike,
                0.0
            )
            i += 1
        return p

    @property
    def value_tree(self):
        """Step Three: Calculate option value at earlier nodes:"""
        p = self.term_value
        v = {}
        v[self.iterations] = p
        j = self.iterations - 1
        while j >= 0:
            i = -j
            out = {}
            while i <= j:
                out[i] = mul(
                    exp(
                        mul(
                            -float(self.rate),
                            self.delta_t
                        )
                    ),
                    mul(
                        self.up_prob,
                        v[j + 1][i + 1]
                    ) +
                    mul(
                        self.down_prob,
                        v[j + 1][i - 1]
                    ) +
                    mul(
                        self.flat_prob,
                        v[j + 1][i]
                    )
                )
                i += 1
            v[j] = out
            j -= 1
        return v

    @property
    def premium(self):
        return self.value_tree[0][0]

def compound_medium(sigma):
    zero = Option(
        strike=20.0,
        term=1.0,
        iterations=4,
        sigma=sigma,
    )
    one = Option(
        strike=zero.premium,
        term=.5,
        iterations=6,
        sigma=sigma,
    )
    two = Option(
        strike=one.premium,
        term=.25,
        iterations=6,
        sigma=sigma,
    )
    three = Option(
        strike=two.premium,
        term=.08,
        iterations=6,
        sigma=sigma,
    )
    return three.premium





def lambda_handler(event, context):
    option = Option(
        event['rate'],
        event['strike'],
        event['term'],
        event['iterations'],
        event['sigma'],
    )
    res = {
        "statusCode": 200,
        "body": option.value_tree,
    }
    return res
