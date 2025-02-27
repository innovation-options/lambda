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

    def __str__(self):
        return f'Premium: {self.premium}'


def get_option(rate=.05, strike=20000000, term=1, iterations=4, sigma=.75):
    option = Option(
        strike=strike,
        term=term,
        iterations=iterations,
        sigma=sigma,
        rate=rate
    )
    return {
        'premium': option.premium,
        'term': option.term,
        'strike': option.strike,
        'iterations': option.iterations,
        'sigma': option.sigma,
        'rate': option.rate,
        'value_tree': option.value_tree
    }


def option_handler(event, context):
    rate = event['rate']
    strike = event['strike']
    term = event['term']
    iterations = event['iterations']
    sigma = event['sigma']
    option = get_option(
        rate=rate,
        strike=strike,
        term=term,
        iterations=iterations,
        sigma=sigma,
    )
    res = {
        "statusCode": 200,
        "body": option,
    }
    return res

def compound_option(rate=.05, strike=20000000, term=1, iterations=4, sigma=.75, tranches=4):
    tranch = 0
    compound = {}
    while tranch > -tranches:
        tranch_number = tranch + tranches
        option = get_option(
            rate=rate,
            strike=strike,
            term=term,
            iterations=iterations,
            sigma=sigma,
        )
        compound[f'tranch-{tranch_number}'] = option
        strike = option['premium']
        term = term / 2
        tranch -= 1
    return compound


def compound_handler(event, context):
    rate = event['rate']
    strike = event['strike']
    term = event['term']
    iterations = event['iterations']
    sigma = event['sigma']
    tranches = event['tranches']
    compound = compound_option(
        rate=rate,
        strike=strike,
        term=term,
        iterations=iterations,
        sigma=sigma,
        tranches=tranches,
    )
    res = {
        "statusCode": 200,
        "body": compound,
    }
    return res
