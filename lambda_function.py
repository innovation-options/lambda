import json
from collections import OrderedDict
from math import exp, log, pow, sqrt
from operator import mul, sub
from operator import truediv as div


class Option:
    def __init__(self, rate, strike, term, iterations, sigma):
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
        v = OrderedDict(
            sorted(
                v.items(),
            )
        )
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
        p = OrderedDict(
            sorted(
                p.items(),
                reverse=True,
            )
        )
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
            out = OrderedDict(
                sorted(
                    out.items(),
                    reverse=True,
                )
            )
            v[j] = out
            j -= 1
        v = OrderedDict(
            sorted(
                v.items(),
            )
        )
        return v

    @property
    def premium(self):
        return self.value_tree[0][0]


    def as_dict(self):
        return {
            'rate': self.rate,
            'strike': self.strike,
            'term': self.term,
            'iterations': self.iterations,
            'sigma': self.sigma,
            'premium': self.premium,
        }

    def as_json(self):
        return json.dumps(self.as_dict())

def lambda_handler(event, context):
    inputs = json.loads(event['body'])
    # inputs = {
    #     'rate': .05,
    #     'strike': 20,
    #     'term': 1,
    #     'iterations': 12,
    #     'sigma': .75,
    # }
    option = Option(
        float(inputs['rate']),
        float(inputs['strike']),
        float(inputs['term']),
        float(inputs['iterations']),
        float(inputs['sigma']),
    )
    res = {
        "statusCode": 200,
        "headers": {
            "Content-Type": "application/json"
        },
        "body": option.as_json()
    }
    return res
