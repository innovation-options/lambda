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

def lambda_handler(event, context):
    body = json.loads(event['body'])
    # option = Option(
    #     float(body['rate']),
    #     float(body['strike']),
    #     float(body['term']),
    #     float(body['iterations']),
    #     float(body['sigma']),
    # )
    premium = str(6.1)
    res = {
        "statusCode": 200,
        "headers": {
            "Content-Type": "*/*"
        },
        "body": premium
    }
    return res
# def lambda_handler(event, context):
#     body = json.loads(event['body'])
#     res = {
#         "statusCode": 200,
#         "headers": {
#             "Content-Type": "*/*"
#         },
#         "body": f'Hello, {body['greeter']}!'
#     }
#     return res
