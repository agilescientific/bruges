# -*- coding: utf-8 -*-
"""
Make random text for Bruges.

"""


def get_bruges(p, n):

    if p < 0:
        p = 0
    if p > 1:
        p = 1

    if n < 0:
        n = 0
    if n > 4:
        n = 4

    b = [
        "box",
        "bevy",
        "buttload",
        "billions",
        "bunch",
        "ball",
        "bucket",
        "bowlful",
        "barrowload",
        ]

    r = [
        "ridiculously",
        "rather",
        "reasonably",
        "regrettably",
        "resolutely",
        "riotously",
        "reassuringly",
        "rarely",
        "risible or"
        ]

    u = [
        "unlikely",
        "untested",
        "unusable",
        "undocumented",
        "uninteresting",
        "unexpected",
        "untrue",
        "unbelievable",
        "useless"
        ]

    s = [
        ", supposedly",
        ", sorta",
        " or something",
        " and shit",
        ", surprisingly",
        " and so on",
        " and such",
        ", see?",
        ", stupid"
        ]

    B = "bag"
    R = "really"
    U = "useful"
    S = " and stuff"

    import random
    if random.random() < p:
        words = random.sample([0, 1, 2, 3], n)
        if 0 in words:
            B = random.choice(b)
        if 1 in words:
            R = random.choice(r)
        if 2 in words:
            U = random.choice(u)
        if 3 in words:
            S = random.choice(s)

    return "{} of {} {} geophysical equations{}".format(B, R, U, S)
