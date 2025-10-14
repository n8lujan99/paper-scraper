from filters import _normalize, score_paper

prefs = {
    "authors": ["Gußmann", "Anantua"],
    "any_keywords": ["black hole"],
}

authors_variants = [
    "T. Gußmann",
    "T. Gussmann",
    "GUßMANN, T.",
]

print("Normalization check:")
for a in authors_variants:
    print(f"{a!r} → {_normalize(a)}")

title = "Black Hole Imaging with the EHT"
abstract = "We present new polarized GRMHD simulations."
for a in authors_variants:
    score, details = score_paper(title, abstract, [a], prefs)
    print(f"{a!r}: score={score}, details={details}")
