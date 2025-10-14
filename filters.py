import fnmatch, re, unicodedata


def _normalize(s):
    # lowercase, unicode normalize, strip accents/diacritics, replace ß with ss
    # invent a more robust approach if needed for other languages
    s = s.lower()
    s = unicodedata.normalize("NFKD", s)
    s = s.replace("ß", "ss")
    # remove combining accents
    s = "".join(c for c in s if not unicodedata.combining(c))
    # collapse whitespace
    s = re.sub(r"\s+", " ", s)
    return s


def match_category(primary_cat, patterns):
    if not patterns: return True
    return any(fnmatch.fnmatch(primary_cat or "", pat) for pat in patterns)


def score_paper(title, abstract, authors, prefs):
    t = _normalize(title); a = _normalize(abstract)
    joined = f"{t} {a}"
    any_kw = prefs.get("any_keywords", [])
    all_kw = prefs.get("all_keywords", [])
    ex_kw  = prefs.get("exclude_keywords", [])
    auth_wl = [s.lower() for s in prefs.get("authors", [])]

    score = 0.0
    details = {}

    # keyword hits
    any_hits = sum(1 for k in any_kw if k.lower() in joined)
    score += any_hits
    details["any_hits"] = any_hits

    if all_kw:
        all_ok = all(k.lower() in joined for k in all_kw)
        if all_ok:
            score += 2.0
        details["all_ok"] = all_ok

    # authors
    al = [x.lower() for x in authors]
    auth_hits = sum(1 for wl in auth_wl if any(wl in au for au in al))
    score += auth_hits
    details["auth_hits"] = auth_hits

    # excludes
    ex_hits = sum(1 for k in ex_kw if k.lower() in joined)
    if ex_hits:
        score -= 2.0 * ex_hits
    details["ex_hits"] = ex_hits

    return score, details
