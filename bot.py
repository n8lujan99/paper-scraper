import os, yaml, feedparser, requests, datetime as dt
from filters import score_paper, match_category
from mailer import send_email

print(f"[astro-ph bot] running: {__file__} SHA={os.environ.get('GITHUB_SHA', 'local')}")

def load_config(path="config.yaml"):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def build_search_query(cfg):
    cats = cfg["arxiv"].get("categories", ["astro-ph*"])
    cat_query = " OR ".join([f"cat:{c}" for c in cats])
    return f"({cat_query})"


def fetch_recent(cfg):
    max_results = cfg["arxiv"].get("max_results", 50)
    query = build_search_query(cfg)
    days_back = cfg["arxiv"].get("days_back", 1)
    since = dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=days_back)

    base = "https://export.arxiv.org/api/query"
    params = {
        "search_query": query,
        "sortBy": "submittedDate",
        "sortOrder": "descending",
        "max_results": str(max_results),
        "start": "0",
    }
    url = base + "?" + "&".join(f"{k}={v}" for k, v in params.items())

    try:
        resp = requests.get(url, timeout=20)
        resp.raise_for_status()
    except Exception as e:
        print(f"[ERROR] could not fetch from arXiv: {e}")
        return []

    feed = feedparser.parse(resp.text)
    results = []
    for entry in feed.entries:
        try:
            pub = dt.datetime.fromisoformat(entry.published.replace("Z", "+00:00"))
        except Exception:
            pub = dt.datetime.now(dt.timezone.utc)
        if pub >= since:
            results.append({
                "title": entry.title.strip(),
                "summary": entry.summary.strip(),
                "published": pub,
                "link": entry.link,
                "authors": [a.name for a in entry.authors],
            })

    print(f"[astro-ph bot] fetched {len(results)} papers (max {max_results})")
    return results


def curate(cfg, results):
    curated = []
    prefs = cfg["preferences"]

    for r in results:
        # extract core fields safely (works for dicts or arxiv.Result objects)
        title = r.get("title") if isinstance(r, dict) else getattr(r, "title", "")
        summary = r.get("summary") if isinstance(r, dict) else getattr(r, "summary", "")
        authors = (
            r.get("authors") if isinstance(r, dict)
            else [a.name for a in getattr(r, "authors", [])]
        )
        url = r.get("url") if isinstance(r, dict) else getattr(r, "entry_id", "")
        pdf_url = r.get("pdf_url") if isinstance(r, dict) else getattr(r, "pdf_url", "")
        category = (
            getattr(r, "primary_category", None)
            if not isinstance(r, dict)
            else cfg["arxiv"]["categories"][0]
        )

        # skip if category doesn’t match user’s config
        if not match_category(category, cfg["arxiv"]["categories"]):
            continue

        # unified scoring
        score, details = score_paper(title, summary, authors, prefs)

        if score >= prefs.get("min_score", 1.0):
            curated.append({
                "title": title,
                "summary": summary,
                "authors": authors,
                "url": url,
                "pdf_url": pdf_url,
                "category": category,
                "score": score,
                "details": details,
            })

    print(f"[astro-ph bot] curated {len(curated)} papers (score ≥ {prefs.get('min_score', 1.0)})")
    return curated


def format_authors(authors, max_authors=5):
    """return 'A. Author et al.' if long list."""
    names = [a.name for a in authors] if authors and hasattr(authors[0], "name") else authors
    if not names:
        return ""
    if len(names) > max_authors:
        first = names[0].split(",")[0].strip()
        return f"{first} et al."
    else:
        return ", ".join(names)


def make_email_body(cfg, curated):
    lines_txt = []
    lines_html = ['<html><body><h2>astro-ph digest</h2><ol>']

    for r in curated:
        # safely extract fields (works for dicts or arxiv.Result objects)
        title = (r.get("title") if isinstance(r, dict) else getattr(r, "title", "")) or ""
        abstract = (r.get("summary") if isinstance(r, dict) else getattr(r, "summary", "")) or ""
        url = (r.get("url") if isinstance(r, dict) else getattr(r, "entry_id", "")) or ""
        authors = (r.get("authors") if isinstance(r, dict) else getattr(r, "authors", [])) or []
        authors_line = ", ".join(authors) if isinstance(authors, list) else str(authors)

        # TEXT block
        lines_txt.append(f"{title}\n{url}\n")
        if authors_line:
            lines_txt.append(f"authors: {authors_line}\n")
        lines_txt.append(abstract + "\n")
        lines_txt.append("-" * 60)

        # HTML block
        lines_html.append("<li>")
        lines_html.append(f'<p><b><a href="{url}">{title}</a></b></p>')
        if authors_line:
            lines_html.append(f"<p><i>{authors_line}</i></p>")
        lines_html.append(f"<p>{abstract}</p>")
        lines_html.append("</li>")

    lines_html.append("</ol></body></html>")
    return "\n".join(lines_txt), "\n".join(lines_html)


def main():
    cfg = load_config()
    papers = fetch_recent(cfg)
    curated = curate(cfg, papers)
    if not curated:
        print("no matches today.")
        subject = f'{cfg["output"]["email"]["subject_prefix"]} {dt.date.today()} — 0 papers'
        send_email(cfg, subject, "no matching papers found today.", "<p>no matches today.</p>")
        return

    text_body, html_body = make_email_body(cfg, curated)
    n = len(curated)
    subject = f'{cfg["output"]["email"]["subject_prefix"]} {dt.date.today()} — {n} paper{"s" if n != 1 else ""}'
    send_email(cfg, subject, text_body, html_body)
    print(f"emailed {n} curated papers.")

if __name__ == "__main__":
    main()
