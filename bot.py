import yaml, feedparser, requests, datetime as dt
from filters import score_paper, match_category
from mailer import send_email
from parsers import try_extract_conclusion


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
    base = "https://export.arxiv.org/api/query"
    params = {
        "search_query": query,
        "sortBy": "submittedDate",
        "sortOrder": "descending",
        "max_results": str(max_results),
        "start": "0"
    }

    # direct single-page request (no arxiv.Client pagination)
    url = base + "?" + "&".join(f"{k}={v}" for k, v in params.items())
    resp = requests.get(url, timeout=30)
    feed = feedparser.parse(resp.text)

    results = []
    for entry in feed.entries:
        results.append({
            "title": entry.title,
            "summary": entry.summary,
            "published": dt.datetime.fromisoformat(entry.published.replace("Z", "+00:00")),
            "link": entry.link,
            "authors": [a.name for a in entry.authors],
        })

    days_back = cfg["arxiv"].get("days_back", 1)
    since = dt.datetime.now(dt.timezone.utc) - dt.timedelta(days=days_back)
    results = [r for r in results if r["published"] >= since]

    print(f"[astro-ph bot] fetched {len(results)} papers (max {max_results})")
    return results


def curate(cfg, papers):
    keep = []
    for r in papers:
        if not match_category(r.primary_category, cfg["arxiv"]["categories"]):
            continue
        score, details = score_paper(
            title=r.title or "",
            abstract=r.summary or "",
            authors=[a.name for a in r.authors],
            prefs=cfg["preferences"],
        )
        if score >= cfg["preferences"].get("min_score", 1.0):
            keep.append((score, details, r))
    keep.sort(key=lambda x: (-x[0],))
    return keep


def format_authors(authors, max_authors=5):
    """return 'A. Author et al.' if long list."""
    names = [a.name for a in authors]
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
    for score, details, r in curated:
        url = r.entry_id
        title = r.title.strip().replace("\n", " ")
        abstract = (r.summary or "").strip()
        authors_line = format_authors(r.authors)
        conclusion = ""
        if cfg["output"].get("include_conclusion", True):
            conclusion = try_extract_conclusion(r.pdf_url)

        # TEXT block
        lines_txt.append(f"{title}\n{url}\n")
        if authors_line:
            lines_txt.append(f"Authors: {authors_line}\n")
        lines_txt.append(abstract + "\n")
        if conclusion:
            lines_txt.append("Conclusion:\n" + conclusion + "\n")
        lines_txt.append("-" * 60)

        # HTML block
        lines_html.append("<li>")
        lines_html.append(f'<p><b><a href="{url}">{title}</a></b></p>')
        if authors_line:
            lines_html.append(f"<p><i>{authors_line}</i></p>")
        lines_html.append(f"<p>{abstract}</p>")
        if conclusion:
            lines_html.append(f"<p><i>Conclusion:</i> {conclusion}</p>")
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
