---
name: update-public-surfaces
description: Use when Will has a new paper/preprint, a preprint got published in a journal, or on the weekly sync — updates every public surface (website publications bib, web CV, master LaTeX CV, news feed, preview image) so they stay consistent. Also invoked with an arXiv ID or DOI ("add 2610.12616 to my site").
---

# Update Will Thompson's public surfaces

Repo: `WillHWThompson/WillHWThompson.github.io` (al-folio Jekyll site, branch `main`,
locally at `$CSDS/website.github.io`). Pushing to `main` deploys via
`.github/workflows/deploy.yml`, which also compiles the master CV PDF.

## Surfaces (all must agree)

| Surface | File | Notes |
|---|---|---|
| Publications page + homepage "selected" | `_bibliography/papers3.bib` | Will's own papers. Biblatex style (`date`, `journaltitle`, `eprint`). Newest first. |
| CCM collaboration papers | `_bibliography/ccm.bib` | Coherent CAPTAIN-Mills papers (first author A. A. Aguilar-Arevalo). |
| Web CV | `_data/cv.yml` | Sections: `Publications (Preprints)`, `Publications (Peer-Reviewed Journals)`, `Coherent CAPTAIN-Mills Collaboration`, `Conference Proceedings`. Newest first. |
| **Master CV (PDF)** | `_cv/cv_master.tex` | The "everything" CV. Same lists as cv.yml. Will's name wrapped in `\me{...}`. |
| News feed | `_news/<slug>.md` | One short inline item per new preprint / publication. |
| Preview image | `assets/img/publication_preview/arxiv-<id>.png` | 470x180 PNG, referenced by `preview = {...}` in the bib. |

Not automatable — list them in the final report as manual TODOs: **ORCID** (add the work at
orcid.org), **joint-lab site** (the lab has its own weekly bot; don't touch it),
**Google Scholar** (usually auto-picks up arXiv; confirm the entry merged).

## 1. Find candidates

Existing IDs: collect every `eprint`, `doi`, and arXiv URL in both bib files.

**New arXiv papers.** Query the arXiv API (no key needed; be polite, 3s between calls):
```
https://export.arxiv.org/api/query?search_query=au:Thompson_W&sortBy=submittedDate&sortOrder=descending&max_results=100
```
"W. Thompson" is a very common name. Accept a paper as Will's ONLY if one of:
- an author is "William Thompson", "Will Thompson", or "W. H. W. Thompson"/"William H. W. Thompson", AND at least one co-author already appears in `papers3.bib`/`ccm.bib` or is a known collaborator (Young, Hébert-Dufresne, Landry, Danforth, Dodds, Chodrow, Brooks, Chawla, Boudreau, St-Onge, Allard, Lokhov, DeDeo, Koehler, Tian, Huynh, Wiedemann, Modisette, Tang, Kaiser); or
- it is a Coherent CAPTAIN-Mills collaboration paper (title/abstract mentions CCM / Coherent CAPTAIN-Mills and the author list includes Aguilar-Arevalo and a Thompson).

Anything else that looks plausible → do NOT add; list it under "needs Will's confirmation" in the report.

**Preprint → published.** For every bib entry that is `@online` / has an `eprint` but no `doi`:
query Crossref `https://api.crossref.org/works?query.bibliographic=<title>&rows=3` and accept a match only if the title matches (case/punctuation-insensitive) and the first author's family name matches. Also check the arXiv API entry's `<arxiv:doi>` / `<arxiv:journal_ref>` fields.

If an arXiv ID / DOI was given explicitly (manual invocation), skip discovery and use it.

## 2. Apply updates

Match the formatting of neighbouring entries exactly — copy an existing entry as a template.

**New preprint**
- `papers3.bib` (or `ccm.bib`): `@online{<lastname><FirstTitleWords><year>, title, author (full names, "Last, First" joined by " and "), date = {YYYY}, eprint, eprinttype = {arXiv}, url = {https://arxiv.org/abs/<id>}, preview = {arxiv-<id>.png}, bibtex_show = {true}, html = {https://arxiv.org/pdf/<id>}}` at the top. Add `selected = {true}` if Will is first author.
- `_data/cv.yml`: prepend to `Publications (Preprints)` (or the CCM section): `- A. B. Last, W. H. W. Thompson, ... "[Title](https://arxiv.org/abs/<id>)." arXiv:<id> (YYYY).`
- `_cv/cv_master.tex`: prepend under the matching `\resumeSubheading{...}` list: `\resumeItem{A.B. Last, \me{W.H.W. Thompson}, ..., ``Title.'' arXiv:<id>, YYYY. \href{https://arxiv.org/abs/<id>}{arxiv.org/abs/<id>}}`. Use the name form the paper uses (`W. Thompson` vs `W.H.W. Thompson`). Escape LaTeX specials (`& % $ # _ { } ~ ^`), accents as `H\'ebert`, quotes as ``` ``...'' ```.
- `_news/<short-slug>-<year>.md`:
  ```
  ---
  layout: post
  date: <arXiv submission date> 10:00:00-0400
  inline: true
  ---

  New preprint: [“Title”](https://arxiv.org/abs/<id>), with <co-authors' first+last names>. <One sentence on the main result, from the abstract.>
  ```
  For CCM papers use "The CCM collaboration paper ...".
- Preview image: download the PDF, render Figure 1 (or the main schematic figure — not the title block) with PyMuPDF (`pip install pymupdf pillow`), crop to the figure, fit onto a white 470x180 canvas, save as `assets/img/publication_preview/arxiv-<id>.png`. Look at the result before committing. If no usable figure, crop the title + abstract top region instead.

**Preprint now published**
- bib: `@online` → `@article`; add `date = {YYYY-MM-DD}`, `journaltitle`, `shortjournal`, `volume`, `number`, `pages`, `doi`, `url = {https://doi.org/<doi>}`, `html = {https://doi.org/<doi>}`; keep `eprint`, `preview`, `bibtex_show`, `selected`. Use the journal's title casing.
- cv.yml: move the line from Preprints to `Publications (Peer-Reviewed Journals)` (CCM stays in its section) as `"[Title](https://doi.org/<doi>)." *Journal Abbrev.* **vol** (issue) pages (year).`
- cv_master.tex: same move: `..., ``Title.'' J. Abbrev., vol. V, no. N, PAGES, YYYY. \href{https://doi.org/<doi>}{doi:<doi>}`.
- News item: `[“Title”](https://doi.org/<doi>) is now published in *Journal* **vol**, pages.`

**Author-list conventions (apply to every new entry)**
- Will's name is **bold** everywhere. Web: automatic via `site.scholar.first_name`/`last_name` in `_config.yml` + `_layouts/bib.html` — if he publishes under a new name form, add it to `first_name`. cv.yml: wrap as `**W. H. W. Thompson**`. tex: `\me{...}`.
- Equal contribution: check the PDF's author block for ∗/† footnotes ("contributed equally"). Mark each equal author: bib `Last*, First` (the `*` on the last name; the layout strips it for matching), cv.yml `Name\*`, tex `Name$^*$`. The "Equal contribution" legends already exist (publications.md, cv.yml after the preprints, tex Preprints heading) — if an equal-contribution paper lands in another section, add a legend there too.
- cv.yml YAML gotchas: any list item that starts with `*` or contains `: ` (e.g. a title with a colon) MUST be wrapped in single quotes (double any `'` inside), otherwise YAML reads it as an alias or a key/value map and the CV renders `=>`. After editing, check with `ruby -ryaml -e 'YAML.load_file("_data/cv.yml").each{|x| (x["contents"]||[]).each{|c| abort c.inspect unless x["type"]!="list" || c.is_a?(String)}}'`.

Never delete or rewrite other content. Keep cv.yml and cv_master.tex lists in the same order with the same papers — after editing, diff the two lists (titles) and fix any drift you find, in either direction.

## 3. Verify

- Bib parses: `python3 -c` with `bibtexparser` or at least check braces balance per entry.
- YAML parses: `python3 -c "import yaml;yaml.safe_load(open('_data/cv.yml'))"`.
- LaTeX compiles — a broken `.tex` breaks the whole deploy. Run `latexmk -pdf -interaction=nonstopmode -halt-on-error -outdir=/tmp/cv _cv/cv_master.tex`. If no TeX is installed, try `apt-get install -y --no-install-recommends texlive-latex-recommended texlive-latex-extra texlive-fonts-recommended latexmk`; if that's impossible, re-read every line you added for unescaped specials.
- Don't run `bundle exec jekyll build` — it takes 15+ minutes locally; CI builds it.

## 4. Land it

- Nothing found → change nothing, no commit.
- Otherwise one commit on `main`: `Add <short title> to publications/CV` (or `Sync publications: ...` for several), then `git push origin main`.
- If the push is rejected (non-fast-forward), `git pull --rebase origin main` and push again; never force-push.

## 5. Report

End with a short report: what was added/moved (titles + IDs), the commit hash, candidates that need Will's confirmation, and the manual TODOs (ORCID, Scholar check). When running as the weekly cloud routine and something changed or needs confirmation, email this report to Will with the Gmail connector (subject `[site sync] ...`); if nothing changed and nothing is pending, don't email.
