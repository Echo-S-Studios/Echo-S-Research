# Course-booklet authoring brief (shared by all module agents)

You are authoring **one self-contained HTML lesson-booklet** for the Echo S Research
training site. It must **mirror exactly** the design of the existing booklet at
`papers/2026-07-exact-arithmetic/exact-arithmetic-thirteen-lessons.html` — a horizontal
deck of full-page lessons in a light scholarly theme. Read that booklet before you start;
it is both your design reference and, for the lessons it overlaps, a source of
already-verified content you may adapt.

## Use the shared skeleton — do not reinvent the chrome

1. Read `site/courses/_skeleton.html`. It already contains the complete `<head>`/CSS, the
   fixed header rail + menu, `<main id="deck">`, and the navigation `<script>`. **Copy it
   verbatim** as your file's shell.
2. Replace the three placeholders: `@TITLE@` (the `<title>`), `@BRAND@` (short bold name in
   the header, e.g. "The Two Walls"), `@SUBTITLE@` (the series line, e.g. "Mahler, Salem &
   Lehmer").
3. Replace the `<!-- … @LESSONS@ … -->` marker with your lesson pages. **Change nothing
   else** — not the `<head>`, CSS, header, or `<script>`. The header ticks, the menu, the
   prev/next buttons and keyboard nav are **auto-built from each page's `data-title`**, so you
   only add pages.

## Markup — reuse the booklet's classes, invent no new CSS

Read the exact-arithmetic booklet for the precise markup of every component:
`.page` / `.inner`, `.eyebrow`, `h1` / `h2` / `.dek`, `h3` (§ section heads), `.eq`
equation strips (each with a short `.cid` check-id in the corner and `<b>` on the key term),
`.mat` matrices (`.mat.c2` / `.mat.c3`), `.tag` epistemic chips
(`forced` / `computed` / `declared` / `plausible` / `open` / `estab`), the `.thread` + `.chips`
"inherited-from" block, the `.handoff` "hand-on-to-next" block, the `.run` verification
badge, `details.cp` checkpoint question cards (`.qn` / `.qtext` / `.reveal` / `.ans`),
`.nextup`, and `.endnote`. Use these classes; add no `<style>` of your own.

## Page order — fill `<main id="deck">` with:

1. one `<section class="page cover" data-title="Cover">` — cover: an `.eyebrow`, `h1` = the
   module title, a one-sentence `.dek` thesis, a short "what you'll learn / how to read this",
   and a nav line: `<a href="../learn.html">← The course</a> · <a href="../index.html">The archive</a>`.
2. one `<section class="page" data-title="N · Short Title">` **per lesson**. Each lesson:
   `.eyebrow` = "Lesson N · <series>", `<h2 tabindex="-1" data-pagehead>` title, a `.dek`,
   then `h3` § sections of taught prose with the real formulae inline, `.eq` strips for the
   load-bearing equations, `.tag` chips on load-bearing claims, at least one `details.cp`
   checkpoint question, and a `.run`/verify note where the claim is machine-checked. End each
   lesson with a `.nextup` pointing to the next.
3. one final `<section class="page" data-title="Ledger">` — a compact claim-ledger table
   (what is **forced / declared / open** in this module) and the module's papers with links.

## Content rules — accuracy is the whole point

- **Teach only from this module's papers.** Read each paper's LaTeX source (paths in your
  task) and restate its definitions, theorems, and the real numbers **faithfully**. Do **not**
  invent theorems, constants, or claims. Every load-bearing statement must be traceable to a
  paper, and must carry that paper's own epistemic tag. When unsure, cite/quote rather than
  assert.
- It is a *course*, not a transcript: explain intuition, connect the papers, build simple →
  hard. But every stated result is the paper's actual result.
- **Link each paper** by its served PDF `../papers/<shortname>.pdf` (LaTeX source on GitHub:
  `https://github.com/Echo-S-Studios/Echo-S-Research/tree/main/papers/<shortname>` if useful).
  The exact-arithmetic booklet is at `../papers/2026-07-exact-arithmetic.html`. Point
  "verify it" notes at the repo's `tests/` or `code/` trees on GitHub.
- **Self-contained**: no external scripts, fonts, stylesheets, or images. **Unicode math only**
  (φ, √5, ℤ/4ℤ, ⊕, ⊗, ψ², ∏, ∫, ρ, κ, λ, …) — **no MathJax**.
- Voice: match the booklet — precise, unshowy, the "read the source → run an exact check →
  teach only what passed" ethos. The epistemic tags are the mathematics, not decoration.

## Output

Write the finished file to `site/courses/<slug>.html` (slug in your task). Prefer **fewer,
genuinely-teaching lessons** over many shallow ones. When done, report: the file path, the
lesson titles in order, and one sentence on what a reader can do after finishing.
