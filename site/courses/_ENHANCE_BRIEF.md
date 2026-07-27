# Course-enhancement brief (equal-depth pass)

You are **enhancing one existing course booklet** so a motivated non-expert can actually
follow the mathematics — and so that **every lesson in the booklet has equal depth**. You
are improving what is already there, not starting over. Read `_AGENT_BRIEF.md` first (the
design/fidelity contract still applies in full), then this.

## Hard constraints (unchanged)

- **Do not touch the chrome.** The `<head>`/CSS and the navigation `<script>` must stay
  **byte-identical** to `_skeleton.html`. Edit only the lesson `<section class="page">`
  bodies inside `<main id="deck">`. Leave the already-filled `@TITLE@`/`@BRAND@`/`@SUBTITLE@`,
  the cover, and the ledger structure intact (you may enrich the cover/ledger lightly, but
  the deep work is the lessons).
- **Reuse the existing CSS classes only** (`.eq`/`.cid`, `.mat`, `.tag`, `.thread`/`.chips`,
  `.handoff`, `.run`, `details.cp` with `.qn`/`.qtext`/`.reveal`/`.ans`, `.nextup`, `.gapfig`,
  …). Add no `<style>`.
- **Self-contained**, Unicode math only (matching the booklet's existing convention — HTML
  numeric entities like `&#966;`=φ and/or literal Unicode), **no MathJax**, no external assets.
- **Faithful.** Teach only from this module's papers (read their `.tex`). Enhancements are
  *pedagogical* — clearer explanations, intuition, worked examples that follow from the
  paper — **never new theorems, numbers, or claims**. Every load-bearing statement stays
  traceable to a paper and keeps its epistemic tag. Do not weaken a tag to make a story
  cleaner.

## What "improve understanding" means — apply to EVERY lesson

1. **Motivate before formalizing.** Begin each lesson (and ideally each `§`) with *why* this
   object/result matters and what question it answers, in plain words, before the notation.
2. **Build intuition around every load-bearing formula.** After an `.eq` strip, add a
   sentence or two saying what it *means* — the mechanism, a picture in prose, an analogy —
   and define notation the first time it appears. A reader should never meet a symbol cold.
3. **Work at least one concrete instance per lesson, in full.** Take the paper's own smallest
   example (a specific matrix, a specific measure, a specific polynomial/number) and run the
   machinery once, step by step, so the reader sees it actually compute. Put the worked
   instance behind a `details.cp` checkpoint where natural: pose the question, reveal the
   steps.
4. **Make the checks concrete.** State what the cited exact-arithmetic check actually
   computes and why passing it certifies the claim — not just a badge count.
5. **Thread the arc.** Each lesson names what it inherits (`.thread`/`.chips`) and what it
   hands forward (`.handoff`/`.nextup`), so the whole booklet reads as one line of thought.

## Equal depth — the requirement, made concrete

Every lesson must end at **comparable, substantial depth**. Target for each lesson:

- **4–5 `§` sections**, each doing real teaching (not a single sentence);
- **≥ 1 fully worked concrete instance**;
- **≥ 3 checkpoint cards** (`details.cp`), at least one of which is a genuine worked exercise;
- a **motivation opening** and a **`.nextup` close**;
- **~1500–1800 words** of taught prose.

**Procedure:** first measure your booklet — count the `§`s, checkpoints, and rough word count
of each lesson (your task lists the current per-lesson word counts). Then **level the thin
lessons up to the richest**, and enrich the rest to the target. **Never thin a rich lesson to
match a poor one.** Deepen *genuinely* — more motivation, more intuition, more worked steps,
more "why" — and never pad with restatement or filler. If a lesson is short because its topic
is genuinely small, give it a fuller worked example and a sharper intuition rather than empty
words.

## Output

Edit `site/courses/<slug>.html` in place. When done, confirm (and report):
- the `<head>`/CSS and `<script>` are still byte-identical to `_skeleton.html` (diff them);
- the file is self-contained and has no leftover placeholders;
- **per-lesson depth, before → after** (§ counts and rough word counts), showing the lessons
  are now within a tight band;
- one sentence on what a reader can now do that they couldn't before.
