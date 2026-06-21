# lessons.md — your prep log

One file for the whole prep. Keep two kinds of entry in **separate subsections each week** — don't mix them in one paragraph.

| Subsection | What goes here | How much detail |
|------------|----------------|-----------------|
| **From the materials** | Notes while watching or reading; answers to each week's reflection exercises | Usually one sentence per video chunk or paper section; reflection exercises can be a short paragraph each |
| **Surprises** | Moments an LLM or agent surprised you — good or bad — in chat or in the IDE | Concrete: tool/model, what you asked, what came back, optional takeaway |

Commit and push weekly. By week 4 this file is one of the most useful artifacts you bring to Brno. (`reflection.md` in week 4 is separate — one final paragraph for assessment.)

---

## From the materials — what to write

While watching or reading, stop every ~20 minutes (or after each major section) and add a line answering:

- *Video:* **What's the one thing I'd want to test from what I just heard?**
- *Paper:* **What claim would I most want to verify on my own data?**

Each week may also assign a **reflection exercise** (structured thinking using the week's mental model). Put those answers here too — they are not required to be personal chat logs.

---

## Surprises — what to write

Add an entry whenever an LLM or agent catches you off guard. Include enough detail that you (or a classmate) could understand the moment months later.

- **When** — approximate date
- **Tool / model** — e.g. ChatGPT (free), Claude, Antigravity agent, Cursor, …
- **What you asked** — paste or paraphrase the prompt; name any file or data involved
- **What happened** — the surprising part
- **Takeaway** (optional) — one line on what you'd do differently

**Bad (too vague):** *"ChatGPT hallucinated something."*

**Good:**

> **2026-05-26 · ChatGPT (free, no browsing)** — Asked: *"What is the Ensembl ID for human BRCA1?"* Answered confidently with `ENSG00000012048` — correct — then cited a made-up paper (*Smith et al., Nature 2019*) and a DOI that 404s. **Takeaway:** right gene, invented provenance; never trust citations without checking.

> **2026-06-03 · Antigravity agent** — Asked it to filter a BED file to chr21. Code ran, printed 1,842 lines, looked plausible. Checked: 0-based coordinates on a file the header said was 1-based. **Takeaway:** spot-check coordinate conventions before trusting counts.

---

## Your entries

(Add below. Newest at the bottom is fine — stay consistent.)

### Week 1

#### From the materials

- **Task:** Counting characters or letters (e.g., "How many 'r's in strawberry?")

- **Why it's hard:** The model doesn't see characters; it sees tokens (atoms of text). Because it isn't "looking" at the character level, it lacks a direct visual sense of the word's structure, and its internal "mental arithmetic" is prone to failing at simple counting tasks.

- **Domain Check:** Always treat string manipulation or character counts as untrustworthy. Verify by asking the model to write and execute a Python script to count them.

#### Surprises

Model/Agent: Gemini 1.5 Pro

- **Task:** Sorting a list of CSV file names by date embedded in the filename (e.g., data_20260105.csv, data_20260104.csv).

- **The Vanilla Gap:** The vanilla model got the order wrong because it treated the numbers as strings in a way that ignored the actual date logic, or hallucinated a sort order based on file length rather than the date sequence.

- **The Code Execution Gap:** Using code execution, the model wrote a simple Python script using datetime objects to parse the strings. It was 100% accurate because it offloaded the logic to a deterministic interpreter rather than using its "mental" statistical patterns.

Scientific Validity Check:

- **Assumption made:** The model assumed the format YYYYMMDD was consistent across all files.

- **Surprise:** The model correctly identified that a file missing a date suffix would break the script and added an if statement to handle the error—a level of "defensive programming" I didn't explicitly prompt for, likely because it saw similar patterns in its training data for robust code.

### Week 2

#### From the materials

- **Software 3.0:** Engineering shifts from manual coding to designing robust evaluations and strict data invariants for AI agents.
- **ReAct Protocol:** Interleaving reasoning and tool use enables agent self-correction and minimizes hallucinations.
- **Architecture:** Predictable, orchestrated workflows (e.g., prompt chaining) outperform highly autonomous agents. 
- **Validation Scaling:** Programmatic invariant checking is essential to catch silent agent logic failures (e.g., missing reverse operations, version mismatches).

#### Surprises

Model/Agent: Gemini 3.5 Flash (High)

- **Silent Logic Bugs:** Agents frequently introduce off-by-one errors across different indexing standards (0-based vs. 1-based). Structural invariant checks are the primary defense.
- **Delegation:** LLM proficiency in boilerplate and frontend design isolates developer focus to core backend logic.

---

### Week 3

#### From the materials

- **Learned Representations:** Deep learning replaces hand-crafted rules with direct coordinate predictions, though accuracy degrades on novel or low-context data.
- **Joint Modeling:** Multi-modal foundation models enable conditional generation, utilizing embeddings that capture deep functional data over basic patterns.
- **Agent Orchestration:** Foundation models function as pipeline tools. System reliability depends on composition design, intermediate error inspection, and maintaining strict input/output contracts.

#### Surprises

Model/Agent: Antigravity Agent (Gemini 3.1 Low)

- **Embedding Configuration:** Default token/layer choices often yield poor results. Optimal layer selection and pooling strategies require explicit instruction and empirical comparison.
- **Padding Dilution:** Failing to mask padding tokens silently degrades embeddings for variable-length inputs. Plotting data length against embedding norms is a critical diagnostic hook.
- **Quality Blindspots:** Agents assess pipelines purely on execution success (pass/fail) and ignore output statistical confidence unless explicitly programmed to evaluate it.


### Week 4

#### From the materials

- **MCP (Model Context Protocol):** Wrapping tools as typed, named interfaces lets agents call them with structured args and get structured responses — no fragile stdout parsing. The leap from "agent runs shell commands" to "agent calls structured tools" is what makes workflows composable and auditable.
- **Three Modes Heuristic:** One-off exploration → direct commands; reusable analysis → agent-written code; production/repeated workflow → MCP. Choosing the wrong mode is where most friction comes from.
- **BixBench (Mitchener et al., 2025):** Frontier agents hit ~17% accuracy on real computational biology tasks. Useful as calibrated tools, not autonomous bioinformaticians. The honest numbers confirm the calibration lesson from weeks 1–3: verify everything, trust selectively.

#### Surprises

<!-- BioTerm-Bench, MCP demo, failure modes -->
