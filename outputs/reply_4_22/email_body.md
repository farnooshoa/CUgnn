# Reply email — 2026-04-22

**To:** Dr. Schwartz-Duval
**Subject:** Re: Histones, Au compounds, cartoon updates
**Attachments:** `copper_cell_cartoon_network.html`, `glossary.md`

---

Dear Dr. Schwartz-Duval,

Thanks for the 4 references. I've worked through them and updated the cartoon — new version attached.

**PRNP and ENOX2 are already nodes in the 58-gene graph**, so no retraining was needed. Their current ranks in the existing GAT readout:

- Saliency: PRNP rank 26 / 58, ENOX1 rank 40 / 58, ENOX2 rank 43 / 58.
- Attention edges (148 total): ATP7B ↔ PRNP ranks 15 (top 10%), PRNP ↔ SNCA / PRNP ↔ APP rank 33 and 34, ATP7A ↔ PRNP rank 122. ENOX1 / ENOX2 show strong self-attention, which the model does when a node has few informative neighbours — adding Au lines to both helps fill that gap.

**Save / Load layout added.** In the new HTML there are two buttons in the sidebar:

- **Save layout** downloads the current node positions as a JSON file.
- **Load layout** reads a JSON file and applies its positions.

Once you have the figure arranged the way you want it, please send me the saved JSON and I will bake it in as the new default layout for the cartoon (so the built-in Reset goes to your arrangement, not mine). That way any future version keeps your figure layout without you redoing the drag-work.

**Au compounds node added**, placed extracellularly on the right. Connecting lines:

- **Solid (published)**: ATOX1 (Gabbiani 2012, Marzo 2015), ATP7A and ATP7B (Spreckelmeyer 2018), LOX (Georgiou 2011). Hovering each line shows the citation and one-sentence mechanism.
- **Dashed (provisional)**: H3-3A, H3-3B, H3C1, H4C1, PRNP, ENOX2. The 4 references you sent do not cover these targets, so I drew them dashed. If you can share the unpublished data you mentioned for histones and for PRNP / ENOX2, I'll promote them to solid and add the references.

Two things I want to flag while you look at the figure:

- **LOX naming.** The 2011 paper tests lipoxygenase (the ALOX class, linoleic-acid oxygenase). Our graph LOX is lysyl oxidase, a different Cu enzyme (ECM crosslinker, family LOX / LOXL1–4). I drew the Au–LOX line provisionally because you listed LOX in your email. If you meant lipoxygenase, we should add ALOX5 / ALOX12 / ALOX15 as separate nodes.
- **SLC31A1 / CTR1.** The Spreckelmeyer 2018 paper explicitly rules out CTR1 as an uptake route for their Au(III) compound, so I did not draw an Au–CTR1 line. I can add it back if the scope is broader than that paper.

Histone coverage. The current 58-node graph contains H3 and H4 representatives (H3-3A, H3-3B, H3C1, H4C1). You mentioned H1 and H2A / H2B — those are not in the current graph. Happy to add them if your data supports Au interaction.

**Video call.** Yes, happy to meet. Please send a few time windows that work for you this week or next. To make the call faster I've attached a short glossary (`glossary.md`) covering the ML and statistical terms from the earlier materials — GAT, saliency, attention, StratifiedGroupKFold, and so on. Please glance through it beforehand and flag anything that still looks opaque, and I'll plan the call around those.

Best,
Ruiheng
