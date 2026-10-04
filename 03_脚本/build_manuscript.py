#!/usr/bin/env python3
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Inches, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "latex" / "mdpi" / "figures"
OUT = ROOT / "latex" / "mdpi" / "manuscript.docx"


def font(run, size=12, bold=False, italic=False):
    run.bold = bold
    run.italic = italic
    run.font.size = Pt(size)
    run.font.color.rgb = RGBColor(0, 0, 0)
    rpr = run._element.get_or_add_rPr()
    fonts = rpr.find(qn("w:rFonts"))
    if fonts is None:
        fonts = OxmlElement("w:rFonts")
        rpr.insert(0, fonts)
    for key in ("w:ascii", "w:hAnsi", "w:cs", "w:eastAsia"):
        fonts.set(qn(key), "Times New Roman")


def add_marked(paragraph, text, size=12, bold=False):
    for i, part in enumerate(text.split("*")):
        if not part:
            continue
        run = paragraph.add_run(part)
        font(run, size=size, bold=bold, italic=i % 2 == 1)


def block(doc, text, size=12, bold=False, center=False, before=0, after=6, indent=False):
    paragraph = doc.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER if center else WD_ALIGN_PARAGRAPH.JUSTIFY
    fmt = paragraph.paragraph_format
    fmt.space_before = Pt(before)
    fmt.space_after = Pt(after)
    fmt.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    fmt.line_spacing = 1.15
    fmt.first_line_indent = Cm(0.74) if indent else Cm(0)
    add_marked(paragraph, text, size=size, bold=bold)
    return paragraph


def heading(doc, text):
    return block(doc, text, bold=True, before=12, after=6)


def edge(name, val, sz="8"):
    el = OxmlElement(f"w:{name}")
    el.set(qn("w:val"), val)
    el.set(qn("w:sz"), "0" if val == "nil" else sz)
    el.set(qn("w:space"), "0")
    el.set(qn("w:color"), "auto")
    return el


def paint(cell, kind):
    tc_pr = cell._tc.get_or_add_tcPr()
    old = tc_pr.find(qn("w:tcBorders"))
    if old is not None:
        tc_pr.remove(old)
    borders = OxmlElement("w:tcBorders")
    borders.append(edge("left", "nil"))
    borders.append(edge("right", "nil"))
    if kind == "header":
        borders.append(edge("top", "single", "12"))
        borders.append(edge("bottom", "single", "6"))
    elif kind == "last":
        borders.append(edge("top", "nil"))
        borders.append(edge("bottom", "single", "12"))
    else:
        borders.append(edge("top", "nil"))
        borders.append(edge("bottom", "nil"))
    tc_pr.append(borders)


def add_table(doc, headers, rows):
    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.autofit = True
    tbl_pr = table._tbl.tblPr
    jc = OxmlElement("w:jc")
    jc.set(qn("w:val"), "center")
    tbl_pr.append(jc)
    borders = OxmlElement("w:tblBorders")
    for name, val, sz in (
        ("top", "single", "12"),
        ("left", "nil", "0"),
        ("bottom", "single", "12"),
        ("right", "nil", "0"),
        ("insideH", "nil", "0"),
        ("insideV", "nil", "0"),
    ):
        borders.append(edge(name, val, sz))
    tbl_pr.append(borders)
    for j, header in enumerate(headers):
        paragraph = table.rows[0].cells[j].paragraphs[0]
        paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
        run = paragraph.add_run(header)
        font(run, size=9, bold=True)
        paint(table.rows[0].cells[j], "header")
    for i, row in enumerate(rows):
        kind = "last" if i == len(rows) - 1 else "body"
        for j, value in enumerate(row):
            paragraph = table.rows[i + 1].cells[j].paragraphs[0]
            paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
            run = paragraph.add_run(value)
            font(run, size=9)
            paint(table.rows[i + 1].cells[j], kind)
    doc.add_paragraph().paragraph_format.space_after = Pt(6)


def add_figure(doc, path, caption):
    paragraph = doc.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.space_before = Pt(6)
    paragraph.paragraph_format.space_after = Pt(2)
    paragraph.add_run().add_picture(str(path), width=Inches(5.8))
    block(doc, caption, size=10, after=8)


def main():
    doc = Document()
    section = doc.sections[0]
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    for margin in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
        setattr(section, margin, Cm(2.5))
    normal = doc.styles["Normal"]
    normal.font.name = "Times New Roman"
    normal.font.size = Pt(12)
    core = doc.core_properties
    core.title = "Peritoneal C3 in public gastric-cancer samples"
    core.author = ""

    block(
        doc,
        "Peritoneal C3 Is Higher in Public Gastric-Cancer Samples, without a Replicated Ligand–Receptor Excess",
        size=16,
        bold=True,
        center=True,
        after=8,
    )
    block(doc, "Author names and affiliations to be completed.", size=11, center=True, after=10)
    heading(doc, "Abstract")
    block(
        doc,
        "Nine primary Visium sections did not retain a fibroblast–macrophage ligand–receptor excess. "
        "In GSE183904, two post-hoc C3 × C3AR1 class-mean products were higher in 3 peritoneal tumors than in 26 primary tumors (q = 0.027 within seven products). "
        "The increase was carried by C3, not by C3AR1. The SPP1 × CD44 macrophage product was not higher (q = 0.490). "
        "Macrophage STAT3 transcript was higher and SOCS3 was not. The same seven products were not higher in GSE308231 (all q > 0.35). "
        "Three peritoneal tumors do not establish a stable immunosuppressive axis.",
        indent=True,
    )
    block(doc, "Keywords: gastric cancer; peritoneal metastasis; Visium; C3; C3AR1", size=11, after=8)

    heading(doc, "1. Introduction")
    block(doc, "Peritoneal metastasis limits treatment in gastric cancer [1]. Public data do not contain a large paired primary–peritoneal spatial series. GSE251950 provides Visium sections from primary tumors and one metastasis that is not recorded as peritoneum [2,3]. GSE183904 has 3 peritoneal tumors and 26 primary tumors, including one paired donor [4].", indent=True)
    block(doc, "A locked 150 μm test of the nine primary sections retained no fibroblast–macrophage pair. The question asked here is different: whether patient-level products of C3 or SPP1 with their receptors are higher in peritoneal tumor. The seven products were named after the class means had been seen. GSE308231 was a separate 3-versus-3 test. Peritoneal fluid was not treated as solid tumor.", indent=True)

    heading(doc, "2. Results")
    heading(doc, "2.1. Primary Visium")
    block(doc, "Contact of macrophage-high spots with a fibroblast-high neighbor was 0.81–0.94, below the label-shuffle 95th percentile (0.95–0.97) on every section. None of 1000 single-gene pairs met both the 5-of-9 rule and leave-one-out. The five prespecified pairs had a nine-section median of zero.", indent=True)
    block(doc, "A proportion-weighted mean removed the all-zero summaries but not the negative excess (Figure 1). At 150 μm, median excesses were −0.068 for APP–CD74, −0.005 for C3–C3AR1 (below the null on 8 of 9), −0.018 for SPP1–CD44 (9 of 9), and −0.001 for CCL2–CCR2. Excess medians at 100, 200, and 300 μm stayed at or below zero, except TGFB1–TGFBR1 at 100 μm (+0.0001). On the unspecified metastasis, the C3–C3AR1 weighted mean was 0.922 against a null 95th percentile of 0.925.", indent=True)
    add_figure(
        doc,
        FIG / "fig_interface.png",
        "Figure 1. Primary Visium sections at 150 μm. (A) Weighted mean divided by the ligand-shift null 95th percentile. The black bar is the nine-section median. The dashed line is 1. (B) Blue points, fibroblast–macrophage correlation. Red circles, null 95th percentile.",
    )

    heading(doc, "2.2. GSE183904")
    block(doc, "Primary sample 22 lacked enough fibroblasts and epithelial cells, so those denominators were 25. Macrophage and T-cell denominators were 26. Figure 2 and Table 1 give the seven products. Both C3 products use macrophage C3AR1. Epithelial C3, fibroblast C3, and macrophage C3AR1 had median ratios of about 30.5, 5.4, and 1.7. Within 25 primary tumors, fibroblast C3 and macrophage C3AR1 had a Spearman correlation of 0.01.", indent=True)
    add_figure(
        doc,
        FIG / "fig_communication.png",
        "Figure 2. GSE183904 class-mean products. Blue, primary tumors. Red, peritoneal tumors. The black bar is the primary median.",
    )
    block(doc, "Table 1. GSE183904. One-sided label permutation. q is Benjamini–Hochberg within these seven products.", size=10, bold=True, after=2)
    add_table(
        doc,
        ["Product", "Peritoneal", "Primary", "Above", "p", "q"],
        [
            ["C3 fib × C3AR1 mac", "0.730", "0.132", "3/3", "0.015", "0.027"],
            ["C3 epi × C3AR1 mac", "0.542", "0.007", "3/3", "0.004", "0.027"],
            ["SPP1 mac × CD44 mac", "0.505", "0.242", "2/3", "0.420", "0.490"],
            ["SPP1 mac × CD44 T", "0.964", "0.172", "3/3", "0.094", "0.131"],
            ["CCL2 epi × CCR2 T", "0.021", "0.001", "3/3", "0.015", "0.027"],
            ["CCL2 fib × CCR2 mac", "0.009", "0.015", "1/3", "0.797", "0.797"],
            ["TGFB1 T × TGFBR2 mac", "0.332", "0.136", "3/3", "0.014", "0.027"],
        ],
    )
    block(doc, "Macrophage STAT3 medians were 0.641 and 0.481 (p = 0.010). SOCS3 was 0.695 versus 0.852. Their mean was 0.668 versus 0.623 (p = 0.47). Epithelial STAT3 was also higher (p = 0.015). For the paired donor, sample 38 minus sample 36 was 1.177 and 1.200 for the two C3 products.", indent=True)

    heading(doc, "2.3. Three Other Public Series")
    block(doc, "GSE308231 did not reproduce the products (Figure 3, Table 2). All three peritoneal C3–epithelial products were above all three primary products, but another labeling gave the same median difference, so p = 0.10 rather than 0.05. Fibroblast C3 was higher in every peritoneal sample than in every primary sample. Macrophage C3AR1 medians were 0.395 and 0.648. Macrophage STAT3 was 0.697 versus 0.864, and every peritoneal SOCS3 value was below every primary value.", indent=True)
    add_figure(
        doc,
        FIG / "fig_public_pm.png",
        "Figure 3. GSE308231. Blue, primary tumors. Red, peritoneal tumors.",
    )
    block(doc, "Table 2. GSE308231, three versus three. The smallest one-sided p is 0.05.", size=10, bold=True, after=2)
    add_table(
        doc,
        ["Product", "Peritoneal", "Primary", "p", "q"],
        [
            ["C3 fib × C3AR1 mac", "0.375", "0.410", "0.70", "0.70"],
            ["C3 epi × C3AR1 mac", "0.119", "0.026", "0.10", "0.35"],
            ["SPP1 mac × CD44 mac", "0.293", "0.200", "0.50", "0.58"],
            ["SPP1 mac × CD44 T", "0.298", "0.200", "0.30", "0.42"],
            ["CCL2 epi × CCR2 T", "0.006", "0.001", "0.10", "0.35"],
            ["CCL2 fib × CCR2 mac", "0.036", "0.016", "0.20", "0.35"],
            ["TGFB1 T × TGFBR2 mac", "0.664", "0.492", "0.20", "0.35"],
        ],
    )
    block(doc, "GSE163558 has one peritoneal sample. Its C3–epithelial product was 0.294 against a primary median of 0.042. Its macrophage STAT3, 0.590, was the lowest of the four samples. No p value was assigned. In GSE228598, 15 of 28 fluid samples had fewer than 20 fibroblasts. Among the other 13, median fibroblast C3 was 2.426. Across all 28, median epithelial C3 was 0.924 and macrophage C3AR1 was 0.829. Fluid was not entered into the solid-tumor tests.", indent=True)

    heading(doc, "3. Discussion")
    block(doc, "The Visium result follows the locked rule. A 150 μm neighborhood on spots about 100 μm apart, with a top-quartile label, puts the contact null near 0.97. High absolute products such as APP–CD74 stayed high after the coordinates were shuffled.", indent=True)
    block(doc, "In GSE183904 the two C3 products have q = 0.027 inside a list chosen after the class means were known. GSE308231 repeats the ligand pattern and not the product. SPP1 × CD44 fails in both series. STAT3 transcript and SOCS3 do not move together, and the STAT3 difference is not confined to macrophages. This is not evidence that C3–C3AR1 drives phosphorylated STAT3. Fluid and a solid implant are different materials, and no protein stain or blockade was done.", indent=True)

    heading(doc, "4. Materials and Methods")
    heading(doc, "4.1. Data")
    block(doc, "Nine primary Visium sections from GSE251950 entered the spatial summary. GC6-PM (GSM7990479) was scored alone and was not called peritoneum [3]. GSE183904 contributed peritoneal tumors GSM5573484, GSM5573485, and GSM5573503, and 26 primary tumors [4]. Sample 38 and sample 36 are the paired donor. GSE308231 has three primary and three peritoneal tumors and is not stated as paired [5]. GSE163558 contributes one peritoneal sample against three primary tumors [6]. GSE228598 has 28 ascites or washing samples and no primary arm [7].", indent=True)
    heading(doc, "4.2. Spatial Readouts")
    block(doc, "Spots were in tissue, had at least 200 genes, and had mitochondrial fraction at most 25%. Counts were normalized to 10,000 and log-transformed. Fibroblast markers were DCN, LUM, COL1A1, and FAP. Macrophage markers were CD68, CSF1R, C1QA, and C1QB. High spots were at or above the section 75th percentile. The spot signal was receptor expression times the mean ligand on fibroblast-high neighbors within 150 μm. The section summary was the median. Contact and coordinate nulls used 200 shuffles, seed 20261002. Retention required at least 5 of 9 primary sections above the null, and the same bar after leaving out any one section. The pair list was the CellChat human table, single-gene symbols only, 1000 pairs [8].", indent=True)
    block(doc, "The weighted readout decomposed each spot by non-negative least squares into fibroblast, macrophage, epithelial, and T-cell programs [9]. The weight was the product of the fibroblast and macrophage proportions. The section score was the weighted mean. The null shifted the ligand map by 200–500 μm, 200 times, seed 20261004. This readout did not replace the retention rule.", indent=True)
    heading(doc, "4.3. Class-Mean Products")
    block(doc, "A cell was assigned to the highest of the four programs when that mean was positive and strictly above the other three. A class with fewer than 20 cells was missing. The seven products were C3 fibroblast or epithelial × C3AR1 macrophage, SPP1 macrophage × CD44 macrophage or T cell, CCL2 epithelial × CCR2 T cell, CCL2 fibroblast × CCR2 macrophage, and TGFB1 T cell × TGFBR2 macrophage. The statistic was the peritoneal median minus the primary median. q values used the Benjamini–Hochberg procedure within these seven tests [10]. STAT3 and SOCS3 were transcript means and were not entered into that q value.", indent=True)

    heading(doc, "5. Conclusions")
    block(doc, "Primary Visium sections did not retain a fibroblast–macrophage excess. In GSE183904, C3 products with macrophage C3AR1 were higher in three peritoneal tumors, mainly because C3 was higher. Those products were not higher in GSE308231.", indent=True)

    heading(doc, "Author Contributions")
    block(doc, "To be completed by the authors.")
    heading(doc, "Funding")
    block(doc, "To be completed by the authors.")
    heading(doc, "Institutional Review Board Statement")
    block(doc, "Public, de-identified matrices from the Gene Expression Omnibus were used. No new patients were enrolled.")
    heading(doc, "Informed Consent Statement")
    block(doc, "Not applicable.")
    heading(doc, "Data Availability Statement")
    block(doc, "GSE251950, GSE183904, GSE308231, GSE163558, and GSE228598 are available from the Gene Expression Omnibus. Summary tables are available from the authors. The analysis script is not in a public repository.")
    heading(doc, "Conflicts of Interest")
    block(doc, "To be completed by the authors.")

    heading(doc, "References")
    refs = [
        "1. Manzanedo, I.; Pereira, F.; Pérez-Viejo, E.; Serrano, Á. Gastric cancer with peritoneal metastases: Current status and prospects for treatment. Cancers 2023, 15, 1777.",
        "2. Ståhl, P.L.; Salmén, F.; Vickovic, S.; Lundmark, A.; Navarro, J.F.; Magnusson, J.; et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 2016, 353, 78–82.",
        "3. Lee, S.; Lee, D.; Choi, J.; Oh, H.; Ham, I.; Ryu, D.; et al. Spatial dissection of tumour microenvironments in gastric cancers reveals the immunosuppressive crosstalk between CCL2+ fibroblasts and STAT3-activated macrophages. Gut 2025, 74, 714–727.",
        "4. Kumar, V.; Ramnarayanan, K.; Sundar, R.; Padmanabhan, N.; Srivastava, S.; Koiwa, M.; et al. Single-cell atlas of lineage states, tumor microenvironment, and subtype-specific expression programs in gastric cancer. Cancer Discov. 2022, 12, 670–691.",
        "5. Gene Expression Omnibus. GSE308231. Available online: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308231 (accessed on 4 October 2026).",
        "6. Jiang, H.; Yu, D.; Yang, P.; Guo, R.; Kong, M.; Gao, Y.; Yu, X.; Lu, X.; Fan, X. Revealing the transcriptional heterogeneity of organ-specific metastasis in human gastric cancer using single-cell RNA sequencing. Clin. Transl. Med. 2022, 12, e730.",
        "7. Sullivan, K.M.; Li, H.; Yang, A.; Zhang, Z.; Munoz, R.R.; Mahuron, K.M.; et al. Tumor and peritoneum-associated macrophage gene signature as a novel molecular biomarker in gastric cancer. Int. J. Mol. Sci. 2024, 25, 4117.",
        "8. Jin, S.; Guerrero-Juarez, C.F.; Zhang, L.; Chang, I.; Ramos, R.; Kuan, C.H.; et al. Inference and analysis of cell-cell communication using CellChat. Nat. Commun. 2021, 12, 1088.",
        "9. Lawson, C.L.; Hanson, R.J. Solving Least Squares Problems; Prentice-Hall: Englewood Cliffs, NJ, USA, 1974.",
        "10. Benjamini, Y.; Hochberg, Y. Controlling the false discovery rate: A practical and powerful approach to multiple testing. J. R. Stat. Soc. B 1995, 57, 289–300.",
    ]
    for ref in refs:
        paragraph = block(doc, ref, size=10, after=3)
        paragraph.paragraph_format.left_indent = Cm(0.75)
        paragraph.paragraph_format.first_line_indent = Cm(-0.75)
        paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
    doc.save(OUT)
    print(OUT)


if __name__ == "__main__":
    main()
