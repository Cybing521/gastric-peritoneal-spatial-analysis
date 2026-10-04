#!/usr/bin/env python3
import csv
from itertools import combinations
from pathlib import Path
from statistics import median

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Cm, Inches, Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "latex" / "mdpi" / "figures"
RES = ROOT / "06_真实分析结果"
OUT = ROOT / "latex" / "mdpi" / "manuscript.docx"
FONT = "Palatino"


def font(run, size=11, bold=False, italic=False):
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
        fonts.set(qn(key), FONT)


def add_runs(paragraph, parts, size=11):
    for text, bold, italic in parts:
        if not text:
            continue
        run = paragraph.add_run(text)
        font(run, size=size, bold=bold, italic=italic)


def block(doc, text, size=11, bold=False, center=False, before=0, after=6, indent=False, left=False, page_break=False):
    paragraph = doc.add_paragraph()
    if center:
        paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    elif left:
        paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
    else:
        paragraph.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    fmt = paragraph.paragraph_format
    fmt.space_before = Pt(before)
    fmt.space_after = Pt(after)
    fmt.line_spacing_rule = WD_LINE_SPACING.MULTIPLE
    fmt.line_spacing = 1.15
    fmt.first_line_indent = Cm(0.75) if indent else Cm(0)
    fmt.page_break_before = page_break
    add_runs(paragraph, [(text, bold, False)], size=size)
    return paragraph


def heading(doc, text, size=12):
    return block(doc, text, size=size, bold=True, before=12, after=4)


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


def shade_off(cell):
    tc_pr = cell._tc.get_or_add_tcPr()
    for tag in ("w:shd", "w:tcMar"):
        old = tc_pr.find(qn(tag))
        if old is not None:
            tc_pr.remove(old)


def add_table(doc, headers, rows, widths=None):
    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.autofit = False
    table.allow_autofit = False
    tbl_pr = table._tbl.tblPr
    jc = OxmlElement("w:jc")
    jc.set(qn("w:val"), "center")
    tbl_pr.append(jc)
    bidi = tbl_pr.find(qn("w:bidiVisual"))
    if bidi is not None:
        tbl_pr.remove(bidi)
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
        cell = table.rows[0].cells[j]
        shade_off(cell)
        paragraph = cell.paragraphs[0]
        paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
        paragraph.paragraph_format.space_before = Pt(1)
        paragraph.paragraph_format.space_after = Pt(1)
        run = paragraph.add_run(header)
        font(run, size=8, bold=True)
        paint(cell, "header")
    tr_pr = table.rows[0]._tr.get_or_add_trPr()
    tr_pr.append(OxmlElement("w:tblHeader"))
    for i, row in enumerate(rows):
        kind = "last" if i == len(rows) - 1 else "body"
        for j, value in enumerate(row):
            cell = table.rows[i + 1].cells[j]
            shade_off(cell)
            paragraph = cell.paragraphs[0]
            paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
            paragraph.paragraph_format.space_before = Pt(1)
            paragraph.paragraph_format.space_after = Pt(1)
            run = paragraph.add_run(value)
            font(run, size=8)
            paint(cell, kind)
    if widths:
        total = sum(widths)
        tbl_w = OxmlElement("w:tblW")
        tbl_w.set(qn("w:w"), str(int(total * 567)))
        tbl_w.set(qn("w:type"), "dxa")
        tbl_pr.append(tbl_w)
        grid = table._tbl.find(qn("w:tblGrid"))
        if grid is not None:
            table._tbl.remove(grid)
        grid = OxmlElement("w:tblGrid")
        for width in widths:
            col = OxmlElement("w:gridCol")
            col.set(qn("w:w"), str(int(width * 567)))
            grid.append(col)
        table._tbl.insert(1, grid)
        for row in table.rows:
            for cell, width in zip(row.cells, widths):
                cell.width = Cm(width)
                tc_w = cell._tc.get_or_add_tcPr().find(qn("w:tcW"))
                if tc_w is None:
                    tc_w = OxmlElement("w:tcW")
                    cell._tc.get_or_add_tcPr().append(tc_w)
                tc_w.set(qn("w:w"), str(int(width * 567)))
                tc_w.set(qn("w:type"), "dxa")
    for row in table.rows:
        tr_pr = row._tr.get_or_add_trPr()
        tr_pr.append(OxmlElement("w:cantSplit"))
    doc.add_paragraph().paragraph_format.space_after = Pt(6)


def add_figure(doc, path, caption):
    paragraph = doc.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.space_before = Pt(8)
    paragraph.paragraph_format.space_after = Pt(2)
    run = paragraph.add_run()
    run.add_picture(str(path), width=Inches(6.3))
    cap = block(doc, caption, size=9, after=8)
    cap.alignment = WD_ALIGN_PARAGRAPH.LEFT
    return cap


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def f3(value):
    return f"{float(value):.3f}"


def f2(value):
    return f"{float(value):.2f}"


def one_sided(pm, pri):
    obs = median(pm) - median(pri)
    pool = list(pm) + list(pri)
    k = len(pm)
    hit = 0
    total = 0
    for idx in combinations(range(len(pool)), k):
        total += 1
        draw = [pool[i] for i in idx]
        rest = [pool[i] for i in range(len(pool)) if i not in idx]
        if median(draw) - median(rest) >= obs - 1e-12:
            hit += 1
    return obs, hit / total, total


def person(raw):
    suffix = ""
    for token in (" 2nd", " 3rd", " 4th"):
        if raw.endswith(token):
            suffix = "," + token
            raw = raw[: -len(token)]
            break
    parts = raw.split()
    last = parts[-1]
    if last.isalpha() and last.isupper() and len(last) <= 4 and len(parts) > 1:
        given = ".".join(last) + "."
        return " ".join(parts[:-1]) + ", " + given + suffix
    return raw + suffix


def author_line(names):
    shown = [person(name) for name in names[:10]]
    text = "; ".join(shown)
    if len(names) > 10:
        text += "; et al."
    return text


def add_ref(doc, number, authors, title, journal, year, volume, pages, doi):
    paragraph = doc.add_paragraph()
    paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
    fmt = paragraph.paragraph_format
    fmt.space_before = Pt(0)
    fmt.space_after = Pt(2)
    fmt.line_spacing = 1.08
    fmt.left_indent = Cm(0.75)
    fmt.first_line_indent = Cm(-0.75)
    bits = [(f"{number}. ", False, False), (authors + " ", False, False), (title + " ", False, False)]
    if journal:
        bits.append((journal + " ", False, True))
        bits.append((f"{year}", True, False))
        bits.append((f", {volume}, {pages}. ", False, False))
    else:
        bits.append((title and "" or "", False, False))
    if doi:
        bits.append((doi, False, False))
    add_runs(paragraph, bits, size=9)


REFS = [
    (["Bray F", "Laversanne M", "Sung H", "Ferlay J", "Siegel RL", "Soerjomataram I", "Jemal A"],
     "Global cancer statistics 2022: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries.",
     "CA Cancer J. Clin.", "2024", "74", "229–263", "https://doi.org/10.3322/caac.21834"),
    (["Smyth EC", "Nilsson M", "Grabsch HI", "van Grieken NC", "Lordick F"],
     "Gastric cancer.",
     "Lancet", "2020", "396", "635–648", "https://doi.org/10.1016/S0140-6736(20)31288-5"),
    (["Manzanedo I", "Pereira F", "Pérez-Viejo E", "Serrano Á"],
     "Gastric cancer with peritoneal metastases: Current status and prospects for treatment.",
     "Cancers", "2023", "15", "1777", "https://doi.org/10.3390/cancers15061777"),
    (["Li Z", "Wang J", "Wang Z", "Xu Y"],
     "Towards an optimal model for gastric cancer peritoneal metastasis: Current challenges and future directions.",
     "eBioMedicine", "2023", "92", "104601", "https://doi.org/10.1016/j.ebiom.2023.104601"),
    (["Gwee YX", "Chia DKA", "So J", "Ceelen W", "Yong WP", "Tan P", "Ong CJ", "Sundar R"],
     "Integration of genomic biology into therapeutic strategies of gastric cancer peritoneal metastasis.",
     "J. Clin. Oncol.", "2022", "40", "2830", "https://doi.org/10.1200/JCO.21.02745"),
    (["Yao X", "Ajani JA", "Song S"],
     "Molecular biology and immunology of gastric cancer peritoneal metastasis.",
     "Transl. Gastroenterol. Hepatol.", "2020", "5", "57", "https://doi.org/10.21037/tgh.2020.02.08"),
    (["Cancer Genome Atlas Research Network"],
     "Comprehensive molecular characterization of gastric adenocarcinoma.",
     "Nature", "2014", "513", "202–209", "https://doi.org/10.1038/nature13480"),
    (["Binnewies M", "Roberts EW", "Kersten K", "Chan V", "Fearon DF", "Merad M", "Coussens LM", "Gabrilovich DI", "Ostrand-Rosenberg S", "Hedrick CC", "Vonderheide RH", "Pittet MJ"],
     "Understanding the tumor immune microenvironment (TIME) for effective therapy.",
     "Nat. Med.", "2018", "24", "541–550", "https://doi.org/10.1038/s41591-018-0014-x"),
    (["Sahai E", "Astsaturov I", "Cukierman E", "DeNardo DG", "Egeblad M", "Evans RM", "Fearon D", "Greten FR", "Hingorani SR", "Hunter T", "Hynes RO", "Jain RK"],
     "A framework for advancing our understanding of cancer-associated fibroblasts.",
     "Nat. Rev. Cancer", "2020", "20", "174–186", "https://doi.org/10.1038/s41568-019-0238-1"),
    (["Ham IH", "Lee D", "Hur H"],
     "Role of cancer-associated fibroblast in gastric cancer progression and resistance to treatments.",
     "J. Oncol.", "2019", "2019", "6270784", "https://doi.org/10.1155/2019/6270784"),
    (["Li X", "Sun Z", "Peng G", "Xiao Y", "Guo J", "Wu B", "Li X", "Zhou W", "Li J", "Li Z", "Bai C", "Zhao L"],
     "Single-cell RNA sequencing reveals a pro-invasive cancer-associated fibroblast subgroup associated with poor clinical outcomes in patients with gastric cancer.",
     "Theranostics", "2022", "12", "620–638", "https://doi.org/10.7150/thno.60540"),
    (["Mantovani A", "Marchesi F", "Malesci A", "Laghi L", "Allavena P"],
     "Tumour-associated macrophages as treatment targets in oncology.",
     "Nat. Rev. Clin. Oncol.", "2017", "14", "399–416", "https://doi.org/10.1038/nrclinonc.2016.217"),
    (["Roumenina LT", "Daugan MV", "Petitprez F", "Sautès-Fridman C", "Fridman WH"],
     "Context-dependent roles of complement in cancer.",
     "Nat. Rev. Cancer", "2019", "19", "698–715", "https://doi.org/10.1038/s41568-019-0210-0"),
    (["Bill R", "Wirapati P", "Messemaker M", "Roh W", "Zitti B", "Duval F", "Kiss M", "Park JC", "Saal TM", "Hoelzl J", "Tarussio D", "Benedetti F"],
     "CXCL9:SPP1 macrophage polarity identifies a network of cellular programs that control human cancers.",
     "Science", "2023", "381", "515–524", "https://doi.org/10.1126/science.ade2292"),
    (["Liu F", "Zhang J", "Gu X", "Guo Q", "Guo W"],
     "Single-cell transcriptome sequencing reveals SPP1-CD44-mediated macrophage–tumor cell interactions drive chemoresistance in TNBC.",
     "J. Cell. Mol. Med.", "2024", "28", "e18525", "https://doi.org/10.1111/jcmm.18525"),
    (["Yang Y", "Sun H", "Yu H", "Wang L", "Gao C", "Mei H", "Jiang X", "Ji M"],
     "Tumor-associated fibrosis and active collagen-CD44 axis characterize a poor-prognosis subtype of gastric cancer and contribute to tumor immunosuppression.",
     "J. Transl. Med.", "2025", "23", "123", "https://doi.org/10.1186/s12967-025-06070-9"),
    (["Zhang F", "Yan Y", "Cao X", "Guo C", "Wang K", "Lv S"],
     "TGF-β-driven LIF expression influences neutrophil extracellular traps (NETs) and contributes to peritoneal metastasis in gastric cancer.",
     "Cell Death Dis.", "2024", "15", "218", "https://doi.org/10.1038/s41419-024-06594-w"),
    (["Yu H", "Lee H", "Herrmann A", "Buettner R", "Jove R"],
     "Revisiting STAT3 signalling in cancer: New and unexpected biological functions.",
     "Nat. Rev. Cancer", "2014", "14", "736–746", "https://doi.org/10.1038/nrc3818"),
    (["Ashrafizadeh M", "Zarrabi A", "Orouei S", "Zarrin V", "Rahmani Moghadam E", "Zabolian A", "Mohammadi S", "Hushmandi K", "Gharehaghajlou Y", "Makvandi P", "Najafi M", "Mohammadinejad R"],
     "STAT3 pathway in gastric cancer: Signaling, therapeutic targeting and future prospects.",
     "Biology", "2020", "9", "126", "https://doi.org/10.3390/biology9060126"),
    (["Su P", "Yu T", "Zhang Y", "Huang H", "Chen M", "Cao C", "Kang W", "Liu Y", "Yu J"],
     "Upregulation of MELK promotes chemoresistance and induces macrophage M2 polarization via CSF-1/JAK2/STAT3 pathway in gastric cancer.",
     "Cancer Cell Int.", "2024", "24", "287", "https://doi.org/10.1186/s12935-024-03453-8"),
    (["Cui Y", "Chang Y", "Ma X", "Sun M", "Huang Y", "Yang F", "Li S", "Zhuo W", "Liu W", "Yang B", "Lin A", "Ou G"],
     "Ephrin A1 stimulates CCL2 secretion to facilitate premetastatic niche formation and promote gastric cancer liver metastasis.",
     "Cancer Res.", "2025", "85", "263–276", "https://doi.org/10.1158/0008-5472.CAN-24-1254"),
    (["Lee SH", "Lee D", "Choi J", "Oh HJ", "Ham IH", "Ryu D", "Lee SY", "Han DJ", "Kim S", "Moon Y", "Song IH", "Song KY"],
     "Spatial dissection of tumour microenvironments in gastric cancers reveals the immunosuppressive crosstalk between CCL2+ fibroblasts and STAT3-activated macrophages.",
     "Gut", "2025", "74", "714–727", "https://doi.org/10.1136/gutjnl-2024-332901"),
    (["Ståhl PL", "Salmén F", "Vickovic S", "Lundmark A", "Navarro JF", "Magnusson J", "Giacomello S", "Asp M", "Westholm JO", "Huss M", "Mollbrink A", "Linnarsson S"],
     "Visualization and analysis of gene expression in tissue sections by spatial transcriptomics.",
     "Science", "2016", "353", "78–82", "https://doi.org/10.1126/science.aaf2403"),
    (["Williams CG", "Lee HJ", "Asatsuma T", "Vento-Tormo R", "Haque A"],
     "An introduction to spatial transcriptomics for biomedical research.",
     "Genome Med.", "2022", "14", "68", "https://doi.org/10.1186/s13073-022-01075-1"),
    (["Moses L", "Pachter L"],
     "Museum of spatial transcriptomics.",
     "Nat. Methods", "2022", "19", "534–546", "https://doi.org/10.1038/s41592-022-01409-2"),
    (["Maynard KR", "Collado-Torres L", "Weber LM", "Uytingco C", "Barry BK", "Williams SR", "Catallini JL 2nd", "Tran MN", "Besich Z", "Tippani M", "Chew J", "Yin Y"],
     "Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex.",
     "Nat. Neurosci.", "2021", "24", "425–436", "https://doi.org/10.1038/s41593-020-00787-0"),
    (["Kleshchevnikov V", "Shmatko A", "Dann E", "Aivazidis A", "King HW", "Li T", "Elmentaite R", "Lomakin A", "Kedlian V", "Gayoso A", "Jain MS", "Park JS"],
     "Cell2location maps fine-grained cell types in spatial transcriptomics.",
     "Nat. Biotechnol.", "2022", "40", "661–671", "https://doi.org/10.1038/s41587-021-01139-4"),
    (["Kumar V", "Ramnarayanan K", "Sundar R", "Padmanabhan N", "Srivastava S", "Koiwa M", "Yasuda T", "Koh V", "Huang KK", "Tay ST", "Ho SWT", "Tan ALK"],
     "Single-cell atlas of lineage states, tumor microenvironment, and subtype-specific expression programs in gastric cancer.",
     "Cancer Discov.", "2022", "12", "670–691", "https://doi.org/10.1158/2159-8290.CD-21-0683"),
    (["Jiang H", "Yu D", "Yang P", "Guo R", "Kong M", "Gao Y", "Yu X", "Lu X", "Fan X"],
     "Revealing the transcriptional heterogeneity of organ-specific metastasis in human gastric cancer using single-cell RNA sequencing.",
     "Clin. Transl. Med.", "2022", "12", "e730", "https://doi.org/10.1002/ctm2.730"),
    (["Sullivan KM", "Li H", "Yang A", "Zhang Z", "Munoz RR", "Mahuron KM", "Yuan YC", "Paz IB", "Von Hoff D", "Han H", "Fong Y", "Woo Y"],
     "Tumor and peritoneum-associated macrophage gene signature as a novel molecular biomarker in gastric cancer.",
     "Int. J. Mol. Sci.", "2024", "25", "4117", "https://doi.org/10.3390/ijms25074117"),
    (["Huang XZ", "Pang MJ", "Li JY", "Chen HY", "Sun JX", "Song YX", "Ni HJ", "Ye SY", "Bai S", "Li TH", "Wang XY", "Lu JY"],
     "Single-cell sequencing of ascites fluid illustrates heterogeneity and therapy-induced evolution during gastric cancer peritoneal metastasis.",
     "Nat. Commun.", "2023", "14", "822", "https://doi.org/10.1038/s41467-023-36310-9"),
    (["Kim H", "Kwon M", "Lee SK", "Son SM", "Lee OJ", "Man Yoon S", "Kim HK", "Yang Y", "Lee KH", "Han HS"],
     "Distinct immunosuppressive tumor microenvironment in gastric cancer with peritoneal metastasis.",
     "J. Gastric Cancer", "2025", "25", "605–620", "https://doi.org/10.5230/jgc.2025.25.e46"),
    (["Peng H", "Jiang L", "Yuan J", "Wu X", "Chen N", "Liu D", "Liang Y", "Xie Y", "Jia K", "Li Y", "Feng X", "Li J"],
     "Single-cell characterization of differentiation trajectories and drug resistance features in gastric cancer with peritoneal metastasis.",
     "Clin. Transl. Med.", "2024", "14", "e70054", "https://doi.org/10.1002/ctm2.70054"),
    ([], "", "", "", "", "", ""),
    (["Benjamini Y", "Hochberg Y"],
     "Controlling the false discovery rate: A practical and powerful approach to multiple testing.",
     "J. R. Stat. Soc. Ser. B (Stat. Methodol.)", "1995", "57", "289–300", "https://doi.org/10.1111/j.2517-6161.1995.tb02031.x"),
    (["Wolf FA", "Angerer P", "Theis FJ"],
     "SCANPY: Large-scale single-cell gene expression data analysis.",
     "Genome Biol.", "2018", "19", "15", "https://doi.org/10.1186/s13059-017-1382-0"),
    (["Jin S", "Guerrero-Juarez CF", "Zhang L", "Chang I", "Ramos R", "Kuan CH", "Myung P", "Plikus MV", "Nie Q"],
     "Inference and analysis of cell-cell communication using CellChat.",
     "Nat. Commun.", "2021", "12", "1088", "https://doi.org/10.1038/s41467-021-21246-9"),
]


def class_median(rows, gene, cls, col):
    for row in rows:
        if row["gene"] == gene and row["cell_class"] == cls:
            return float(row[col])
    raise KeyError((gene, cls))


def main():
    source = read_csv(RES / "界面深化" / "scrna_source_summary.csv")
    primary = read_csv(RES / "界面深化" / "primary_summary.csv")
    pairs = read_csv(RES / "pair_summary.csv")
    perm = read_csv(RES / "腹膜通讯" / "communication_permutation.csv")
    by_perm = {row["pair"]: row for row in perm}
    paired_c3_fib = float(by_perm["C3_fib_C3AR1_mac"]["paired_38_minus_36"])
    paired_c3_epi = float(by_perm["C3_epi_C3AR1_mac"]["paired_38_minus_36"])
    paired_ccl2 = float(by_perm["CCL2_epi_CCR2_t"]["paired_38_minus_36"])
    n_comb_fib = int(float(by_perm["C3_fib_C3AR1_mac"]["n_comb"]))
    n_comb_mac = int(float(by_perm["TGFB1_t_TGFBR2_mac"]["n_comb"]))
    if n_comb_fib != 3276 or n_comb_mac != 3654:
        raise SystemExit((n_comb_fib, n_comb_mac))
    perm308 = read_csv(RES / "公开腹膜" / "gse308231_permutation.csv")
    stat = read_csv(RES / "腹膜通讯" / "scrna_stat3_class.csv")
    fluid = read_csv(RES / "公开腹膜" / "gse228598_fluid_summary.csv")
    fluid_map = {row["measure"]: row for row in fluid}
    public = read_csv(RES / "公开腹膜" / "sample_class_means.csv")
    g308_rows = [row for row in public if row["dataset"] == "GSE308231"]
    c3_pm = [float(r["C3__fib"]) for r in g308_rows if r["role"] == "peritoneal_tumor"]
    c3_pri = [float(r["C3__fib"]) for r in g308_rows if r["role"] == "primary_tumor"]
    if not (min(c3_pm) > max(c3_pri)):
        raise SystemExit("GSE308231 fibroblast C3 is not separated")
    g163 = {row["pair"]: row for row in read_csv(RES / "公开腹膜" / "gse163558_direction.csv")}

    five = [row for row in pairs if int(row["n_sections_above_null"]) >= 5]
    retained = [row for row in pairs if row["retained"] == "True"]
    mac_pm = [float(r["STAT3__mac"]) for r in stat if r["role"] == "peritoneal_tumor"]
    mac_pri = [float(r["STAT3__mac"]) for r in stat if r["role"] == "primary_tumor"]
    soc_pm = [float(r["SOCS3__mac"]) for r in stat if r["role"] == "peritoneal_tumor"]
    soc_pri = [float(r["SOCS3__mac"]) for r in stat if r["role"] == "primary_tumor"]
    epi_pm = [float(r["STAT3__epi"]) for r in stat if r["role"] == "peritoneal_tumor" and int(r["n_epi"]) >= 20]
    epi_pri = [float(r["STAT3__epi"]) for r in stat if r["role"] == "primary_tumor" and int(r["n_epi"]) >= 20]
    st_diff, st_p, st_n = one_sided(mac_pm, mac_pri)
    so_diff, so_p, so_n = one_sided(soc_pm, soc_pri)
    ep_diff, ep_p, ep_n = one_sided(epi_pm, epi_pri)
    below = []
    for row in primary:
        if float(row["median_excess_high"]) < 0:
            below.append(row["pair"])

    doc = Document()
    section = doc.sections[0]
    section.page_width = Cm(21.0)
    section.page_height = Cm(29.7)
    for margin in ("top_margin", "bottom_margin", "left_margin", "right_margin"):
        setattr(section, margin, Cm(2.5))
    sect = section._sectPr
    for child in list(sect):
        if child.tag == qn("w:bidi"):
            sect.remove(child)
    normal = doc.styles["Normal"]
    normal.font.name = FONT
    normal.font.size = Pt(11)
    core = doc.core_properties
    core.title = "C3 ligand in one peritoneal gastric single-cell cohort"
    core.author = ""
    core.subject = "Article"

    block(doc, "Article", size=11, bold=True, center=True, after=8)
    block(doc, "C3 Ligand Is Higher in One Peritoneal Single-Cell Cohort, without a Replicated Ligand–Receptor Excess", size=18, bold=True, center=True, after=8)
    block(doc, "Author names, affiliations, and correspondence are to be completed by the authors.", size=11, center=True, after=10)

    heading(doc, "Abstract")
    block(doc, "This is an exploratory reanalysis of public matrices, not a confirmatory mechanism study. On nine primary Visium sections the screening universe was 1000 single-gene pairs. Five pairs were named before the summaries were read and were judged by the same retention rule. None of the 1000 was retained. A weighted score still had a negative median excess at 150 μm. A second question was written only after class means from GSE183904 had been seen: whether seven class-mean products were higher in three peritoneal tumors than in primary tumors. Four of those seven had q = 0.027 inside the list. Both C3 products were among the four, and the increase was carried by C3, not by macrophage C3AR1 (Spearman 0.01 within 25 primary tumors). That q value does not correct for choosing the list after seeing the data. A class-mean product is a pair of expression averages, not evidence that a ligand binds a receptor or that an axis is active. GSE308231 did not significantly replicate the products. With three versus three samples the smallest one-sided p is 0.05, so p = 0.10 means the test is underpowered, not that the products are absent. The spatial result and the peritoneal single-cell result are separate questions.", indent=True)
    block(doc, "Keywords: gastric cancer; peritoneal metastasis; spatial transcriptomics; Visium; C3; C3AR1; SPP1", size=11, after=6)

    heading(doc, "1. Introduction")
    block(doc, "Gastric cancer remains a leading cause of cancer death worldwide [1]. Peritoneal metastasis is one of the routes that narrows surgery and systemic treatment [2,3]. Model systems and trial designs for this pattern of spread are still incomplete [4,5], and the cellular immunology has been assembled largely from bulk tissue and from single-cell suspensions [6].", indent=True)
    block(doc, "Primary gastric adenocarcinoma is already heterogeneous at the molecular level [7]. The cell types invoked in a stromal hypothesis sit inside a wider tumor immune microenvironment [8]. Cancer-associated fibroblasts are not one state [9], and gastric tumors contain fibroblast subsets linked to invasion and treatment resistance [10,11]. Tumor-associated macrophages are established targets in other cancers [12]. Those reviews describe compartments. They do not measure the products tested below.", indent=True)
    block(doc, "Complement can promote or restrain tumor growth depending on the tissue and the receptor [13]. SPP1-high macrophages form one pole of a program described across human cancers [14]. An SPP1–CD44 interaction between macrophages and tumor cells has been reported in triple-negative breast cancer [15], and a collagen–CD44 axis has been described in a poor-prognosis gastric subtype [16]. TGF-β-linked neutrophil traps have been tied to gastric peritoneal metastasis [17]. STAT3 has many functions in cancer [18] and is an active node in gastric cancer [19]. A separate gastric study connected MELK to macrophage polarization through STAT3 [20]. CCL2 has been implicated in a gastric premetastatic niche in the liver [21]. A Visium study of primary gastric cancer reported CCL2-positive fibroblasts next to STAT3-activated macrophages [22]. That published spatial result is background. It is not re-derived here, and the prespecified CCL2–CCR2 score in the present screen was not retained.", indent=True)
    block(doc, "Spatial transcriptomics places gene expression back onto a tissue section [23,24]. Platforms differ in resolution, and a Visium spot is not a cell [25,26]. Mixture models such as cell2location can assign finer cell types inside a spot [27]. That model was not run. Spots were summarized with marker scores and with non-negative least squares of four programs.", indent=True)
    block(doc, "Public single-cell cohorts already separate primary gastric cancer from metastatic sites, but not in the same way. GSE183904 contains many primary tumors and three peritoneal tumors [28]. GSE163558 profiles organ-specific metastases and contributes one peritoneal sample to the contrast used here [29]. GSE228598 profiles macrophages from tumor and from peritoneal fluid [30]. Ascites single-cell profiles record therapy-associated change and are not solid implants [31]. An immunosuppressive peritoneal microenvironment has been described in another surgical series [32], and a further single-cell study followed differentiation and drug-resistance programs in peritoneal metastasis [33]. GSE308231 is a public count matrix without an accompanying journal article [34]. No public series supplies eight or more paired primary–peritoneal Visium sections.", indent=True)
    block(doc, "Two questions were kept separate. The first was locked before the spatial summaries were interpreted: on nine primary Visium sections, does a fibroblast–macrophage ligand–receptor score exceed a coordinate null often enough to be retained? It does not. The second question was specified after class means from GSE183904 had been seen: are seven class-mean products higher in peritoneal tumor than in primary tumor? Multiplicity for that second question was handled only inside the list of seven, by the Benjamini–Hochberg procedure [35]. Count matrices were processed in Scanpy [36]. The spatial pair list was the single-gene subset of the CellChat human table [37].", indent=True)

    heading(doc, "2. Results")
    heading(doc, "2.1. Primary Visium Sections Do Not Retain a Fibroblast–Macrophage Excess", size=12)
    block(doc, "Table 1 lists the cohorts and the contrast each one can support. These contrasts are not one mechanism. The Visium test uses primary sections only. The peritoneal statements use single-cell class means. Nine primary sections from GSE251950 entered the spatial ranking. GC6-PM (GSM7990479) is named as a metastasis in the series, but the GEO tissue field does not name the peritoneum, and the published Visium analysis of this series used the nine primary sections [22]. It was scored with the same locked rules and was not added to the nine-section median. On that section the C3–C3AR1 weighted score was 0.922 against a null 95th percentile of 0.925, so the section stayed out of the nine-section median. No marker-set or normalization sweep was run after the rules were locked.", indent=True)
    block(doc, "Table 1. Cohorts and the contrast used in this study.", size=10, bold=True, after=2)
    add_table(doc, ["Accession", "Material used here", "Contrast", "What it can support"], [
        ["GSE251950", "9 primary Visium sections; GC6-PM scored apart", "Section-level spatial null", "Primary interface only"],
        ["GSE183904", "3 peritoneal tumors, 26 primary tumors", "Label permutation of samples", "Patient-level products"],
        ["GSE308231", "3 peritoneal and 3 primary tumors", "Separate 3-versus-3 permutation", "Replication, not pooled"],
        ["GSE163558", "1 peritoneal sample, 3 primary tumors", "Direction only", "No p value"],
        ["GSE228598", "28 ascites or washing samples", "Description", "Not a solid-implant test"],
    ])
    block(doc, "Macrophage-high spots touched a fibroblast-high neighbor in 0.81–0.94 of cases. The label-shuffle 95th percentile was 0.95–0.97 on every section, including GC6-PM (Figure 1A). Observed contact sat below the null. The formal screen was the full list of 1000 single-gene pairs, not a second test limited to five pairs. The five prespecified pairs sit inside those 1000 and used the same 5-of-9 and leave-one-out bar. None of the 1000 was retained. Three pairs reached 5 of 9 and then failed leave-one-out: " + ", ".join(row["pair"] for row in five) + f". Retained pairs: {len(retained)}. CCL2–CCR2, C3–C3AR1, SPP1–CD44, TGFB1–TGFBR1, and TGFB1–TGFBR2 did not clear the bar (Figure 1B).", indent=True)

    heading(doc, "2.2. A Weighted Score Removes Zeros and Does Not Create a Positive Excess", size=12)
    block(doc, "The first section summary was a median. Many spots contributed a structural zero, so several prespecified medians were zero even when a minority of spots carried signal. The weighted readout replaces that median with a mean weighted by the product of fibroblast and macrophage proportions. At 150 μm the median excess, defined as the section score minus the 95th percentile of a ligand-shift null, was below zero for all nine pairs (Figure 1C, Table 2). Bootstrap intervals of that median lay entirely below zero for " + ", ".join(below) + ". APP–CD74 had the largest absolute score and a wide interval that crossed zero. Changing the radius to 100, 200, or 300 μm did not produce a positive median excess across the panel (Figure 1D). The single exception of a positive point estimate was TGFB1–TGFBR1 at 100 μm, on the order of 0.0001, and it was not selected after the fact. C3–C3AR1 had a median zero-weight fraction of about three quarters, and CCL2–CCR2 was higher still, so those weighted scores are carried by a minority of spots. Table 2 contains all nine pairs in this weighted panel.", indent=True)
    add_figure(doc, FIG / "fig1_spatial.png", "Figure 1. Primary Visium sections. (A) Contact between a macrophage-high spot and a fibroblast-high neighbor. Blue, observed. Gray, 95th percentile of 200 label shuffles. GC6-PM is shown and is not counted as peritoneum. (B) Primary sections, out of 9, whose ligand–receptor score exceeded the coordinate null. The dashed line is 5. (C) Median excess at 150 μm. The point is the nine-section median and the bar is its bootstrap interval. The dashed line is 0. (D) The same median excess at 100, 150, 200, and 300 μm. Color is centered at 0. Blue is a negative excess.")
    block(doc, "Table 2. Weighted fibroblast–macrophage score at 150 μm on nine primary sections. Excess is the section score minus the ligand-shift null 95th percentile. The interval is the bootstrap interval of the nine-section median. The last column is the weighted-score permutation, not the retention count in Figure 1B.", size=10, bold=True, after=2, page_break=True)
    order = ["CCL2–CCR2", "C3–C3AR1", "SPP1–CD44", "TGFB1–TGFBR1", "TGFB1–TGFBR2", "TGFB1–TGFBR_MIN", "APP–CD74", "CD99–CD99", "SEMA4D–PLXNB2"]
    by_pair = {row["pair"]: row for row in primary}
    add_table(doc, ["Pair", "Median score", "Median excess", "Bootstrap interval", "Sections with p ≤ 0.05"], [
        [name.replace("TGFB1–TGFBR_MIN", "TGFB1–TGFBR min"), f3(by_pair[name]["median_score"]), f3(by_pair[name]["median_excess"]), f"{f3(by_pair[name]['median_excess_low'])} to {f3(by_pair[name]['median_excess_high'])}", f"{int(float(by_pair[name]['n_sections_p_le_0.05']))}/9"]
        for name in order
    ], widths=[3.4, 2.4, 2.6, 3.8, 3.4])

    heading(doc, "2.3. In GSE183904 the C3 Increase Is Mostly the Ligand", size=12)
    c3_fib_pm = class_median(source, "C3", "fib", "pm_median")
    c3_fib_pri = class_median(source, "C3", "fib", "primary_median")
    c3_epi_pm = class_median(source, "C3", "epi", "pm_median")
    c3_epi_pri = class_median(source, "C3", "epi", "primary_median")
    c3ar_pm = class_median(source, "C3AR1", "mac", "pm_median")
    c3ar_pri = class_median(source, "C3AR1", "mac", "primary_median")
    spp_pm = class_median(source, "SPP1", "mac", "pm_median")
    spp_pri = class_median(source, "SPP1", "mac", "primary_median")
    block(doc, f"Cells were labeled by the highest of four marker programs when that program was positive and strictly above the other three. A class with fewer than 20 cells was left missing. Primary sample 22 fell below that count for fibroblasts and epithelial cells, so those denominators are 25 rather than 26. Figure 2A shows the log2 ratio of peritoneal to primary class medians. Fibroblast C3 was {c3_fib_pm:.2f} versus {c3_fib_pri:.2f}, epithelial C3 was {c3_epi_pm:.2f} versus {c3_epi_pri:.2f}, and macrophage C3AR1 was {c3ar_pm:.2f} versus {c3ar_pri:.2f}. Macrophage SPP1 was {spp_pm:.2f} versus {spp_pri:.2f}. The selected bars in Figure 2B make the same comparison for the classes that enter the products: the C3 ratios are large, and the macrophage C3AR1 ratio is not.", indent=True)
    names = {
        "C3_fib_C3AR1_mac": "C3 fib × C3AR1 mac",
        "C3_epi_C3AR1_mac": "C3 epi × C3AR1 mac",
        "SPP1_mac_CD44_mac": "SPP1 mac × CD44 mac",
        "SPP1_mac_CD44_t": "SPP1 mac × CD44 T",
        "CCL2_epi_CCR2_t": "CCL2 epi × CCR2 T",
        "CCL2_fib_CCR2_mac": "CCL2 fib × CCR2 mac",
        "TGFB1_t_TGFBR2_mac": "TGFB1 T × TGFBR2 mac",
    }
    block(doc, f"The seven products in Figure 2C and Table 3 are exploratory. They were named after the class means had been inspected. Multiplying two class means does not show that the ligand and receptor occupy the same neighborhood, that a communication model was fit, or that the receptor was perturbed. Both C3 products use macrophage C3AR1, so they are not two independent receptor findings. Within 25 primary tumors, fibroblast C3 and macrophage C3AR1 had a Spearman correlation of 0.01, and epithelial C3 and macrophage C3AR1 had a Spearman correlation of −0.07 (Figure 3A,B). The receptor does not rise with the ligand inside primary tumors. One donor is paired: sample 38 minus sample 36 was {paired_c3_fib:.3f} for the fibroblast C3 product and {paired_c3_epi:.3f} for the epithelial C3 product. The epithelial CCL2 product differed by {paired_ccl2:.3f} in that same donor. These are differences, not p values.", indent=True)
    add_figure(doc, FIG / "fig2_source.png", "Figure 2. GSE183904 class means. (A) log2 ratio of the peritoneal median to the primary median. Rows are genes and columns are the four programs. Color is centered at 0 and is not clipped. (B) The same ratio for the classes used in the products. (C) Class-mean products. Blue, primary tumors. Orange, peritoneal tumors. The black bar is the primary median.")
    block(doc, f"Table 3. GSE183904 products. The statistic is the peritoneal median minus the primary median. p is a one-sided exact label permutation. q is Benjamini–Hochberg within these seven rows only. Fibroblast and epithelial products use {n_comb_fib} combinations. Macrophage and T-cell products use {n_comb_mac}.", size=10, bold=True, after=2, page_break=True)
    add_table(doc, ["Product", "Peritoneal", "Primary", "Above", "p", "q"], [
        [names[row["pair"]], f3(row["pm_median"]), f3(row["primary_median"]), f"{int(float(row['n_pm_above_primary_median']))}/{int(float(row['pm_n']))}", f3(row["perm_p"]), f3(row["bh_q"])]
        for row in perm
    ], widths=[4.6, 2.3, 2.2, 1.8, 2.2, 2.2])
    block(doc, f"Macrophage STAT3 transcript medians were {median(mac_pm):.3f} and {median(mac_pri):.3f} (difference {st_diff:.3f}, one-sided p = {st_p:.3f}). The p value is the fraction of ways to label 3 of the {len(mac_pm) + len(mac_pri)} macrophage samples as peritoneal in which the median difference is at least as large as observed. That count is C({len(mac_pm) + len(mac_pri)}, 3) = {st_n}. Macrophage SOCS3 was {median(soc_pm):.3f} versus {median(soc_pri):.3f} and was not higher (one-sided p = {so_p:.3f}). Epithelial STAT3 was also higher ({median(epi_pm):.3f} versus {median(epi_pri):.3f}, p = {ep_p:.3f}, n = {len(epi_pm)} peritoneal and {len(epi_pri)} primary samples with at least 20 epithelial cells). STAT3 and SOCS3 were not entered into the seven-row q value, and this p value was not adjusted further. The measurement is a transcript mean, not phosphorylation and not a blockade (Figure 3C,D).", indent=True)
    add_figure(doc, FIG / "fig4_ligand.png", "Figure 3. Ligand, receptor, and two transcripts in GSE183904. (A, B) Primary tumors with at least 20 cells in each class. (C, D) Macrophage transcript means. Blue, primary tumors. Orange, peritoneal tumors. The black bar is the group median.")

    heading(doc, "2.4. The Products Are Not Significantly Higher in a Second Solid-Tumor Series", size=12)
    block(doc, "GSE308231 used the same seven products, the same class rule, and its own permutation. It was not pooled with GSE183904, and it was not treated as a prespecified validation cohort. The platform is SeekOne, the series is not described as paired patients, and three versus three can produce a one-sided p no smaller than 0.05. All seven q values were above 0.35 (Figure 4A,B, Table 4). Failure to cross that floor is not proof that the products are equal. The epithelial C3 product was higher in all three peritoneal samples than in all three primary samples, but one swap of labels preserved the median difference, so p = 0.10. Fibroblast C3 itself was higher in every peritoneal sample than in every primary sample. That ligand pattern was not a second locked endpoint. Macrophage C3AR1 and macrophage STAT3 were not higher. A cross-cohort meta-analysis was not done, because the platforms and the pairing are not the same.", indent=True)
    add_figure(doc, FIG / "fig3_replication.png", "Figure 4. Replication cohorts, kept separate. (A) GSE308231 class-mean products. Blue, primary tumors. Orange, peritoneal tumors. The black bar is the primary median. (B) Peritoneal median minus primary median. Circles, GSE183904. Squares, GSE308231. No interval is drawn. (C) GSE163558. Blue, median of three primary tumors. Orange, the single peritoneal sample. The axis is log(1 + product). No p value was calculated.")
    block(doc, "Table 4. GSE308231, three versus three. q is Benjamini–Hochberg within these seven rows. The permutation floor is 0.05.", size=10, bold=True, after=2, page_break=True)
    add_table(doc, ["Product", "Peritoneal", "Primary", "Above", "p", "q"], [
        [names[row["pair"]], f3(row["pm_median"]), f3(row["primary_median"]), f"{int(float(row['n_pm_above_primary_median']))}/3", f3(row["perm_p"]), f3(row["bh_q"])]
        for row in perm308
    ], widths=[4.6, 2.3, 2.2, 1.8, 2.2, 2.2])
    fib_n = int(float(fluid_map["C3__fib"]["n"]))
    block(doc, f"GSE163558 contributes one peritoneal sample, P1, against three primary tumors. Its epithelial C3 product was {float(g163['C3_epi_C3AR1_mac']['pm_value']):.3f} against a primary median of {float(g163['C3_epi_C3AR1_mac']['primary_median']):.3f} (Figure 4C). No p value was assigned. In GSE228598, fibroblast C3 could be summarized for {fib_n} of 28 fluid samples (median {float(fluid_map['C3__fib']['median']):.2f}). Across all 28, median epithelial C3 was {float(fluid_map['C3__epi']['median']):.2f} and median macrophage C3AR1 was {float(fluid_map['C3AR1__mac']['median']):.2f}. The series has no primary-tumor arm, so these medians were not entered into the solid-tumor permutation. Ascites and a solid implant are different materials [30,31].", indent=True)

    heading(doc, "3. Discussion")
    block(doc, "The spatial result is the one the locked rule can support. Contact was common in absolute terms and still below a label shuffle, because a 150 μm neighborhood on this lattice often reaches a fibroblast-high spot once macrophage-high spots are defined as the top quartile. A large score such as APP–CD74 can remain large after the ligand map is shifted. Absolute abundance is not spatial excess. The weighted mean was added so that a structural zero would not be mistaken for a negative biological result. It removed the zero medians and left the excess negative at the prespecified radius and at the neighboring radii. That is a stable negative readout on nine primary sections. It is not a test of peritoneal tissue.", indent=True)
    block(doc, "The peritoneal products are a different, exploratory question. The Benjamini–Hochberg adjustment was computed only inside the seven rows. It does not undo the fact that the seven rows were chosen after the class means were known, so q = 0.027 is a description of that list, not a confirmatory p value. Four rows reach that q: both C3 products, the epithelial CCL2 product, and the T-cell TGFB1 product (Table 3). The two C3 products share macrophage C3AR1. The fold from primary to peritoneal tumor is much larger for C3 than for C3AR1, and within primary tumors the ligand and the receptor do not move together. The supported statement is that the C3 ligand is higher in this one cohort. It is not that a C3–C3AR1 axis is active. The CCL2 epithelial product is numerically higher and absolutely small. The TGFB1 product uses T cells, not the fibroblast–macrophage pair in the spatial screen. SPP1 × CD44 on macrophages is the comparison suggested by macrophage polarity and by a breast-cancer interaction [14,15], and it is not higher here.", indent=True)
    block(doc, "GSE308231 does not significantly replicate the products, and its size does not let a negative call stand in for a precise null. It does show a higher fibroblast C3 ligand in all three peritoneal samples. Calling that ligand pattern a new significant test would change the endpoint after seeing the split. The epithelial product is higher in every peritoneal sample and still has p = 0.10. GSE163558 is one sample. GSE228598 is fluid. The ascites series records peritoneal fluid [31]. The other peritoneal series [32,33] were not reanalyzed here. None of them is a paired solid implant. A negative primary-interface result cannot be joined to a higher peritoneal C3 ligand and read as one mechanism: the first measurement is spatial and limited to primary sections, and the second is a single-cell class mean. Lee and colleagues have already reported CCL2–STAT3 spatial crosstalk in GSE251950 [22]. The present prespecified CCL2–CCR2 score was not retained, macrophage SOCS3 was not higher, and epithelial STAT3 moved with macrophage STAT3. Those observations do not add phosphorylated STAT3, and they do not make the earlier paper a peritoneal C3 result. NicheNet was not run. cell2location was not run [27].", indent=True)
    block(doc, "What is missing is also specific. A return to a spatial claim about peritoneal implants needs paired primary and peritoneal sections, with more than a handful of patients. Protein-level C3 and C3AR1, and a receptor perturbation, are outside this study. cell2location was not fit [27]. Marker weights are a coarse mixture model, which is why the weighted excess is reported against a coordinate null rather than as a cell-resolution map. Visium spots are on the order of the 55 μm spot and the 100 μm pitch of that assay, not a 20–50 μm cell neighborhood [26].", indent=True)

    heading(doc, "4. Materials and Methods")
    heading(doc, "4.1. Cohorts", size=12)
    block(doc, "GSE251950 supplied the Visium sections. Nine primary sections entered the ranking. GC6-PM was excluded from that ranking. GSE183904 supplied peritoneal tumors GSM5573484, GSM5573485, and GSM5573503 (samples 19, 20, and 38) and 26 primary tumors. Sample 38 and sample 36 are the paired donor. Peritoneal normal sample 37 was excluded. GSE308231 supplied GSM9240670–GSM9240672 as primary tumors and GSM9240673–GSM9240675 as peritoneal tumors. GSE163558 supplied P1 (GSM5004187) against primary tumors PT1–PT3 (GSM5004180–GSM5004182); adjacent mucosa, lymph node, ovary, and liver were not part of the contrast. GSE228598 supplied 28 peritoneal-fluid samples. The GEO record does not split malignant ascites from peritoneal washing, and patient 15 is absent. Cohorts were not pooled. No new patients were enrolled.", indent=True)
    heading(doc, "4.2. Visium Scores", size=12)
    block(doc, "Spots were kept when they were in tissue, detected at least 200 genes, and had a mitochondrial fraction of at most 25%. Counts were normalized to 10,000 and log-transformed with a pseudocount of 1 [36]. Fibroblast markers were DCN, LUM, COL1A1, and FAP. Macrophage markers were CD68, CSF1R, C1QA, and C1QB. Those marker genes were not used as the ligands or receptors under test. High spots were at or above the section’s 75th percentile. An interface spot was a macrophage-high spot with at least one fibroblast-high neighbor inside 150 μm. The spot score was receptor expression times the mean ligand on those neighbors. The section summary was the median. Sections with fewer than 10 interface spots would have been dropped; none were. Contact used 200 label shuffles, seed 20261002, and the 95th percentile of that null. The ligand–receptor null shuffled coordinates 200 times. The section p value was (1 + the number of nulls at least as large as the observed score) / 201. Retention required at least 5 of 9 primary sections above the null, and at least 5 of the remaining 8 after each section was left out. The pair universe was the CellChat human table restricted to single-gene symbols, 1000 pairs [37]. Underscore complexes were excluded. Five pairs were prespecified: CCL2–CCR2, C3–C3AR1, SPP1–CD44, TGFB1–TGFBR1, and TGFB1–TGFBR2.", indent=True)
    heading(doc, "4.3. Weighted Interface", size=12)
    block(doc, "Each spot was decomposed by non-negative least squares into fibroblast, macrophage, epithelial, and T-cell programs, using the markers above plus EPCAM, KRT8, KRT18, and KRT19 for epithelium and CD3D, CD3E, and TRAC for T cells. Coefficients were scaled to sum to 1. The weight was the product of the fibroblast and macrophage proportions. The section score was the weighted mean, not the median. The null translated the ligand map by 200, 300, 400, or 500 μm at a random angle, 200 times, seed 20261004. A shifted coordinate counted only when the nearest real spot after the inverse shift was within 60 μm. Radii of 100, 150, 200, and 300 μm were summarized. The primary radius remained 150 μm. This readout did not replace the retention rule, and it was not used to choose a friendlier radius. cell2location was not run [27].", indent=True)
    heading(doc, "4.4. Class-Mean Products", size=12)
    block(doc, "For the single-cell cohorts, integer counts were library-normalized to 10,000 and log-transformed. Duplicate gene symbols were summed before normalization. A cell took the class with the highest marker mean when that mean was positive and strictly above the second class; otherwise it was unlabeled. Class means with fewer than 20 cells were missing. The seven products were C3 on fibroblasts or epithelial cells times C3AR1 on macrophages, SPP1 on macrophages times CD44 on macrophages or T cells, CCL2 on epithelial cells times CCR2 on T cells, CCL2 on fibroblasts times CCR2 on macrophages, and TGFB1 on T cells times TGFBR2 on macrophages. The statistic was the peritoneal median minus the primary median, one-sided for a higher peritoneal value. q values were computed inside those seven tests only [35]. Products that use a fibroblast or epithelial mean omit primary sample 22, so the permutation count is C(28, 3) = 3276. Products that use macrophage and T-cell means use all 26 primary samples, so the count is C(29, 3) = 3654. The paired donor was reported as sample 38 minus sample 36, which is a difference and not a p value. STAT3 and SOCS3 used the full-transcriptome library size and were not part of the seven-row adjustment. GSE308231 used the same products and a separate permutation. GSE163558 was descriptive. GSE228598 was not compared with a primary arm.", indent=True)

    heading(doc, "5. Conclusions")
    block(doc, "Nine primary Visium sections did not retain a fibroblast–macrophage ligand–receptor excess, and a weighted score did not turn that excess positive. In one single-cell cohort, GSE183904, exploratory C3 products were higher in three peritoneal tumors because the C3 ligand was higher. That list was chosen after the class means were seen, and the within-list q value is not a confirmatory test. A second cohort of three versus three did not significantly replicate the products and was too small to exclude them. The two results stay separate. They do not establish a C3–C3AR1 axis.", indent=True)

    heading(doc, "Author Contributions")
    block(doc, "To be completed by the authors.")
    heading(doc, "Funding")
    block(doc, "To be completed by the authors.")
    heading(doc, "Institutional Review Board Statement")
    block(doc, "Not applicable. The study used de-identified public matrices from the Gene Expression Omnibus. No new patients were enrolled.")
    heading(doc, "Informed Consent Statement")
    block(doc, "Not applicable.")
    heading(doc, "Data Availability Statement")
    block(doc, "GSE251950, GSE183904, GSE308231, GSE163558, and GSE228598 are available from the Gene Expression Omnibus. Raw count matrices were not redeposited. Analysis scripts, seeds, and processed summary tables are in the public repository https://github.com/Cybing521/gastric-peritoneal-spatial-analysis. Scanpy 1.12.3 was used for the matrix objects. The permutation counts and seeds are stated in the Methods.")
    heading(doc, "Acknowledgments")
    block(doc, "To be completed by the authors.")
    heading(doc, "Conflicts of Interest")
    block(doc, "To be completed by the authors.")

    heading(doc, "References")
    for i, (authors, title, journal, year, volume, pages, doi) in enumerate(REFS, start=1):
        if i == 34:
            paragraph = doc.add_paragraph()
            paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT
            fmt = paragraph.paragraph_format
            fmt.left_indent = Cm(0.75)
            fmt.first_line_indent = Cm(-0.75)
            fmt.space_after = Pt(2)
            add_runs(paragraph, [("34. Gene Expression Omnibus. GSE308231. Available online: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308231 (accessed on 4 October 2026).", False, False)], size=9)
            continue
        add_ref(doc, i, author_line(authors), title, journal, year, volume, pages, doi)

    assert not retained
    assert len(five) == 3
    assert "C3–C3AR1" in below and "SPP1–CD44" in below
    assert abs(st_p - 0.010) < 0.002
    doc.save(OUT)
    print(OUT)
    print("stat3", round(st_p, 4), "socs3", round(so_p, 4), "epi", round(ep_p, 4), "below", below)


if __name__ == "__main__":
    main()
