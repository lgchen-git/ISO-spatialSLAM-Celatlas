#!/usr/bin/env python3
"""生成样本前缀、自包含、配置驱动的中文长读长空间 SLAM 报告。"""

from __future__ import annotations

import argparse
import base64
import html
import json
from pathlib import Path

import pandas as pd


def read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def read_key_value_tsv(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    try:
        frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    except (OSError, pd.errors.ParserError):
        return {}
    if not {"field", "value"}.issubset(frame.columns):
        return {}
    return dict(zip(frame["field"].astype(str), frame["value"].astype(str)))


def image_uri(path: Path) -> str:
    if not path.is_file() or path.stat().st_size == 0:
        return ""
    mime = "image/svg+xml" if path.suffix.lower() == ".svg" else "image/png"
    return f"data:{mime};base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def display_title(sample: str, title: str) -> str:
    clean = title.strip()
    return clean if clean.lower().startswith(sample.lower()) else f"{sample} · {clean}"


def figure_card(
    path: Path,
    title: str,
    note: str = "",
    layout: str = "standard",
) -> str:
    uri = image_uri(path)
    if not uri:
        return f'<div class="empty">未生成：{html.escape(title)}</div>'
    note_html = f"<span>{html.escape(note)}</span>" if note else ""
    return (
        f'<figure class="figure-card {html.escape(layout)}">'
        f'<button class="expand" type="button" data-title="{html.escape(title)}">'
        f'<img src="{uri}" alt="{html.escape(title)}"></button>'
        f"<figcaption><strong>{html.escape(title)}</strong>{note_html}</figcaption></figure>"
    )


def fmt(value, mode: str = "number") -> str:
    if value is None or value == "":
        return "NA"
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return html.escape(str(value))
    if mode == "percent":
        return f"{numeric:.2%}"
    if mode == "float":
        return f"{numeric:.3f}"
    return f"{int(numeric):,}" if numeric.is_integer() else f"{numeric:,.2f}"


def cards(items: list[tuple[str, object, str]]) -> str:
    return '<div class="metric-grid">' + "".join(
        f'<div class="metric"><span>{html.escape(label)}</span>'
        f"<strong>{fmt(value)}</strong><small>{html.escape(note)}</small></div>"
        for label, value, note in items
    ) + "</div>"


def table_html(
    path: Path,
    title: str,
    columns: dict[str, str],
    frame: pd.DataFrame | None = None,
    top: int = 50,
) -> str:
    if frame is None:
        if not path.is_file() or path.stat().st_size == 0:
            return f'<div class="empty">未生成：{html.escape(title)}</div>'
        frame = pd.read_csv(path, sep="\t")
    available = [column for column in columns if column in frame.columns]
    if not available or frame.empty:
        return f'<div class="empty">无可展示记录：{html.escape(title)}</div>'
    display = frame.loc[:, available].head(top).copy()
    display.columns = [columns[column] for column in available]
    body = []
    for row in display.itertuples(index=False, name=None):
        cells = []
        for value in row:
            if isinstance(value, float):
                text = f"{value:.4g}" if pd.notna(value) else "NA"
            else:
                text = str(value)
            cells.append(f"<td>{html.escape(text)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    table_id = "table-" + str(abs(hash((str(path), title))))
    head = "".join(f"<th>{html.escape(column)}</th>" for column in display.columns)
    return (
        f'<details class="table-details"><summary>{html.escape(title)}'
        f"（Top {len(display)}）</summary><div class=\"table-card\">"
        f'<div class="table-scroll"><table id="{table_id}" class="data-table" data-page-size="10">'
        f"<thead><tr>{head}</tr></thead><tbody>{''.join(body)}</tbody></table></div>"
        f'<div class="pager" data-table="{table_id}"><button type="button" data-prev>上一页</button>'
        f"<span></span><button type=\"button\" data-next>下一页</button></div></div></details>"
    )


def filtered_markers(path: Path) -> pd.DataFrame:
    if not path.is_file():
        return pd.DataFrame()
    frame = pd.read_csv(path, sep="\t")
    if "passes_display_filter" in frame:
        frame = frame[frame["passes_display_filter"].astype(str).str.lower().eq("true")]
    return frame.sort_values(["cluster", "rank"]).groupby("cluster", observed=True).head(5)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--chip-id", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--chemistry", required=True)
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--report-title", default="长读长空间转录组与 SLAM 分析报告")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()
    root = Path(args.sample_dir)
    outdir = Path(args.outdir)
    meta = outdir / "05_metadata"
    outdir.mkdir(parents=True, exist_ok=True)
    meta.mkdir(parents=True, exist_ok=True)

    read_qc = dict(
        zip(
            *(pd.read_csv(root / "01_read_qc/01_tables/qc_metrics.tsv", sep="\t").values.T)
        )
    )
    structure = read_json(root / "02_barcode_cdna/05_metadata/structure_summary.json")
    barcode = read_json(root / "02_barcode_cdna/05_metadata/barcode_summary.json")
    extract = read_json(root / "02_barcode_cdna/05_metadata/extract_summary.json")
    gate = read_json(root / "03_align/05_metadata/quality_gate.json")
    isoform = read_json(root / "04_full_length_transcriptome/05_metadata/isoform_class_summary.json")
    slam = read_json(root / "05_slam/05_metadata" / f"{args.sample_id}.slam_summary.json")
    gene_slam = read_json(root / "05_slam/05_metadata/gene_slam_summary.json")
    spatial = read_json(root / "06_whole_transcriptome/05_metadata/long_read_spatial_clustering.json")
    spatial_qc = read_json(root / "06_whole_transcriptome/05_metadata/clustering_qc.json")
    joint = read_json(root / "07_slam_joint/05_metadata/slam_joint_summary.json")
    sample_design = read_key_value_tsv(root / "01_read_qc/05_metadata/sample_design.tsv")
    spatial_summary = read_json(root / "06_whole_transcriptome/05_metadata/spatial_qc_summary.json")
    counts = slam.get("class_counts", {})

    read_fig = root / "01_read_qc/02_figures"
    barcode_fig = root / "02_barcode_cdna/02_figures"
    align_fig = root / "03_align/02_figures"
    isoform_fig = root / "04_full_length_transcriptome/02_figures"
    slam_fig = root / "05_slam/02_figures"
    spatial_fig = root / "06_whole_transcriptome/02_figures"
    joint_fig = root / "07_slam_joint/02_figures"
    marker_path = root / "06_whole_transcriptome/01_tables/long_read_cluster_markers.tsv"
    tissue_image = Path(sample_design.get("tissue_image", "")) if sample_design.get("tissue_image", "").strip() else Path()
    design_table = table_html(
        root / "01_read_qc/05_metadata/sample_design.tsv",
        "实验设计元数据",
        {"field": "字段", "value": "值"},
    )

    nav = [
        ("s1", "1 样本与输入"),
        ("s2", "2 原始数据质控"),
        ("s3", "3 条形码与 cDNA"),
        ("s4", "4 比对结果"),
        ("s5", "5 全长转录组"),
        ("s6", "6 SLAM 分析"),
        ("s7", "7 长读长空间转录组"),
        ("s8", "8 长读长空间转录组与 SLAM 联合"),
        ("s9", "9 方法与参数"),
    ]
    nav_html = "".join(f'<a href="#{anchor}">{label}</a>' for anchor, label in nav)
    title = display_title(args.sample_id, args.report_title)
    css = """
    :root{--navy:#253746;--cyan:#4DBBD5;--red:#E64B35;--teal:#00A087;--ink:#27323a;--muted:#667681;--line:#dce3e7;--soft:#f5f8f9}
    *{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:#f2f5f6;color:var(--ink);font-family:Arial,"Noto Sans SC","Microsoft YaHei",sans-serif;font-size:13px;line-height:1.55}
    .sidebar{position:fixed;inset:0 auto 0 0;width:228px;background:var(--navy);color:#fff;padding:25px 18px;overflow:auto;z-index:20}.brand{font-size:11px;font-weight:700;letter-spacing:.12em;color:#91D1C2}.sidebar h1{font-size:17px;line-height:1.35;margin:8px 0}.sidebar p{color:#c8d2d8;margin:0 0 16px}.sidebar a{display:block;color:#dce7eb;text-decoration:none;padding:7px 9px;margin:2px 0;border-left:2px solid transparent}.sidebar a.active,.sidebar a:hover{background:#344b5d;border-left-color:var(--cyan);color:#fff}
    .main{margin-left:228px;width:min(1380px,calc(100vw - 228px));padding:26px 32px 58px}.title-card,section{background:#fff;border:1px solid var(--line);border-radius:7px}.title-card{padding:23px 26px 21px;border-top:4px solid var(--navy)}.eyebrow{font-size:10px;font-weight:700;letter-spacing:.08em;text-transform:uppercase;color:var(--red)}.title-card h1{font-size:25px;line-height:1.3;margin:4px 0 7px}.title-card p{margin:3px 0;color:var(--muted)}.chips{display:flex;gap:6px;flex-wrap:wrap;margin-top:12px}.chip{border:1px solid var(--line);padding:3px 9px;border-radius:999px;font-size:11px;background:var(--soft)}
    section{margin-top:16px;padding:21px 23px 24px;scroll-margin-top:12px}h2{display:flex;align-items:center;gap:9px;margin:0 0 7px;font-size:19px;color:var(--navy)}h2:after{content:"";height:1px;background:var(--line);flex:1}h3{margin:20px 0 7px;font-size:14px;color:#354b5c}.method{margin:12px 0;padding:11px 13px;background:var(--soft);border-left:3px solid var(--cyan)}.empty{margin:12px 0;padding:10px;border:1px dashed var(--line);color:var(--muted)}
    .metric-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:9px;margin:12px 0}.metric{border-top:2px solid var(--cyan);background:var(--soft);padding:10px 11px;min-width:0}.metric span,.metric small{display:block;color:var(--muted);font-size:10px}.metric strong{display:block;margin:3px 0;font-size:18px;color:var(--navy);overflow-wrap:anywhere}
    .figure-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;margin:12px 0 16px;align-items:start}.figure-card{counter-increment:figure;margin:0;border:1px solid var(--line);min-width:0;background:#fff}.figure-card.wide,.figure-card.panorama{grid-column:1/-1}.expand{display:block;width:100%;border:0;background:#fff;padding:10px;cursor:zoom-in}.expand img{display:block;width:100%;height:auto;max-height:620px;object-fit:contain;margin:auto}.figure-card.panorama .expand img{max-height:520px}.figure-card figcaption{display:flex;flex-direction:column;min-height:49px;padding:8px 10px;border-top:1px solid var(--line);background:#fafbfb}.figure-card figcaption strong:before{content:"图 " counter(figure) " · "}.figure-card figcaption strong{font-size:11px;color:#354b5c}.figure-card figcaption span{margin-top:2px;color:var(--muted);font-size:10px}
    .table-details,details{margin-top:11px;border:1px solid var(--line);background:#fafbfb}.table-details>summary,details>summary{cursor:pointer;padding:9px 11px;font-weight:700;color:#455d70}.table-card{border-top:1px solid var(--line);padding:0 10px 9px;background:#fff}.table-scroll{overflow-x:auto;max-height:520px}.data-table{width:100%;border-collapse:collapse;font-size:11px}.data-table th,.data-table td{padding:7px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:normal;overflow-wrap:anywhere}.data-table th{position:sticky;top:0;background:var(--soft);z-index:1}.data-table tbody tr:nth-child(even){background:#fafbfb}.pager{display:flex;justify-content:flex-end;align-items:center;gap:7px;margin-top:8px;color:var(--muted);font-size:10px}.pager button{border:1px solid var(--line);background:#fff;padding:3px 8px}.pager button:disabled{opacity:.4}
    .footer{margin-top:20px;padding:16px 4px;color:var(--muted);font-size:10px}.modal{display:none;position:fixed;inset:0;background:rgba(10,20,28,.88);z-index:100;padding:24px}.modal.open{display:flex;align-items:center;justify-content:center}.modal-box{position:relative;max-width:96vw;max-height:94vh;background:#fff;padding:40px 14px 14px}.modal-box img{display:block;max-width:92vw;max-height:82vh}.modal-title{position:absolute;left:15px;top:10px;font-weight:700}.modal-close{position:absolute;right:10px;top:5px;border:0;background:none;font-size:27px;cursor:pointer}
    @media(max-width:1120px){.main{padding:22px}.figure-grid{grid-template-columns:1fr}.figure-card.wide,.figure-card.panorama{grid-column:auto}.metric-grid{grid-template-columns:repeat(2,minmax(0,1fr))}}@media(max-width:760px){.sidebar{position:static;width:auto}.sidebar nav{display:flex;gap:5px;overflow-x:auto}.sidebar a{flex:0 0 auto;border-left:0;border-bottom:2px solid transparent}.main{margin:0;width:100%;padding:10px}.title-card,section{padding:17px}}
    """
    js = """
    document.querySelectorAll('table[data-page-size]').forEach(table=>{const rows=[...table.tBodies[0].rows],size=+table.dataset.pageSize,pager=document.querySelector(`[data-table="${table.id}"]`),label=pager.querySelector('span');let page=0;function draw(){const pages=Math.max(1,Math.ceil(rows.length/size));page=Math.min(page,pages-1);rows.forEach((row,i)=>row.hidden=i<page*size||i>=(page+1)*size);label.textContent=`${page+1} / ${pages}`;pager.querySelector('[data-prev]').disabled=page===0;pager.querySelector('[data-next]').disabled=page===pages-1}pager.querySelector('[data-prev]').onclick=()=>{page--;draw()};pager.querySelector('[data-next]').onclick=()=>{page++;draw()};draw()});
    const modal=document.getElementById('figure-modal'),img=modal.querySelector('img'),mt=modal.querySelector('.modal-title');document.querySelectorAll('.expand').forEach(btn=>btn.onclick=()=>{const source=btn.querySelector('img');img.src=source.src;img.alt=source.alt;mt.textContent=btn.dataset.title;modal.classList.add('open')});function close(){modal.classList.remove('open');img.src=''}modal.querySelector('.modal-close').onclick=close;modal.onclick=e=>{if(e.target===modal)close()};document.onkeydown=e=>{if(e.key==='Escape')close()};
    const links=[...document.querySelectorAll('.sidebar nav a')];if('IntersectionObserver' in window){const observer=new IntersectionObserver(entries=>{entries.filter(x=>x.isIntersecting).forEach(entry=>links.forEach(link=>link.classList.toggle('active',link.getAttribute('href')==='#'+entry.target.id)))},{rootMargin:'-15% 0px -75% 0px'});document.querySelectorAll('main section').forEach(section=>observer.observe(section))}
    """

    isoform_table = table_html(
        root / "04_full_length_transcriptome/01_tables/isoform_class_summary.tsv",
        "转录本结构类别",
        {"class_code": "类别", "description": "说明", "transcripts": "转录本数", "fraction": "比例"},
    )
    mutation_table = table_html(
        root / "05_slam/01_tables/base_substitution_comparison.tsv",
        "12 种单碱基替换统计",
        {"channel": "替换类型", "count": "错配数", "reference_coverage": "参考覆盖", "rate_percent": "错配率（%）", "expected_slam_channel": "SLAM 通道"},
    )
    gene_path = root / "05_slam/01_tables/gene_new_supported_summary.tsv"
    gene_frame = pd.read_csv(gene_path, sep="\t") if gene_path.is_file() else pd.DataFrame()
    if not gene_frame.empty and "passes_min_informative" in gene_frame:
        gene_frame = gene_frame[gene_frame["passes_min_informative"].astype(str).str.lower().eq("true")]
    gene_table = table_html(
        gene_path,
        "基因层面新生支持比例",
        {"gene_name": "基因", "gene_id": "基因 ID", "informative_molecules": "信息分子数", "new_supported": "新生支持", "new_supported_fraction": "新生支持比例", "wilson95_low": "95% CI 下限", "wilson95_high": "95% CI 上限"},
        frame=gene_frame,
    )
    cluster_table = table_html(
        root / "06_whole_transcriptome/01_tables/long_read_cluster_summary.tsv",
        "长读长空间簇汇总",
        {"cluster": "表达簇", "bins": "Bin 数", "total_lr_molecules": "分子数", "median_lr_molecules": "每 Bin 中位分子数", "median_detected_genes": "每 Bin 中位基因数"},
    )
    marker_table = table_html(
        marker_path,
        "各簇标志基因",
        {"cluster": "表达簇", "rank": "排名", "gene_id": "基因 ID", "gene_name": "基因", "logfoldchange": "log2 倍数变化", "pvalue_adj": "校正 P 值"},
        frame=filtered_markers(marker_path),
    )
    joint_table = table_html(
        root / "07_slam_joint/01_tables/cluster_slam_summary.tsv",
        "表达簇与 SLAM 联合汇总",
        {"cluster": "表达簇", "new_supported": "新生支持", "no_conversion_detected": "未检出转换", "ambiguous": "不明确", "informative_molecules": "信息分子", "new_supported_fraction": "新生支持比例", "wilson95_low": "95% CI 下限", "wilson95_high": "95% CI 上限"},
    )
    joint_gene_path = root / "07_slam_joint/01_tables/cluster_gene_slam.tsv"
    joint_gene = pd.read_csv(joint_gene_path, sep="\t") if joint_gene_path.is_file() else pd.DataFrame()
    if not joint_gene.empty:
        joint_gene = joint_gene[joint_gene["passes_min_informative"].astype(str).str.lower().eq("true")]
    joint_gene_table = table_html(
        joint_gene_path,
        "基因 × 表达簇联合结果",
        {"cluster": "表达簇", "gene_name": "基因", "gene_id": "基因 ID", "informative_molecules": "信息分子", "new_supported": "新生支持", "new_supported_fraction": "新生支持比例", "wilson95_low": "95% CI 下限", "wilson95_high": "95% CI 上限"},
        frame=joint_gene,
    )
    isoform_summary_path = root / "04_full_length_transcriptome/01_tables/isoform_summary.tsv"
    if isoform_summary_path.is_file():
        isoform_summary_frame = pd.read_csv(isoform_summary_path, sep="\t")
        isoform_summary_metrics = dict(zip(isoform_summary_frame["metric"], isoform_summary_frame["value"]))
    else:
        isoform_summary_metrics = {}

    body = f"""<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title><style>{css}</style></head><body>
    <aside class="sidebar"><div class="brand">LR-Celatlas</div><h1>{html.escape(title)}</h1><p>样本 {html.escape(args.sample_id)}</p><nav>{nav_html}</nav></aside>
    <main class="main"><header class="title-card"><div class="eyebrow">Long-read spatial transcriptome + SLAM</div><h1>{html.escape(title)}</h1><p>芯片：{html.escape(args.chip_id)}　·　物种：{html.escape(args.species)}　·　化学版本：{html.escape(args.chemistry)}</p><div class="chips"><span class="chip">Bin {fmt(spatial.get("bin_size_um"))} µm</span><span class="chip">离线自包含</span><span class="chip">配置驱动</span></div></header>

    <section id="s1"><h2>1 样本与输入</h2>{cards([("样本 ID", args.sample_id, "样本表"), ("芯片 ID", args.chip_id, "样本表"), ("物种", args.species, "样本表"), ("化学版本", args.chemistry, "样本表"), ("空间单位", f"{fmt(spatial.get('bin_size_um'))} µm", "手工配置")])}{design_table}{figure_card(tissue_image, "组织图", "来自样本表 tissue_image；未提供时不生成组织学图")}</section>

    <section id="s2"><h2>2 原始数据质控</h2>{cards([("读段数", read_qc.get("n_reads"), "输入读段"), ("碱基数", read_qc.get("bases"), "输入碱基总量"), ("读长 N50", read_qc.get("N50"), "bp"), ("平均质量值", read_qc.get("mean_q"), "平均 Q 值")])}<div class="figure-grid">{figure_card(read_fig/"raw_read_length.png","原始读长分布")}{figure_card(read_fig/"raw_mean_q.png","原始读段质量分布")}</div></section>

    <section id="s3"><h2>3 条形码与 cDNA</h2>{cards([("结构候选读段", structure.get("candidate_found"), "结构扫描"), ("精确匹配条形码", barcode.get("assign_counts",{}).get("exact"), "白名单匹配"), ("校正条形码", barcode.get("assign_counts",{}).get("rescued"), "允许范围内校正"), ("保留 cDNA", extract.get("counts",{}).get("kept"), "进入比对")])}<div class="figure-grid">{figure_card(barcode_fig/"structure_status.png","读段结构识别")}{figure_card(barcode_fig/"cdna_length.png","提取 cDNA 长度分布")}{figure_card(barcode_fig/"barcode_status.png","条形码分配状态")}{figure_card(barcode_fig/"barcode_edit_distance.png","条形码编辑距离")}</div></section>

    <section id="s4"><h2>4 比对结果</h2>{cards([("质量门控", gate.get("overall","NA"), "配置阈值"), ("MAPQ≥20 比例", fmt(gate.get("mapq20_rate"),"percent"), "保留 cDNA 读段"), ("坐标范围", "通过" if gate.get("coordinates_valid") else "未通过", "芯片边界")])}<div class="figure-grid">{figure_card(align_fig/"alignment_summary.png","剪接比对汇总")}{figure_card(align_fig/"top50_barcodes.png","Top 50 空间条形码")}</div></section>

    <section id="s5"><h2>5 全长转录组</h2>{cards([("候选转录本", isoform.get("transcripts"), "FLAMES"), ("精确参考匹配", isoform.get("class_counts",{}).get("="), "gffcompare 类别 ="), ("新剪接组合", isoform.get("class_counts",{}).get("j"), "gffcompare 类别 j"), ("多外显子转录本", isoform_summary_metrics.get("multi_exon_transcripts"), "结构汇总")])}<div class="method">FLAMES 用于全长转录本发现，gffcompare 用于参考注释结构比较；本节报告候选转录本结构，不将 gffcompare 的零值 TPM/FPKM 解释为表达定量。</div><div class="figure-grid">{figure_card(isoform_fig/"isoform_class_distribution.png","候选转录本结构类别")}{figure_card(isoform_fig/"isoform_structure_distribution.png","单/多外显子结构")}</div>{isoform_table}</section>

    <section id="s6"><h2>6 SLAM 分析</h2>{cards([("分子总数", slam.get("molecules_total"), "条形码、UMI 与基因聚合"), ("新生支持", counts.get("new_supported",0), "手工阈值"), ("未检出转换", counts.get("no_conversion_detected",0), "覆盖满足阈值"), ("不明确", counts.get("ambiguous",0), "未满足分类条件"), ("达基因阈值", gene_slam.get("genes_passing_threshold"), "信息分子阈值")])}<h3>碱基替换统计</h3><div class="figure-grid">{figure_card(slam_fig/"base_substitution_comparison.png","12 种单碱基替换率与错配数",layout="wide")}</div>{mutation_table}<h3>分子与基因汇总</h3><div class="figure-grid">{figure_card(slam_fig/"slam_class_composition.png","SLAM 分子类别构成")}{figure_card(slam_fig/"gene_new_supported_fraction.png","基因层面新生支持比例","误差线为 Wilson 95% CI")}</div>{gene_table}</section>

    <section id="s7"><h2>7 长读长空间转录组</h2>{cards([("输入分子", spatial.get("input_molecules"), "全部长读长分子"), ("有效 Bin", spatial.get("clustered_bins"), "全视野空间 Bin"), ("检出基因", spatial.get("genes_after_filter"), "过滤后"), ("Leiden 簇", spatial.get("clusters"), "手工配置"), ("Bin 中位 UMI", spatial_summary.get("median_umi_per_bin"), "bin50"), ("Bin 中位 detected gene", spatial_summary.get("median_detected_genes_per_bin"), "bin50")])}<div class="method">长读长空间转录组矩阵使用全部有效长读长分子，不按 SLAM 类别筛选；聚类参数直接读取配置文件，不在流程中自动扫描或推荐。空间 QC 中的 UMI 定义为每个 Bin 内唯一 barcode+UMI；样本级 isoform 数来自 FLAMES，不代表每个 Bin 的 isoform 数。</div><h3>核心空间结果</h3><div class="figure-grid">{figure_card(spatial_fig/"long_read_spatial_clusters.png","空间表达簇")}{figure_card(spatial_fig/"long_read_umap_clusters.png","空间表达 UMAP")}</div>{cluster_table}{table_html(root / "06_whole_transcriptome/01_tables/spatial_bin_qc.tsv", "bin50 空间 QC", {"bin_id": "Bin", "cluster": "表达簇", "umi_count": "UMI 数", "lr_molecule_count": "长读长分子数", "detected_gene_count": "检出基因数", "new_supported_fraction": "探索性 dynamic proxy"})}{table_html(root / "06_whole_transcriptome/01_tables/spatial_qc_summary.tsv", "空间 QC 汇总", {"metric": "指标", "value": "值"})}<h3>表达深度与标志基因</h3><div class="figure-grid">{figure_card(spatial_fig/"long_read_molecules.png","每 Bin 分子数")}{figure_card(spatial_fig/"cluster_marker_dotplot.png","标志基因点图",layout="wide")}{figure_card(spatial_fig/"cluster_depth_qc.png","各簇测序深度与检出基因数",layout="panorama")}{figure_card(spatial_fig/"spatial_coherence.png","空间连续性",layout="panorama")}</div>{marker_table}</section>

    <section id="s8"><h2>8 长读长空间转录组与 SLAM 联合</h2>{cards([("联合分子", joint.get("molecules_joined_to_whole_transcriptome_bins"), "落入有效空间 Bin"), ("表达簇", joint.get("clusters"), "复用长读长空间表达簇"), ("联合基因", joint.get("genes"), "基因 × 表达簇"), ("信息分子阈值", joint.get("min_informative_molecules"), "手工配置")])}<div class="method">SLAM 结果叠加到同一套长读长空间表达簇上；old/new 不分别重新聚类。当前联合层级为 Bin、表达簇、基因和基因 × 表达簇。没有未标记对照时，new-supported 及其空间比例统一解释为探索性 dynamic proxy，不作为绝对 nascent RNA 结论。</div><div class="figure-grid">{figure_card(joint_fig/"cluster_slam_composition.png","各表达簇 SLAM 类别构成")}{figure_card(joint_fig/"cluster_new_supported_fraction.png","各表达簇探索性 dynamic proxy")}{figure_card(joint_fig/"cluster_gene_slam_dotplot.png","基因表达支持与探索性 dynamic proxy",layout="wide")}</div>{joint_table}{joint_gene_table}</section>

    <section id="s9"><h2>9 方法与参数</h2><div class="method">全长转录本由 FLAMES 发现并用 gffcompare 与参考注释比较。空间矩阵按配置的 Bin 大小聚合，依次执行总量归一化、log1p、HVG、标准化、PCA、深度线性校正、表达图与空间四邻域图加权、Leiden 聚类、小簇合并和 UMAP。SLAM 分类阈值、聚类分辨率、空间图权重、Marker 阈值及随机种子均由配置文件手工给定，并记录在阶段元数据中。</div><details><summary>统计定义与范围</summary><p>Marker 使用 t-test_overestim_var 进行簇对其余 Bin 的描述性排序；联合分析的新生支持比例使用信息分子为分母，并给出 Wilson 95% 置信区间。当前分子表不含 read→isoform 映射，因此不输出转录本级 SLAM。</p></details></section>

    <footer class="footer">LR-Celatlas 长读长空间转录组与 SLAM · 离线自包含报告</footer></main>
    <div id="figure-modal" class="modal" role="dialog" aria-modal="true"><div class="modal-box"><span class="modal-title"></span><button class="modal-close" type="button">×</button><img alt=""></div></div><script>{js}</script></body></html>"""

    report = (outdir / f"{args.sample_id}_report.html").resolve()
    legacy = outdir / "report.html"
    if legacy.exists() and legacy.resolve() != report:
        legacy.unlink()
    report.write_text(body, encoding="utf-8")
    (meta / "delivery_manifest.tsv").write_text(
        "description\tsource\ttarget\tstatus\n"
        f"{args.sample_id} 最终离线中文 HTML 报告\t{report}\t{report}\tup_to_date\n",
        encoding="utf-8",
    )
    print(report)


if __name__ == "__main__":
    main()
