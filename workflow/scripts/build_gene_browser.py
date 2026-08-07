#!/usr/bin/env python3
"""Self-contained interactive gene browser.

One HTML file with an embedded JSON payload: search a gene or fragmentform id, see every
fragmentform's exon structure drawn to scale, click an exon (or drag-select a region) and the site
table filters to the modification sites inside it -- per-fragmentform stoichiometry, differential
results, APA motif class and poly(A) tail length.

The payload is built per gene so the browser stays responsive: genes are indexed by name and only
the selected gene's detail is rendered. No frameworks, no CDN -- it opens anywhere, offline.
"""
import argparse
import json
import os
import re

import numpy as np
import pandas as pd


def parse_args():
    ap = argparse.ArgumentParser(description="Build the interactive gene/fragmentform browser HTML.")
    ap.add_argument("--gtf", required=True, help="Assembled fragmentform GTF (exon structures)")
    ap.add_argument("--sites-long", default="", help="*_FILTERED_sites_long.tsv (per-site x transcript x sample)")
    ap.add_argument("--diff-results", default="", help="*__ZN_site_diff_results.tsv")
    ap.add_argument("--classification-summary", default="", help="*_classification_summary.tsv")
    ap.add_argument("--apa-motifs", default="", help="*_apa_motifs.tsv")
    ap.add_argument("--polya-fragmentform", default="", help="*_polya_fragmentform.tsv")
    ap.add_argument("--condition-mod-diffs", default="", help="between_conditions *_mod_diffs.tsv")
    ap.add_argument("--hierarchical-stoich", default="", help="*_hierarchical_stoich.tsv")
    ap.add_argument("--out-html", required=True)
    ap.add_argument("--max-genes", type=int, default=4000, help="Cap genes embedded (largest by read support)")
    ap.add_argument("--title", default="modulator gene browser")
    ap.add_argument("--verbose", action="store_true")
    return ap.parse_args()


def _read(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df.columns = [str(c).lstrip("#") for c in df.columns]
    return df


def load_structures(gtf_path):
    """gene -> [{zt, chrom, strand, exons:[[s,e]..]}]  (transcript_id -> gene via the GTF attrs)."""
    tx = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            g = re.search(r'gene_name "([^"]+)"', f[8]) or re.search(r'gene_id "([^"]+)"', f[8])
            t = tx.setdefault(m.group(1), {"zt": m.group(1), "chrom": f[0], "strand": f[6],
                                           "gene": g.group(1) if g else m.group(1).split(".")[0],
                                           "exons": []})
            t["exons"].append([int(f[3]) - 1, int(f[4])])
    genes = {}
    for t in tx.values():
        t["exons"].sort()
        genes.setdefault(t["gene"], []).append(t)
    return genes


def main():
    args = parse_args()
    genes = load_structures(args.gtf)
    if args.verbose:
        print(f"[browser] {sum(len(v) for v in genes.values()):,} fragmentforms in {len(genes):,} genes", flush=True)

    sites = _read(args.sites_long)
    diffs = _read(args.diff_results)
    summ = _read(args.classification_summary)
    apa = _read(args.apa_motifs)
    tails = _read(args.polya_fragmentform)
    cond = _read(args.condition_mod_diffs)
    hier = _read(args.hierarchical_stoich)

    # per-fragmentform annotations, keyed by the zt_label used in the GTF (gene name stripped)
    def _gtfkey(zt):
        p = str(zt).split(".")
        return ".".join(p[1:]) if len(p) > 3 else str(zt)

    ff_ann = {}
    if not summ.empty and "zt_label" in summ.columns:
        for r in summ.itertuples(index=False):
            ff_ann[_gtfkey(r.zt_label)] = {
                "classification": getattr(r, "classification", ""),
                "reads": int(getattr(r, "read_support", 0) or 0),
                "tes": int(getattr(r, "iso_tes", 0) or 0),
            }
    if not apa.empty and "zt_label" in apa.columns:
        for r in apa.itertuples(index=False):
            ff_ann.setdefault(_gtfkey(r.zt_label), {}).update({
                "pas": getattr(r, "apa_motif_class", ""), "pas_motif": getattr(r, "pas_motif", "") or ""})
    if not tails.empty and "ZT" in tails.columns:
        for r in tails.itertuples(index=False):
            ff_ann.setdefault(_gtfkey(r.ZT), {}).update({
                "tail": float(getattr(r, "median_tail", float("nan"))),
                "tail_n": int(getattr(r, "n_reads", 0) or 0)})

    # per-gene site stoichiometry: collapse samples -> per (site, transcript) modified fraction
    site_by_gene = {}
    if not sites.empty:
        need = {"gene_name", "chrom", "start0", "strand", "mod_code", "Nvalid_cov", "Nmod"}
        if need.issubset(sites.columns):
            s = sites.copy()
            keys = ["gene_name", "chrom", "start0", "strand", "mod_code"]
            if "ZN_transcript_index" in s.columns:
                keys_zt = keys + ["ZN_transcript_index"]
                g = s.groupby(keys_zt, sort=False)[["Nvalid_cov", "Nmod"]].sum().reset_index()
            else:
                g = s.groupby(keys, sort=False)[["Nvalid_cov", "Nmod"]].sum().reset_index()
                g["ZN_transcript_index"] = -1
            g["frac"] = (g["Nmod"] / g["Nvalid_cov"].replace(0, np.nan)).round(4)
            for gene, gg in g.groupby("gene_name", sort=False):
                site_by_gene[str(gene)] = gg

    diff_by_gene = {str(k): v for k, v in diffs.groupby("gene_name", sort=False)} if "gene_name" in diffs.columns else {}
    cond_by_gene = {str(k): v for k, v in cond.groupby("gene_name", sort=False)} if "gene_name" in cond.columns else {}
    hier_by_gene = {str(k): v for k, v in hier.groupby("gene_name", sort=False)} if "gene_name" in hier.columns else {}

    payload = {}
    for gene, forms in genes.items():
        ann_reads = sum(ff_ann.get(f["zt"], {}).get("reads", 0) for f in forms)
        rec = {
            "gene": gene, "chrom": forms[0]["chrom"], "strand": forms[0]["strand"],
            "reads": ann_reads,
            "forms": [{"zt": f["zt"], "exons": f["exons"], **ff_ann.get(f["zt"], {})} for f in forms],
            "sites": [], "diffs": [], "cond": [], "hier": [],
        }
        sg = site_by_gene.get(gene)
        if sg is not None:
            for r in sg.itertuples(index=False):
                rec["sites"].append({"pos": int(r.start0), "mod": str(r.mod_code),
                                     "zn": int(getattr(r, "ZN_transcript_index", -1)),
                                     "cov": int(r.Nvalid_cov), "frac": (None if pd.isna(r.frac) else float(r.frac))})
        dg = diff_by_gene.get(gene)
        if dg is not None:
            for r in dg.itertuples(index=False):
                rec["diffs"].append({"pos": int(getattr(r, "start0", 0)), "mod": str(getattr(r, "mod_code", "")),
                                     "effect": float(getattr(r, "effect_max_abs_frac_diff", float("nan")) or 0),
                                     "padj": float(getattr(r, "p_adj_bh", float("nan")) or 1)})
        cg = cond_by_gene.get(gene)
        if cg is not None:
            for r in cg.itertuples(index=False):
                rec["cond"].append({"pos": int(getattr(r, "start0", 0)), "mod": str(getattr(r, "mod_code", "")),
                                    "delta": float(getattr(r, "delta", 0) or 0),
                                    "padj": float(getattr(r, "p_adj_bh", 1) or 1),
                                    "contrast": str(getattr(r, "contrast", ""))})
        hg = hier_by_gene.get(gene)
        if hg is not None:
            for r in hg.itertuples(index=False):
                rec["hier"].append({"pos": int(getattr(r, "site_pos", 0)),
                                    "a": _gtfkey(getattr(r, "fragmentform_a", "")),
                                    "b": _gtfkey(getattr(r, "fragmentform_b", "")),
                                    "delta": float(getattr(r, "delta", 0) or 0),
                                    "padj": float(getattr(r, "p_adj_bh", 1) or 1),
                                    "ninf": int(getattr(r, "n_informative", 0) or 0),
                                    "div3p": int(getattr(r, "divergence_from_3p_nt", 0) or 0)})
        payload[gene] = rec

    ranked = sorted(payload.values(), key=lambda r: (-r["reads"], r["gene"]))[:args.max_genes]
    data = {"genes": {r["gene"]: r for r in ranked},
            "index": [{"g": r["gene"], "n": len(r["forms"]), "r": r["reads"], "c": r["chrom"]} for r in ranked]}
    os.makedirs(os.path.dirname(args.out_html) or ".", exist_ok=True)
    with open(args.out_html, "w") as fh:
        fh.write(_HTML.replace("__TITLE__", args.title)
                      .replace("__DATA__", json.dumps(data, separators=(",", ":"))))
    if args.verbose:
        print(f"[browser] wrote {len(ranked):,} genes -> {args.out_html} "
              f"({os.path.getsize(args.out_html)/1e6:.1f} MB)", flush=True)


_HTML = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title><style>
:root{--bg:#fbfcfd;--panel:#fff;--ink:#16202b;--muted:#5d6b7a;--line:#e3e9ef;--accent:#2b6a9c;
--accent-soft:#eaf2f8;--exon:#3b6ea5;--exon-alt:#7fa8cc;--hit:#c1121f;--ok:#2f6a4f;--stripe:#f5f8fa}
@media(prefers-color-scheme:dark){:root{--bg:#0e131a;--panel:#141c25;--ink:#e6edf4;--muted:#93a2b2;
--line:#243040;--accent:#5fa8dd;--accent-soft:#16232f;--exon:#5fa8dd;--exon-alt:#3c6c92;--hit:#e0846f;
--ok:#7fc9a3;--stripe:#111922}}
*{box-sizing:border-box}html,body{height:100%}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}
.app{display:grid;grid-template-columns:280px 1fr;height:100vh}
aside{border-right:1px solid var(--line);background:var(--panel);display:flex;flex-direction:column;min-height:0}
.brand{padding:14px 16px;border-bottom:1px solid var(--line)}
.brand b{font-family:Georgia,serif;font-size:17px;letter-spacing:-.01em}
.brand span{display:block;color:var(--muted);font-size:11.5px;margin-top:2px}
#q{width:100%;padding:9px 11px;border:1px solid var(--line);border-radius:8px;background:var(--bg);
color:var(--ink);font-size:14px;outline:none}#q:focus{border-color:var(--accent)}
.search{padding:12px 14px;border-bottom:1px solid var(--line)}
#list{overflow:auto;flex:1;min-height:0}
.gi{padding:8px 14px;border-bottom:1px solid var(--line);cursor:pointer;display:flex;justify-content:space-between;gap:8px}
.gi:hover{background:var(--accent-soft)}.gi.sel{background:var(--accent-soft);box-shadow:inset 3px 0 var(--accent)}
.gi b{font-weight:600;font-size:13.5px}.gi small{color:var(--muted);font-size:11px;font-variant-numeric:tabular-nums}
main{overflow:auto;padding:20px 24px 60px;min-width:0}
h1{font-family:Georgia,serif;font-size:26px;margin:0 0 2px}
.sub{color:var(--muted);font-size:13px;margin-bottom:16px;font-variant-numeric:tabular-nums}
.card{background:var(--panel);border:1px solid var(--line);border-radius:11px;padding:14px 16px;margin-bottom:16px}
.card h2{font-size:13px;text-transform:uppercase;letter-spacing:.08em;color:var(--accent);margin:0 0 10px;font-weight:650}
.hint{color:var(--muted);font-size:12px;margin:-4px 0 10px}
svg{width:100%;display:block;overflow:visible}
.exon{cursor:pointer}.exon:hover rect{stroke:var(--ink);stroke-width:1.5}
.sel-exon rect{stroke:var(--hit)!important;stroke-width:2!important}
table{border-collapse:collapse;width:100%;font-size:12.5px;font-variant-numeric:tabular-nums}
th{text-align:left;background:var(--accent-soft);padding:7px 9px;position:sticky;top:0;font-weight:650;white-space:nowrap}
td{padding:6px 9px;border-top:1px solid var(--line);white-space:nowrap;font-family:ui-monospace,Menlo,monospace;font-size:12px}
tr:nth-child(even) td{background:var(--stripe)}
.tw{max-height:400px;overflow:auto;border:1px solid var(--line);border-radius:8px}
.sig{color:var(--hit);font-weight:700}.up{color:var(--ok)}
.pill{display:inline-block;padding:1px 7px;border-radius:99px;font-size:11px;background:var(--accent-soft);color:var(--accent);margin-left:6px}
.empty{color:var(--muted);font-style:italic;padding:8px 2px}
.legend{display:flex;gap:14px;flex-wrap:wrap;color:var(--muted);font-size:11.5px;margin-top:8px}
.sw{display:inline-block;width:11px;height:11px;border-radius:3px;vertical-align:-1px;margin-right:4px}
button.clr{background:none;border:1px solid var(--line);color:var(--muted);border-radius:7px;padding:3px 9px;font-size:11.5px;cursor:pointer}
button.clr:hover{border-color:var(--accent);color:var(--accent)}
</style></head><body><div class="app">
<aside>
 <div class="brand"><b>modulator</b><span>gene &amp; fragmentform browser</span></div>
 <div class="search"><input id="q" placeholder="Search gene or fragmentform id…" autocomplete="off"></div>
 <div id="list"></div>
</aside>
<main id="main"><div class="empty">Select a gene to begin.</div></main>
</div>
<script>
const DATA=__DATA__;
const $=s=>document.querySelector(s);
let cur=null, selExon=null;
const fmt=n=>n==null||isNaN(n)?"–":(+n).toLocaleString();
const pct=v=>v==null||isNaN(v)?"–":(100*v).toFixed(1)+"%";
const sci=p=>p==null||isNaN(p)?"–":(p<1e-4?p.toExponential(1):p.toFixed(4));

function renderList(f){
  const q=(f||"").trim().toLowerCase();
  const rows=DATA.index.filter(r=>!q||r.g.toLowerCase().includes(q)||
      (DATA.genes[r.g].forms||[]).some(x=>x.zt.toLowerCase().includes(q))).slice(0,400);
  $("#list").innerHTML=rows.map(r=>`<div class="gi${cur===r.g?' sel':''}" data-g="${r.g}">
     <b>${r.g}</b><small>${r.n} ff · ${fmt(r.r)}</small></div>`).join("")
     ||`<div class="empty" style="padding:14px">No match.</div>`;
  document.querySelectorAll(".gi").forEach(e=>e.onclick=()=>select(e.dataset.g));
}
function select(g){cur=g;selExon=null;renderList($("#q").value);draw();}

function draw(){
  const G=DATA.genes[cur]; if(!G){return;}
  const all=G.forms.flatMap(f=>f.exons), lo=Math.min(...all.map(e=>e[0])), hi=Math.max(...all.map(e=>e[1]));
  const W=1000,H=Math.max(60,G.forms.length*30+22),PAD=8,span=Math.max(hi-lo,1);
  const X=p=>PAD+(p-lo)/span*(W-2*PAD);
  const rows=G.forms.map((f,i)=>{
    const y=i*30+16, mid=y+8;
    const introns=`<line x1="${X(Math.min(...f.exons.map(e=>e[0])))}" x2="${X(Math.max(...f.exons.map(e=>e[1])))}"
       y1="${mid}" y2="${mid}" stroke="var(--line)" stroke-width="2"/>`;
    const ex=f.exons.map((e,j)=>`<g class="exon" data-s="${e[0]}" data-e="${e[1]}" data-zt="${f.zt}">
        <rect x="${X(e[0])}" y="${y}" width="${Math.max(X(e[1])-X(e[0]),2)}" height="16" rx="2.5"
          fill="${j%2?'var(--exon-alt)':'var(--exon)'}"><title>${f.zt} exon ${j+1}\n${fmt(e[0])}–${fmt(e[1])} (${fmt(e[1]-e[0])} nt)\nclick to filter sites</title></rect></g>`).join("");
    const lab=`<text x="0" y="${mid+4}" font-size="10.5" fill="var(--muted)" font-family="ui-monospace,monospace">${f.zt.split('.').slice(-2).join('.')}</text>`;
    return `<g>${introns}${ex}</g>`+`<g transform="translate(${W+8},0)">${lab}</g>`;
  }).join("");
  const arrow=G.strand==="+"?"5′ → 3′":"3′ ← 5′";
  const ffTable=G.forms.map(f=>`<tr><td>${f.zt}</td><td>${f.classification||"–"}</td><td>${fmt(f.reads)}</td>
      <td>${f.pas||"–"}</td><td>${f.tail?f.tail.toFixed(0)+" nt":"–"}</td><td>${f.exons.length}</td></tr>`).join("");
  $("#main").innerHTML=`
   <h1>${G.gene}</h1>
   <div class="sub">${G.chrom} · ${G.strand} strand (${arrow}) · ${G.forms.length} fragmentforms · ${fmt(G.reads)} reads</div>
   <div class="card"><h2>Fragmentform structures</h2>
     <div class="hint">Click an exon to filter the tables below to modification sites inside it.
       <span id="selinfo"></span></div>
     <div style="overflow-x:auto"><svg viewBox="0 0 ${W+150} ${H}" style="min-width:640px">${rows}</svg></div>
     <div class="legend"><span><i class="sw" style="background:var(--exon)"></i>exon</span>
       <span><i class="sw" style="background:var(--line)"></i>intron</span>
       <span>hover an exon for coordinates</span></div>
   </div>
   <div class="card"><h2>Fragmentforms</h2><div class="tw"><table>
     <thead><tr><th>fragmentform</th><th>class</th><th>reads</th><th>PAS</th><th>median tail</th><th>exons</th></tr></thead>
     <tbody>${ffTable}</tbody></table></div></div>
   <div id="tables"></div>`;
  document.querySelectorAll(".exon").forEach(el=>el.onclick=()=>{
    const s=+el.dataset.s,e=+el.dataset.e;
    selExon=(selExon&&selExon[0]===s&&selExon[1]===e)?null:[s,e];
    document.querySelectorAll(".exon").forEach(x=>x.classList.remove("sel-exon"));
    if(selExon) document.querySelectorAll(`.exon[data-s="${s}"][data-e="${e}"]`).forEach(x=>x.classList.add("sel-exon"));
    tables();
  });
  tables();
}

function inSel(p){return !selExon || (p>=selExon[0] && p<selExon[1]);}
function tables(){
  const G=DATA.genes[cur];
  $("#selinfo").innerHTML = selExon
    ? `<span class="pill">exon ${fmt(selExon[0])}–${fmt(selExon[1])}</span>
       <button class="clr" onclick="clearSel()">clear</button>`
    : `<span class="pill">all sites</span>`;
  const S=G.sites.filter(s=>inSel(s.pos)).sort((a,b)=>a.pos-b.pos);
  const D=G.diffs.filter(s=>inSel(s.pos)).sort((a,b)=>a.padj-b.padj);
  const C=G.cond.filter(s=>inSel(s.pos)).sort((a,b)=>a.padj-b.padj);
  const Hh=G.hier.filter(s=>inSel(s.pos)).sort((a,b)=>a.padj-b.padj);
  const tbl=(title,hint,head,body)=>`<div class="card"><h2>${title}</h2>${hint?`<div class="hint">${hint}</div>`:""}
     ${body?`<div class="tw"><table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table></div>`
            :`<div class="empty">No rows${selExon?" in the selected exon":""}.</div>`}</div>`;
  $("#tables").innerHTML =
   tbl(`Modification sites <span class="pill">${S.length}</span>`,
       "Pooled stoichiometry per site (per ZN transcript partition where available).",
       "<th>position</th><th>mod</th><th>ZN</th><th>coverage</th><th>modified</th>",
       S.slice(0,400).map(s=>`<tr><td>${fmt(s.pos)}</td><td>${s.mod}</td><td>${s.zn<0?"–":s.zn}</td>
         <td>${fmt(s.cov)}</td><td>${pct(s.frac)}</td></tr>`).join(""))
 + tbl(`Differential between transcripts <span class="pill">${D.length}</span>`,
       "Sites whose modified fraction differs between this gene's transcript partitions.",
       "<th>position</th><th>mod</th><th>effect</th><th>FDR</th>",
       D.slice(0,200).map(d=>`<tr><td>${fmt(d.pos)}</td><td>${d.mod}</td><td>${pct(d.effect)}</td>
         <td class="${d.padj<0.05?'sig':''}">${sci(d.padj)}</td></tr>`).join(""))
 + (C.length||G.cond.length? tbl(`Between conditions <span class="pill">${C.length}</span>`,
       "Replicate-aware differential modification between conditions.",
       "<th>position</th><th>mod</th><th>contrast</th><th>delta</th><th>FDR</th>",
       C.slice(0,200).map(c=>`<tr><td>${fmt(c.pos)}</td><td>${c.mod}</td><td>${c.contrast}</td>
         <td class="${c.delta>0?'up':''}">${(c.delta>0?"+":"")+(100*c.delta).toFixed(1)}%</td>
         <td class="${c.padj<0.05?'sig':''}">${sci(c.padj)}</td></tr>`).join("")):"")
 + (Hh.length||G.hier.length? tbl(`Truncation-aware fragmentform comparison <span class="pill">${Hh.length}</span>`,
       "Only reads that demonstrably span each pair's divergence point — <i>n informative</i> shows the power that survived.",
       "<th>position</th><th>A</th><th>B</th><th>delta</th><th>n inf.</th><th>div. from 3′</th><th>FDR</th>",
       Hh.slice(0,200).map(h=>`<tr><td>${fmt(h.pos)}</td><td>${h.a.split('.').slice(-1)}</td><td>${h.b.split('.').slice(-1)}</td>
         <td>${(h.delta>0?"+":"")+(100*h.delta).toFixed(1)}%</td><td>${fmt(h.ninf)}</td><td>${fmt(h.div3p)}</td>
         <td class="${h.padj<0.05?'sig':''}">${sci(h.padj)}</td></tr>`).join("")):"");
}
function clearSel(){selExon=null;document.querySelectorAll(".exon").forEach(x=>x.classList.remove("sel-exon"));tables();}
$("#q").addEventListener("input",e=>renderList(e.target.value));
renderList("");
if(DATA.index.length) select(DATA.index[0].g);
</script></body></html>"""


if __name__ == "__main__":
    main()
