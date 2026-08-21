#!/usr/bin/env python3
"""
Post-differential granular classification of between-transcript (ZN) m6A sites.

This runs AFTER ``test_stoichiometry_diffs.py``. For every significant ZN site
(by default: BH-FDR < 0.05 AND effect = max |Δ stoichiometry| >= 0.10, i.e. the
">10% absolute" rule) it assigns a single, MECE structural category that
explains *why* the isoforms differ in m6A at that base, anchored to the LONGEST
3'UTR isoform of the gene.

INPUTS
------
--diff-tsv : ``{prefix}__ZN_site_diff_results.tsv`` from test_stoichiometry_diffs
             (cols: gene_name, mod_code, chrom, start0, end0, strand,
             n_tx_tested, effect_max_abs_frac_diff, p_adj_bh, per_transcript_json
             where each entry is {"ZN":int,"Ncov":int,"Nmod":int,"frac":float}).
--gtf      : the assembled, read-backed GTF ``{prefix}.gtf``. Exon features carry
             ref_gene_name, zn_index, tes -> isoform exon models keyed by
             (gene, zn_index). ``zn_index`` is what the ZN partition tag (hence
             ``ZN_transcript_index`` in the diff table) is built from.

EXON SOURCE: exon models come from the GTF ``exon`` features (NOT the
``intron_chain`` attribute, which lists INTRONS). Architecture
(IPA / TANDEM_APA / FULL_LENGTH / DISTAL_EXT / REFERENCE) is recomputed from the
real terminal exon, exactly as in the offline analysis (step7).

POSITION STATUS of a genomic site within one isoform:
  exonic_internal  - in an exon that is NOT the 3'-terminal exon
  exonic_terminal  - in the 3'-terminal exon (incl. 3'UTR)
  intronic         - within the transcript span but spliced out (a gap)
  absent           - outside [first_exon_start, last_exon_end]

CATEGORIES (keyed on the HIGH-m6A isoform's architecture, the site's status in
the ANCHOR (longest) isoform, the hi/lo APA direction, and distance to the
nearest spliced junction, EJC_NT=150):

 -- A is PRIVATE to the high isoform (not in the long mature transcript) --
  IPA_UNIQUE              hi = IPA isoform; site exonic_TERMINAL in hi but
                          intronic/absent in the anchor -> the A only exists in
                          the mature IPA transcript (intron-derived terminal
                          exon). IPA-private; cleavage-dependent.
  SPLICED_EXON_UNIQUE     A in an internal/cassette exon present in hi but spliced
                          out (intronic) / absent in the anchor or low comparator.
  LAST_EXON_DISTAL_ONLY   A is in the anchor but the low (proximal) isoform's TES
                          is upstream of it -> distal-3'UTR-private A.

 -- A is SHARED (exonic in hi, lo and anchor); terminalization / EJC --
  IPA_SHARED_EJC          hi = IPA; A exonic_TERMINAL in hi but exonic_INTERNAL in
                          the long anchor -> same A in both mature RNAs, gains m6A
                          in IPA because the downstream EJC is gone (cleavage-dep).
  SPLICING_EJC            shared A, non-IPA: terminalized, or a junction within
                          EJC_NT in the low/anchor is removed in hi -> EJC relief.

 -- A SHARED in the SAME terminal exon; tandem 3'UTR APA (no cleavage diff) --
  LAST_EXON_PROXIMAL_APA_FAVORED  same last exon (same acceptor), different TES,
                          PROXIMAL (shorter 3'UTR) isoform carries more m6A.
  LAST_EXON_DISTAL_APA_FAVORED    same geometry; DISTAL (longer 3'UTR) favored.

 -- A SHARED via an ALTERNATIVE / SEPARATE terminal exon --
  ALTERNATIVE_LAST_EXON   hi & lo terminal at the site but their terminal exons
                          start at DIFFERENT (but nearby/overlapping) acceptors.
  INTERGENIC_TERMINAL_EXON  (NEW) the site is exonic_TERMINAL in a non-IPA high
                          isoform whose terminal exon is genomically DISJOINT from
                          the comparator's terminal exon and separated by a large
                          gap (>= --intergenic-gap, default 1000 bp). The contrast
                          spans spatially-separated terminal exons (read-through /
                          downstream-independent / intergenic-scale alternative
                          last exons) rather than a local 3'UTR-length (tandem APA)
                          or a near alternative-last-exon. Absorbs a large share of
                          what otherwise lands in SPLICED_EXON_UNIQUE /
                          LAST_EXON_DISTAL_ONLY / ALTERNATIVE_LAST_EXON / residual.

 -- A SHARED with NO local 3' structural difference --
  SHARED_TERMINAL_EXON    hi & lo share the SAME terminal exon AND the SAME 3'
                          cleavage (TES); m6A tracks isoform identity, not APA/EJC.
  SHARED_INTERNAL_EXON    site in a constitutive INTERNAL exon, no junction
                          asymmetry -> not attributable to 3' architecture.

 -- leftover / artifact --
  UNEXPLAINED_SHARED      rare residual (terminal in hi, internal in lo, no nearby
                          differential junction).
  HI_INTRONIC_ARTIFACT    the high-m6A isoform does NOT structurally contain the A
                          (intronic/absent) -> the "high" stoich is intron-read
                          noise (NOT a real isoform-specific site).

OUTPUT: ``{prefix}__ZN_site_classified.tsv`` = the diff rows that pass the
significance + effect gate, augmented with hi/lo isoform identity & architecture,
the site's status in hi/lo/anchor, junction distances, and the category.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import sys
from collections import Counter, defaultdict

csv.field_size_limit(10**9)

# attribute regexes for the GTF (same convention as the assembler output)
_AT = {k: re.compile(k + r' "([^"]*)"') for k in
       ('ref_gene_name', 'gene_id', 'zn_index', 'transcript_index', 'tes',
        'read_support')}


def _attr(s, key):
    m = _AT[key].search(s)
    return m.group(1) if m else None


# ----------------------------- isoform models ----------------------------- #

def classify_arch(tes, exons, strand, ref_tes, ref_exons, tes_tol, inside_tol):
    """architecture of one isoform vs the gene reference (longest), from REAL exons."""
    if not ref_exons:
        return 'AMBIGUOUS'
    if strand == '+':
        upstream = ref_tes - tes
        ref_term_start = ref_exons[-1][0]
        inside_term = (ref_term_start - inside_tol) <= tes <= (ref_tes + inside_tol)
        gmin, gmax = ref_exons[0][0], max(ref_exons[-1][1], ref_tes)
    else:
        upstream = tes - ref_tes
        ref_term_end = ref_exons[0][1]
        inside_term = (ref_tes - inside_tol) <= tes <= (ref_term_end + inside_tol)
        gmin, gmax = min(ref_exons[0][0], ref_tes), ref_exons[-1][1]
    if abs(upstream) <= tes_tol:
        return 'FULL_LENGTH'
    if upstream < 0:
        return 'DISTAL_EXT'
    if inside_term:
        return 'TANDEM_APA'              # proximal TES inside reference terminal exon
    if gmin - inside_tol <= tes <= gmax + inside_tol:
        return 'IPA'                     # proximal TES upstream of terminal exon
    return 'AMBIGUOUS'


def load_isoforms(gtf_path, tes_tol, inside_tol):
    """iso[(gene, zn)] = dict(strand, chrom, tes, exons=[(s,e)...], rs, arch).

    Keyed by (ref_gene_name, zn_index): zn_index is what the ZN partition tag --
    hence the "ZN" key in the diff table per_transcript_json -- is built from.

    Version tolerance: older modulator GTFs (pre-metagene partitioning) carry no
    ``zn_index`` attribute; there the diff table's "ZN" == the per-gene
    ``transcript_index``, so we fall back to ``transcript_index`` when ``zn_index``
    is absent. (A single GTF is internally consistent: every line either has
    ``zn_index`` or none do, so the per-line fallback never mixes the two.)
    """
    exons = defaultdict(list)
    meta = {}
    with open(gtf_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            typ = f[2]; attr = f[8]
            g = _attr(attr, 'ref_gene_name') or _attr(attr, 'gene_id')
            zn = _attr(attr, 'zn_index') or _attr(attr, 'transcript_index')
            if not g or not zn:
                continue
            key = (g, zn)
            if typ == 'exon':
                exons[key].append((int(f[3]), int(f[4])))
            elif typ == 'transcript':
                tes = _attr(attr, 'tes'); rs = _attr(attr, 'read_support')
                meta[key] = dict(strand=f[6], chrom=f[0],
                                 tes=int(tes) if tes and tes.lstrip('-').isdigit() else None,
                                 rs=int(rs) if rs and rs.isdigit() else 0)
    iso = {}
    for key, exl in exons.items():
        exl.sort()
        m = meta.get(key, {})
        strand = m.get('strand', '+')
        tes = m.get('tes')
        if tes is None:
            tes = exl[-1][1] if strand == '+' else exl[0][0]
        iso[key] = dict(strand=strand, chrom=m.get('chrom'), tes=tes, exons=exl,
                        rs=m.get('rs', 0), arch=None)
    # per-gene reference (longest 3'UTR) + architecture
    genes = defaultdict(list)
    for (g, zn) in iso:
        genes[g].append(zn)
    for g, zns in genes.items():
        strand = iso[(g, zns[0])]['strand']
        if strand == '+':
            ref = max(zns, key=lambda z: (iso[(g, z)]['tes'], iso[(g, z)]['rs']))
        else:
            ref = min(zns, key=lambda z: (iso[(g, z)]['tes'], -iso[(g, z)]['rs']))
        rd = iso[(g, ref)]
        for zn in zns:
            d = iso[(g, zn)]
            d['arch'] = 'REFERENCE' if zn == ref else classify_arch(
                d['tes'], d['exons'], strand, rd['tes'], rd['exons'], tes_tol, inside_tol)
    return iso, genes


def status_in(exons, strand, pos):
    lo = exons[0][0]; hi = exons[-1][1]
    if pos < lo or pos > hi:
        return 'absent'
    term_idx = len(exons) - 1 if strand == '+' else 0
    for i, (s, e) in enumerate(exons):
        if s <= pos <= e:
            return 'exonic_terminal' if i == term_idx else 'exonic_internal'
    return 'intronic'


def anchor_of(gene, zns, iso):
    """longest-3'UTR isoform = most distal TES (tie-break read_support)."""
    best = None
    for zn in zns:
        d = iso[(gene, zn)]
        if d['tes'] is None:
            continue
        key = (d['tes'] if d['strand'] == '+' else -d['tes'], d['rs'])
        if best is None or key > best[0]:
            best = (key, zn)
    return best[1] if best else None


def terminal_exon(exons, strand):
    return exons[-1] if strand == '+' else exons[0]


def internal_junctions(exons):
    js = []
    for i in range(len(exons) - 1):
        js.append(exons[i][1]); js.append(exons[i + 1][0])
    return js


def dist_to_junction(exons, pos):
    return min((abs(pos - j) for j in internal_junctions(exons)), default=10**9)


def _intron_inside(exons, lo, hi):
    """True if ``exons`` has an intron (the gap between two consecutive exons) lying ENTIRELY within
    the open interval (lo, hi). Used to tell intron RETENTION (the other form splices out a whole
    intron that this form's single exon spans) from a mere alternative donor/acceptor (a single
    shifted boundary that still rejoins downstream)."""
    for i in range(len(exons) - 1):
        istart = exons[i][1]      # donor: end of exon i
        iend = exons[i + 1][0]    # acceptor: start of exon i+1
        if lo < istart and iend < hi and istart < iend:   # require a real (positive-length) intron
            return True
    return False


def is_proximal(tes_a, tes_b, strand):
    return (tes_a < tes_b) if strand == '+' else (tes_a > tes_b)


# --- Orthogonal stoichiometry axes layered on top of the structural `category` ---
# The primary `category` is the structural MECE label. These add, in separate columns, the
# stoichiometry RELATIONSHIP the structural label does not encode: which fragmentform is more
# modified (direction), how large the gap is (tier), and whether the favored form is itself
# high- or low-stoichiometry (level). This answers "split by high/low stoichiometry" without
# multiplying the primary label.
STRUCTURAL_OF = {
    'IPA_UNIQUE': 'INTRONIC_POLYADENYLATION', 'IPA_SHARED_EJC': 'INTRONIC_POLYADENYLATION',
    'SPLICING_EJC': 'EJC_SPLICING', 'SPLICED_EXON_UNIQUE': 'CASSETTE_EXON',
    'LAST_EXON_PROXIMAL_APA_FAVORED': 'TANDEM_APA', 'LAST_EXON_DISTAL_APA_FAVORED': 'TANDEM_APA',
    'LAST_EXON_DISTAL_ONLY': 'TANDEM_APA', 'ALTERNATIVE_LAST_EXON': 'ALTERNATIVE_LAST_EXON',
    'INTERGENIC_TERMINAL_EXON': 'INTERGENIC_TERMINAL_EXON',
    'SHARED_TERMINAL_EXON': 'SHARED_TERMINAL_EXON', 'SHARED_INTERNAL_EXON': 'SHARED_INTERNAL_EXON',
    'UNEXPLAINED_SHARED': 'UNEXPLAINED', 'HI_INTRONIC_ARTIFACT': 'ARTIFACT', 'UNCLASSIFIED': 'UNCLASSIFIED',
}


def stoich_tier(delta):
    d = abs(delta)
    return ('T1_MARGINAL' if d < 0.25 else 'T2_MODERATE' if d < 0.50
            else 'T3_STRONG' if d < 0.75 else 'T4_NEAR_BINARY')


def hi_stoich_level(f):
    """Absolute modification level of the favored (higher) fragmentform -- the 'is the alt last
    exon itself high- or low-stoichiometry' axis."""
    return 'HI_HYPER' if f >= 0.66 else 'HI_INTERMED' if f >= 0.33 else 'HI_HYPO'


def stoich_direction(gene, hiZN, loZN, iso, cat, tes_tol):
    """(direction, context-aware alias). Direction is the universal 3'UTR-length polarity of the
    MORE-modified fragmentform; the alias renames it per structural mechanism."""
    ihi = iso.get((gene, hiZN)); ilo = iso.get((gene, loZN))
    if ihi is None or ilo is None:
        return '', ''
    ht, lt = ihi['tes'], ilo['tes']
    if ht is None or lt is None or abs(ht - lt) <= tes_tol:
        base = 'CO_TERMINAL'
    elif is_proximal(ht, lt, ihi['strand']):
        base = 'PROXIMAL_HIGHER'          # shorter-3'UTR fragmentform carries more modification
    else:
        base = 'DISTAL_HIGHER'            # longer-3'UTR fragmentform carries more modification
    if cat in ('IPA_UNIQUE', 'IPA_SHARED_EJC'):
        ctx = 'IPA_FORM_HIGHER' if ihi['arch'] == 'IPA' else 'FULLLENGTH_HIGHER'
    elif cat == 'SPLICED_EXON_UNIQUE':
        ctx = 'INCLUDED_ISOFORM_HIGHER'   # hi structurally must contain the base
    elif cat == 'SPLICING_EJC':
        ctx = 'EJC_REMOVED_HIGHER'
    else:
        ctx = base
    return base, ctx


def make_class_key(structural_category, stoich_direction):
    """The single primary classification key = mechanism + which fragmentform is more modified,
    e.g. TANDEM_APA__PROXIMAL_HIGHER. This REPLACES the old fused 14-label `category`. Rows with no
    direction (UNCLASSIFIED, or missing isoform models) key on structural_category alone."""
    return f"{structural_category}__{stoich_direction}" if stoich_direction else structural_category


def same_last_exon_start(th, tl, strand, inside_tol):
    a = th[0] if strand == '+' else th[1]
    b = tl[0] if strand == '+' else tl[1]
    return abs(a - b) <= inside_tol


def exon_gap(a, b):
    """genomic gap between two (start,end) exons; 0 if they overlap."""
    if a[1] < b[0]:
        return b[0] - a[1]
    if b[1] < a[0]:
        return a[0] - b[1]
    return 0


# =============================================================================
# Taxonomy tree (4 top-level buckets):
#   PRIVATE        - the modified base exists (exonically) ONLY in the higher fragmentform; the
#                    difference is trivial (the base is physically absent from the lower form).
#   SHARED_LOCAL   - base exonic in BOTH; the distinguishing structural change is ON / adjacent to
#                    the base's own exon.
#   SHARED_DISTAL  - base exonic in both, in an identical local context; the forms differ only
#                    ELSEWHERE, so the modification tracks isoform identity.
#   UNEXPLAINABLE  - cannot be attributed structurally (5' blind spot, intron-read/soft-clip
#                    artifact, or no usable isoform model).
# class_key = BUCKET__EVENT__DIRECTION
# =============================================================================
TAXONOMY = {
    "PRIVATE":       ["SKIPPED_EXON", "INTRONIC_POLYA", "ALT_LAST_EXON"],
    "SHARED_LOCAL":  ["ALT_DONOR", "ALT_ACCEPTOR", "ALT_POLYA_SITE", "RETAINED_INTRON",
                      "IPA_EXTENSION", "NEAR_ALT_JUNCTION"],
    "SHARED_DISTAL": ["DISTAL_APA", "DISTAL_SPLICING"],
    "UNEXPLAINABLE": ["FIVE_PRIME_UNCERTAIN", "INTRON_READ_ARTIFACT", "NO_MODEL", "UNRESOLVED"],
}
BUCKET_ORDER = ["PRIVATE", "SHARED_LOCAL", "SHARED_DISTAL", "UNEXPLAINABLE"]


def _exon_containing(exons, pos):
    for (s, e) in exons:
        if s <= pos <= e:
            return (s, e)
    return None


def _exon_index(exons, pos):
    for i, (s, e) in enumerate(exons):
        if s <= pos <= e:
            return i
    return None


def _is_last_exon(exons, strand, idx):
    return idx == (len(exons) - 1 if strand == '+' else 0)


def _is_first_exon(exons, strand, idx):
    return idx == (0 if strand == '+' else len(exons) - 1)


def classify_tree(gene, pos, hiZN, loZN, iso, *, tes_tol, ejc_nt):
    """Assign the site to (bucket, event, direction). Uses only transcript-MODEL coordinates -- a
    PRIVATE call means the base is intronic/absent in the low form's assembled model (no spanning-
    read requirement), guarded only against the 5' blind spot where the model is unreliable."""
    ihi = iso[(gene, hiZN)]; ilo = iso[(gene, loZN)]
    strand = ihi['strand']
    sh = status_in(ihi['exons'], strand, pos)
    sl = status_in(ilo['exons'], strand, pos)
    info = dict(status_hi=sh, status_lo=sl,
                jd_hi=dist_to_junction(ihi['exons'], pos),
                jd_lo=dist_to_junction(ilo['exons'], pos), delta_nt='')

    # (0) the high (more-modified) form must structurally contain the base; else it is an intron /
    # soft-clip read artifact, not a real exonic modification.
    if sh not in ('exonic_terminal', 'exonic_internal'):
        return 'UNEXPLAINABLE', 'INTRON_READ_ARTIFACT', '', info

    # ---- LEVEL 1: PRIVATE (base intronic/absent in the low form, by model coordinates) ----
    if sl in ('intronic', 'absent'):
        # 5' guard: an "absent" call 5' of the low model's own 5' end is truncation, not real
        # absence -- we can't trust the model there. Everything else (intronic, or absent past the
        # reliable 3' end) is a genuine private call.
        lo5 = ilo['exons'][0][0] if strand == '+' else ilo['exons'][-1][1]
        truncated_5p = (pos < lo5) if strand == '+' else (pos > lo5)
        if sl == 'absent' and truncated_5p:
            return 'UNEXPLAINABLE', 'FIVE_PRIME_UNCERTAIN', '', info
        if ihi['arch'] == 'IPA' and sh == 'exonic_terminal':
            return 'PRIVATE', 'INTRONIC_POLYA', 'IPA_TRANSCRIPT_HIGHER', info
        if sh == 'exonic_terminal':
            th, tl = terminal_exon(ihi['exons'], strand), terminal_exon(ilo['exons'], strand)
            info['delta_nt'] = exon_gap(th, tl)
            d = 'LONGER_EXON_HIGHER' if (th[1] - th[0]) >= (tl[1] - tl[0]) else 'SHORTER_EXON_HIGHER'
            return 'PRIVATE', 'ALT_LAST_EXON', d, info
        ehi = _exon_containing(ihi['exons'], pos)
        if ehi:
            info['delta_nt'] = ehi[1] - ehi[0]
        return 'PRIVATE', 'SKIPPED_EXON', 'WITH_EXON_HIGHER', info

    # ---- base is SHARED (exonic in both). LEVEL 2: LOCAL vs DISTAL ----
    ehi = _exon_containing(ihi['exons'], pos); elo = _exon_containing(ilo['exons'], pos)
    ihx = _exon_index(ihi['exons'], pos); ilx = _exon_index(ilo['exons'], pos)
    hi_last = _is_last_exon(ihi['exons'], strand, ihx); lo_last = _is_last_exon(ilo['exons'], strand, ilx)
    hi_first = _is_first_exon(ihi['exons'], strand, ihx); lo_first = _is_first_exon(ilo['exons'], strand, ilx)

    # 3'UTR-length polarity of the higher form (used by the shared events)
    ht, lt = ihi['tes'], ilo['tes']
    if ht is None or lt is None or abs(ht - lt) <= tes_tol:
        polar = 'CO_TERMINAL_HIGHER'
    elif is_proximal(ht, lt, strand):
        polar = 'PROXIMAL_HIGHER'
    else:
        polar = 'DISTAL_HIGHER'

    # (a) SAME last exon (shared acceptor), different poly(A) site -> ALT_POLYA_SITE (tandem APA).
    # A DIFFERENT acceptor means a mutually-exclusive alternative last exon, which must fall through
    # to (b) so the acceptor/splicing shift is named rather than mislabeled as tandem APA.
    if (hi_last and lo_last and ht is not None and lt is not None and abs(ht - lt) > tes_tol
            and ehi and elo and same_last_exon_start(ehi, elo, strand, tes_tol)):
        info['delta_nt'] = abs(ht - lt)
        return 'SHARED_LOCAL', 'ALT_POLYA_SITE', polar, info

    # (b) the base's OWN exon has a shifted boundary between the two forms. This stays LOCAL even when
    # the base's exon is terminal in one form -- the change is at the base. We separate three
    # mechanisms rather than lumping every longer exon into "alt donor":
    #   IPA_EXTENSION   the longer form reads into the intron and polyadenylates there (its extended
    #                   donor boundary == its own, more-proximal poly(A) site);
    #   RETAINED_INTRON the shorter form splices out a whole intron that the longer form's single exon
    #                   spans (both intron ends fall inside the longer exon);
    #   ALT_DONOR / ALT_ACCEPTOR  a single shifted splice boundary that still rejoins downstream.
    if ehi and elo:
        if strand == '+':
            acc_hi, acc_lo, don_hi, don_lo = ehi[0], elo[0], ehi[1], elo[1]
        else:
            acc_hi, acc_lo, don_hi, don_lo = ehi[1], elo[1], ehi[0], elo[0]
        len_hi, len_lo = ehi[1] - ehi[0], elo[1] - elo[0]
        hi_longer = len_hi >= len_lo
        long_exon = ehi if hi_longer else elo
        short_exons = ilo['exons'] if hi_longer else ihi['exons']
        long_don = don_hi if hi_longer else don_lo
        long_tes = ht if hi_longer else lt
        short_tes = lt if hi_longer else ht
        long_is_last = hi_last if hi_longer else lo_last
        exon_dir = 'LONGER_EXON_HIGHER' if hi_longer else 'SHORTER_EXON_HIGHER'
        acc_diff = (acc_hi != acc_lo) and not (hi_first or lo_first)   # 5' end is the blind spot
        don_diff = (don_hi != don_lo)
        if acc_diff or don_diff:
            # IPA_EXTENSION: the longer form's base exon IS its terminal exon and ends at its own
            # poly(A) site (reads into the intron and terminates there), NOT co-terminal with the other
            # form. Requiring the base exon to be terminal prevents an internal exon whose splice donor
            # merely happens to sit within tes_tol of a short downstream terminal exon from being
            # mislabeled IPA (that is really an ALT_DONOR).
            not_coterm = long_tes is not None and (short_tes is None or abs(long_tes - short_tes) > tes_tol)
            if (don_diff and long_is_last and long_tes is not None
                    and abs(long_don - long_tes) <= tes_tol and not_coterm):
                info['delta_nt'] = abs(don_hi - don_lo)
                d = 'EXTENDS_TO_PA_HIGHER' if hi_longer else 'EXTENDS_TO_PA_LOWER'
                return 'SHARED_LOCAL', 'IPA_EXTENSION', d, info
            # RETAINED_INTRON: the shorter form has a whole intron inside the longer form's base exon.
            if _intron_inside(short_exons, long_exon[0], long_exon[1]):
                info['delta_nt'] = abs(len_hi - len_lo)
                d = 'INTRON_RETAINED_HIGHER' if hi_longer else 'INTRON_RETAINED_LOWER'
                return 'SHARED_LOCAL', 'RETAINED_INTRON', d, info
            # otherwise a single shifted splice boundary that rejoins downstream.
            if acc_diff:
                info['delta_nt'] = abs(acc_hi - acc_lo)
                return 'SHARED_LOCAL', 'ALT_ACCEPTOR', exon_dir, info
            if don_diff:
                info['delta_nt'] = abs(don_hi - don_lo)
                return 'SHARED_LOCAL', 'ALT_DONOR', exon_dir, info

    # (c) base near a junction present in one form but not the other (EJC window) -> NEAR_ALT_JUNCTION
    if (info['jd_hi'] <= ejc_nt) != (info['jd_lo'] <= ejc_nt):
        d = 'JUNCTION_REMOVED_HIGHER' if info['jd_hi'] > ejc_nt else 'JUNCTION_PRESENT_HIGHER'
        return 'SHARED_LOCAL', 'NEAR_ALT_JUNCTION', d, info

    # (d) the base's exon is identical in both -> the structural difference is DISTAL. Sub-classify
    # by WHERE it is: a different 3' end (distal APA) or a splicing difference elsewhere (same 3' end).
    if ht is not None and lt is not None and abs(ht - lt) > tes_tol:
        return 'SHARED_DISTAL', 'DISTAL_APA', polar, info
    return 'SHARED_DISTAL', 'DISTAL_SPLICING', 'CO_TERMINAL_HIGHER', info


CATEGORY_ORDER = [
    'IPA_UNIQUE', 'SPLICED_EXON_UNIQUE', 'LAST_EXON_DISTAL_ONLY',
    'IPA_SHARED_EJC', 'SPLICING_EJC',
    'LAST_EXON_PROXIMAL_APA_FAVORED', 'LAST_EXON_DISTAL_APA_FAVORED',
    'ALTERNATIVE_LAST_EXON', 'INTERGENIC_TERMINAL_EXON',
    'SHARED_TERMINAL_EXON', 'SHARED_INTERNAL_EXON',
    'UNEXPLAINED_SHARED', 'HI_INTRONIC_ARTIFACT', 'UNCLASSIFIED',
]


def scan_private_sites(zn_long_path, iso, genes, *, min_frac, min_cov):
    """Coverage-INDEPENDENT PRIVATE detection.

    The differential test only reaches a site when >=2 fragmentforms are covered there, so a base
    that is genuinely ABSENT from a fragmentform (the cleanest private events) is never tested. This
    scan bypasses that: for every modified site carried by >=1 fragmentform (frac>=min_frac & pooled
    cov>=min_cov, exonic there), it checks the base's status by transcript-MODEL coordinates in every
    OTHER fragmentform of the gene. If the base is intronic/absent (and not in the 5' blind spot) in
    >=1 other form, the site is PRIVATE -- reported without requiring the absent form to have any
    coverage. Returns (rows, out_cols)."""
    out_cols = ['gene_name', 'mod_code', 'chrom', 'start0', 'end0', 'strand',
                'bucket', 'event', 'direction', 'structural_delta_nt',
                'carry_ZN', 'carry_arch', 'carry_frac', 'carry_cov',
                'n_forms_present', 'n_forms_absent', 'absent_in_ZN']
    agg = defaultdict(lambda: [0, 0])   # (gene,zn,chrom,start0,strand,mod) -> [Nmod, Ncov]
    with open(zn_long_path) as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        for r in rd:
            try:
                nm = int(float(r['Nmod'])); nc = int(float(r['Nvalid_cov']))
            except (KeyError, ValueError, TypeError):
                continue
            k = (r.get('gene_name', ''), str(r['ZN_transcript_index']), r['chrom'],
                 int(r['start0']), r['strand'], r['mod_code'])
            a = agg[k]; a[0] += nm; a[1] += nc
    sites = defaultdict(dict)
    for (g, zn, chrom, s0, strand, mod), (nm, nc) in agg.items():
        if nc > 0:
            sites[(g, chrom, s0, strand, mod)][zn] = (nm / nc, nc)

    rows = []
    for (g, chrom, s0, strand, mod), zdict in sites.items():
        pos = s0 + 1  # 1-based, to match the GTF exon coords used by status_in
        gene_zns = [str(z) for z in genes.get(g, []) if (g, str(z)) in iso]
        if len(gene_zns) < 2:
            continue
        present, absent = [], []
        for z in gene_zns:
            exons = iso[(g, z)]['exons']
            st = status_in(exons, strand, pos)
            if st in ('exonic_terminal', 'exonic_internal'):
                present.append(z)
            elif st in ('intronic', 'absent'):
                blind = pos <= exons[0][1] if strand == '+' else pos >= exons[-1][0]
                if not blind:
                    absent.append(z)
        if not absent:
            continue
        carriers = [(z, zdict[z][0], zdict[z][1]) for z in present
                    if z in zdict and zdict[z][0] >= min_frac and zdict[z][1] >= min_cov]
        if not carriers:
            continue
        carriers.sort(key=lambda t: t[1], reverse=True)   # highest modified fraction first
        cz, cfrac, ccov = carriers[0]
        ihi = iso[(g, cz)]; sh = status_in(ihi['exons'], strand, pos)
        if ihi['arch'] == 'IPA' and sh == 'exonic_terminal':
            event, direction, delta = 'INTRONIC_POLYA', 'IPA_TRANSCRIPT_HIGHER', ''
        elif sh == 'exonic_terminal':
            event, direction = 'ALT_LAST_EXON', ''
            th = terminal_exon(ihi['exons'], strand); tl = terminal_exon(iso[(g, absent[0])]['exons'], strand)
            delta = exon_gap(th, tl)
            direction = 'LONGER_EXON_HIGHER' if (th[1] - th[0]) >= (tl[1] - tl[0]) else 'SHORTER_EXON_HIGHER'
        else:
            e = _exon_containing(ihi['exons'], pos)
            event, direction, delta = 'SKIPPED_EXON', 'WITH_EXON_HIGHER', (e[1] - e[0]) if e else ''
        rows.append([g, mod, chrom, s0, s0 + 1, strand, 'PRIVATE', event, direction, delta,
                     cz, ihi['arch'], f"{cfrac:.4f}", ccov, len(present), len(absent), ','.join(absent)])
    rows.sort(key=lambda x: (x[7], x[0], x[3]))
    return rows, out_cols


def parse_args():
    ap = argparse.ArgumentParser(description="Granular structural classification of ZN diff sites.")
    ap.add_argument("--diff-tsv", required=True, help="{prefix}__ZN_site_diff_results.tsv")
    ap.add_argument("--gtf", required=True, help="assembled read-backed GTF {prefix}.gtf")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--min-effect", type=float, default=0.10,
                    help="min effect = max |Δ stoichiometry| (the >10%% rule). Default 0.10")
    ap.add_argument("--fdr", type=float, default=0.05, help="max BH-FDR (p_adj_bh). Default 0.05")
    ap.add_argument("--mod-filter", nargs="*", default=[],
                    help="mod_code(s) to classify. Default: empty = ALL modifications "
                         "(the diff table already carries every mod_code emitted upstream). "
                         "Pass e.g. --mod-filter a to restrict to m6A only.")
    ap.add_argument("--min-cov", type=int, default=0,
                    help="extra per-isoform Ncov floor on the JSON entries (the diff "
                         "table is already coverage-filtered by test_diffs --min-cov). Default 0")
    # --- coverage-independent PRIVATE-site scan (does not need the differential test; reuses
    # --zn-long defined below for the FILTERED_sites_long table) ---
    ap.add_argument("--private-out-tsv", default="",
                    help="output TSV for the coverage-independent PRIVATE-site scan: every modified "
                         "site is checked structurally against ALL fragmentforms of its gene, so "
                         "private events are found even when the form lacking the base has no coverage "
                         "there (the differential test cannot see those). Needs --zn-long.")
    ap.add_argument("--private-min-frac", type=float, default=0.10,
                    help="min modified fraction for a fragmentform to 'carry' the site in the PRIVATE scan.")
    ap.add_argument("--private-min-cov", type=int, default=20,
                    help="min pooled coverage in the carrying fragmentform for the PRIVATE scan.")
    ap.add_argument("--tes-tol", type=int, default=25,
                    help="TES match tolerance (bp); matches assembler.tes_match_tol. Default 25 "
                         "(was 200, which lumped sub-200bp tandem-APA into SHARED_TERMINAL_EXON).")
    ap.add_argument("--inside-tol", type=int, default=50)
    ap.add_argument("--ejc-nt", type=int, default=150)
    ap.add_argument("--intergenic-gap", type=int, default=1000,
                    help="min genomic gap (bp) between disjoint terminal exons to call "
                         "INTERGENIC_TERMINAL_EXON. Default 1000")
    ap.add_argument("--verbose", action="store_true")
    # --- optional per-category figure rendering (the "same kind of graphs") ---
    ap.add_argument("--zn-long", default="",
                    help="ZN filtered long table ({prefix}__ZN.filtered.long.tsv) with "
                         "per-sample Nvalid_cov/Nmod. When given together with --figs-dir, "
                         "renders the SAME 2-panel per-site stoichiometry figure as "
                         "test_stoichiometry_diffs for the top sites of EACH category.")
    ap.add_argument("--figs-dir", default="",
                    help="output directory for per-category figures (one subdir per "
                         "category). Requires --zn-long. Disabled when empty.")
    ap.add_argument("--figs-per-category", type=int, default=10,
                    help="max sites (ranked by effect) to plot per category. Default 10. "
                         "Applies to BOTH the stoichiometry (--figs-dir) and the "
                         "architecture-map (--arch-figs-dir) figures.")
    ap.add_argument("--arch-figs-dir", default="",
                    help="output directory for per-category ISOFORM ARCHITECTURE-MAP "
                         "figures (the locus track plot: every isoform of the gene drawn "
                         "as exon [blue=internal, orange=terminal/3'UTR] / intron models, "
                         "the site marked on each isoform [green o=exonic in mature RNA, "
                         "red x=intronic/absent], hi/lo/anchor tagged, per-isoform "
                         "stoichiometry labelled). Built from the GTF isoform models + the "
                         "diff table's per_transcript_json -- no --zn-long needed. One "
                         "subdir per category. Disabled when empty.")
    return ap.parse_args()


def render_category_figures(fig_records, zn_long_path, figs_dir, per_category,
                            verbose=False):
    """Render the SAME 2-panel per-site stoichiometry figure used by
    ``test_stoichiometry_diffs.make_plot`` for the top ``per_category`` sites
    (ranked by effect = max |Δ stoichiometry|) of EACH category.

    Figures land under ``{figs_dir}/{CATEGORY}/rankNN__gene__mod__locus.png`` so
    the HTML report can embed them grouped by category. Returns
    ``{category: n_figs_written}``.
    """
    if not fig_records or per_category <= 0:
        return {}
    by_cat = defaultdict(list)
    for rec in fig_records:
        by_cat[rec['class_key']].append(rec)

    # Headless backend MUST be set before the sibling imports pyplot.
    try:
        import matplotlib
        matplotlib.use("Agg")
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from test_stoichiometry_diffs import make_plot, site_key_tuple
        import pandas as pd
    except Exception as e:  # pragma: no cover - environment dependent
        print(f"[classify] WARN cannot set up figure rendering ({e}); skipping figures",
              file=sys.stderr)
        return {}

    df = pd.read_csv(zn_long_path, sep="\t", low_memory=False)
    need = {"gene_name", "mod_code", "chrom", "start0", "end0", "strand",
            "ZN_transcript_index", "sample", "Nvalid_cov", "Nmod"}
    missing = need - set(df.columns)
    if missing:
        print(f"[classify] WARN zn-long missing columns {sorted(missing)}; skipping figures",
              file=sys.stderr)
        return {}
    for c in ["start0", "end0", "Nvalid_cov", "Nmod", "ZN_transcript_index"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["start0", "end0", "Nvalid_cov", "Nmod",
                           "ZN_transcript_index"]).copy()
    df["start0"] = df["start0"].astype(int)
    df["end0"] = df["end0"].astype(int)
    df["ZN_transcript_index"] = df["ZN_transcript_index"].astype(int)
    df["site_key"] = df.apply(site_key_tuple, axis=1)
    df_idx = df.set_index("site_key")

    made = {}
    for cat, recs in by_cat.items():
        recs = sorted(recs, key=lambda d: d['effect'], reverse=True)[:per_category]
        cat_dir = os.path.join(figs_dir, cat)
        n = 0
        for rank, rec in enumerate(recs, 1):
            per_tx = rec['per_tx']
            if not per_tx or len(per_tx) < 2:
                continue
            key = (str(rec['gene']), str(rec['mod']), str(rec['chrom']),
                   int(rec['start0']), int(rec['end0']), str(rec['strand']))
            try:
                df_site = df_idx.loc[[key]].reset_index(drop=True)
            except KeyError:
                continue
            if df_site.empty:
                continue
            g_safe = re.sub(r'[^A-Za-z0-9._-]', '_', str(rec['gene']))
            c_safe = re.sub(r'[^A-Za-z0-9._-]', '_', str(rec['chrom']))
            title = (f"{rec['gene']} | {rec['mod']} | "
                     f"{rec['chrom']}:{rec['start0']}-{rec['end0']}({rec['strand']})\n"
                     f"{cat}   FDR={rec['padj']:.2e}, max|Δfrac|={rec['effect']:.3f}")
            os.makedirs(cat_dir, exist_ok=True)
            out_png = os.path.join(
                cat_dir,
                f"rank{rank:02d}__{g_safe}__{rec['mod']}__"
                f"{c_safe}_{rec['start0']}_{rec['end0']}_{rec['strand']}.png")
            try:
                make_plot(df_site, per_tx, title, out_png)
                n += 1
            except Exception as e:  # pragma: no cover - defensive
                print(f"[classify] WARN figure failed ({cat} {rec['gene']}): {e}",
                      file=sys.stderr)
        if n:
            made[cat] = n
    if verbose:
        tot = sum(made.values())
        print(f"[classify] rendered {tot} category figure(s) over {len(made)} "
              f"categor{'y' if len(made) == 1 else 'ies'} under {figs_dir}",
              file=sys.stderr)
    return made


# ----------------------- isoform architecture-map figures ----------------------- #

# Display names for the standard single-letter modkit codes; any other code
# (e.g. ChEBI ids like 17596) is shown verbatim so labels stay unambiguous.
MOD_DISP = {'a': 'm6A', 'm': '5mC', 'h': '5hmC', 'C': '4mC',
            '17596': 'inosine', '17802': 'pseudoU',
            '69426': 'Am', '19227': 'Um', '19228': 'Cm', '19229': 'Gm'}


def mod_label(code):
    return MOD_DISP.get(str(code), str(code))


def _zn_sort_key(z):
    s = str(z)
    return (0, int(s)) if s.lstrip('-').isdigit() else (1, s)


def plot_locus_arch(rec, iso, genes, out_png):
    """Render ONE isoform architecture-map (locus track) figure for a classified
    site -- the SAME plot style as the offline step7 analysis: every captured
    isoform of the gene is drawn as exon (blue=internal, orange=terminal/3'UTR) /
    intron (grey line) models; the site is marked on each isoform (green o =
    exonic in the mature RNA, red x = intronic/absent); each isoform is tagged
    HIGH/LOW/ANCHOR and labelled with its pooled stoichiometry (cov).

    Pooled per-isoform Nmod/Ncov come straight from the diff table's
    per_transcript_json, so no per-sample long table is required. Returns True if
    a figure was written.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plot_utils import save_figure
    from matplotlib.patches import Rectangle, Patch
    from matplotlib.lines import Line2D

    gene = rec['gene']; chrom = rec['chrom']; pos = int(rec['start0'])
    strand = rec['strand']; cat = rec['class_key']; mod = rec['mod']
    hiZN = str(rec.get('hiZN', '')); loZN = str(rec.get('loZN', ''))
    anchorZN = str(rec.get('anchorZN', '') or '')
    eff = float(rec['effect'])
    hf = float(rec.get('hi_frac', 0.0)); lf = float(rec.get('lo_frac', 0.0))
    mdisp = mod_label(mod)

    # pooled per-isoform [Nmod, Ncov] at the site, from per_transcript_json
    agg_site = {}
    for t in (rec.get('per_tx') or []):
        try:
            agg_site[str(t['ZN'])] = [int(t.get('Nmod', 0)), int(t.get('Ncov', 0))]
        except (KeyError, TypeError, ValueError):
            continue

    zns = sorted((z for z in genes.get(gene, []) if (gene, z) in iso), key=_zn_sort_key)
    if not zns:
        return False

    lo_g = min(iso[(gene, z)]['exons'][0][0] for z in zns)
    hi_g = max(iso[(gene, z)]['exons'][-1][1] for z in zns)
    span = max(hi_g - lo_g, 1)
    pad = max(span * 0.02, 200)
    x0, x1 = lo_g - pad, hi_g + pad

    fig, ax = plt.subplots(figsize=(16.5, 0.85 * len(zns) + 3.4), layout="constrained")
    yh = 0.6
    for row, z in enumerate(zns):
        d = iso[(gene, z)]
        y = len(zns) - row
        # intron line spanning the transcript body
        ax.plot([d['exons'][0][0], d['exons'][-1][1]], [y, y], color='0.6', lw=1, zorder=1)
        # exons: terminal (3'-most for the strand) orange, internal blue
        term_idx = len(d['exons']) - 1 if d['strand'] == '+' else 0
        for i, (s, e) in enumerate(d['exons']):
            col = '#d97706' if i == term_idx else '#2563eb'
            ax.add_patch(Rectangle((s, y - yh / 2), max(e - s, 1), yh,
                                   facecolor=col, edgecolor='none', zorder=2))
        # status of the site within THIS isoform (splice-aware)
        st = status_in(d['exons'], d['strand'], pos)
        st_short = {'exonic_terminal': 'exon-term', 'exonic_internal': 'exon-int',
                    'intronic': 'INTRONIC', 'absent': 'absent'}[st]
        nmod, cov = agg_site.get(str(z), [0, 0])
        intron_caveat = st in ('intronic', 'absent')
        if cov > 0:
            frac = nmod / cov
            lab = f"{mdisp}={frac:.2f} (cov={cov})"
            color = '#b91c1c' if intron_caveat else 'black'
            if intron_caveat:
                lab += "  [intron-read-derived, NOT mature isoform]"
        else:
            lab = "(not covered)"; color = '0.7'
        tag = []
        if z == hiZN: tag.append('HIGH')
        if z == loZN: tag.append('LOW')
        if z == anchorZN: tag.append('ANCHOR/longest')
        tagstr = (' [' + ','.join(tag) + ']') if tag else ''
        ax.text(x1 + pad * 0.3, y, f"ZN{z}/{d['arch']}{tagstr}  @site={st_short}   {lab}",
                va='center', ha='left', fontsize=8, color=color)
        # marker AT the site on this isoform row
        if st in ('exonic_terminal', 'exonic_internal'):
            ax.plot([pos], [y], marker='o', ms=6, mfc='#16a34a', mec='white', mew=0.8, zorder=5)
        else:
            ax.plot([pos], [y], marker='x', ms=7, color='#b91c1c', mew=2, zorder=5)

    # site marker + headline. Stop the red guide line just BELOW the red coordinate label (callout
    # above the line's tip). Offset the label to ONE SIDE of the line rather than centring it, so it
    # never crosses the red line or the axes frame: a site in the left half gets a right-extending
    # label (ha='left'), a site in the right half a left-extending one (ha='right').
    line_top = len(zns) + 0.55
    ax.plot([pos, pos], [0.2, line_top], color='red', ls='--', lw=1.2, zorder=3)
    x_right = x1 + span * 1.7 + pad          # matches the set_xlim() right edge below
    x_off = (x_right - x0) * 0.015
    if pos < 0.5 * (x0 + x_right):
        lab_x, lab_ha = pos + x_off, 'left'
    else:
        lab_x, lab_ha = pos - x_off, 'right'
    ax.text(lab_x, line_top + 0.12,
            f"{mdisp} @ {chrom}:{pos}\nΔ={eff:.2f} (hi {hf:.2f} vs lo {lf:.2f})",
            color='red', ha=lab_ha, va='bottom', fontsize=8)

    ax.set_xlim(x0, x1 + span * 1.7 + pad)   # extra right room for the (enlarged) per-ZN labels
    ax.set_ylim(0.2, len(zns) + 2.6)
    ax.set_yticks([])
    ax.set_xlabel(f"{chrom} genomic position ({strand} strand)")
    # Short headline (renders large under the house style); the colour/marker key is a real legend
    # below the axes (constrained_layout reserves its band) instead of a long, clip-prone title.
    ax.set_title(f"{cat}: {gene}  ({mdisp})", loc="left")
    _key = [
        Patch(facecolor='#2563eb', edgecolor='none', label='internal exon'),
        Patch(facecolor='#d97706', edgecolor='none', label="terminal exon / 3'UTR"),
        Line2D([0], [0], marker='o', color='none', markerfacecolor='#16a34a',
               markeredgecolor='white', markersize=8, label='site: exonic (mature)'),
        Line2D([0], [0], marker='x', color='#b91c1c', markersize=8, mew=2, lw=0,
               label='site: intronic / absent'),
    ]
    fig.legend(handles=_key, loc='outside lower center', ncol=4, frameon=False,
               fontsize=8, handletextpad=0.4, columnspacing=1.4)
    os.makedirs(os.path.dirname(out_png) or '.', exist_ok=True)
    save_figure(fig, out_png, dpi=140)   # PNG + SVG
    plt.close(fig)
    return True


def render_arch_figures(fig_records, iso, genes, figs_dir, per_category, verbose=False):
    """Render the per-category ISOFORM ARCHITECTURE-MAP figures (top
    ``per_category`` sites by effect of EACH category) under
    ``{figs_dir}/{CATEGORY}/rankNN__gene__mod__locus.png``. Built purely from the
    GTF isoform models + per_transcript_json (no --zn-long). Returns
    ``{category: n_figs_written}``.
    """
    if not fig_records or per_category <= 0:
        return {}
    try:
        import matplotlib  # noqa: F401  (fail fast if unavailable)
        matplotlib.use("Agg")
    except Exception as e:  # pragma: no cover - environment dependent
        print(f"[classify] WARN cannot import matplotlib for arch figures ({e}); skipping",
              file=sys.stderr)
        return {}
    by_cat = defaultdict(list)
    for rec in fig_records:
        by_cat[rec['class_key']].append(rec)
    made = {}
    for cat, recs in by_cat.items():
        recs = sorted(recs, key=lambda d: d['effect'], reverse=True)[:per_category]
        cat_dir = os.path.join(figs_dir, cat)
        n = 0
        for rank, rec in enumerate(recs, 1):
            g_safe = re.sub(r'[^A-Za-z0-9._-]', '_', str(rec['gene']))
            c_safe = re.sub(r'[^A-Za-z0-9._-]', '_', str(rec['chrom']))
            out_png = os.path.join(
                cat_dir,
                f"rank{rank:02d}__{g_safe}__{rec['mod']}__"
                f"{c_safe}_{rec['start0']}_{rec['end0']}_{rec['strand']}.png")
            try:
                if plot_locus_arch(rec, iso, genes, out_png):
                    n += 1
            except Exception as e:  # pragma: no cover - defensive
                print(f"[classify] WARN arch figure failed ({cat} {rec['gene']}): {e}",
                      file=sys.stderr)
        if n:
            made[cat] = n
    if verbose:
        tot = sum(made.values())
        print(f"[classify] rendered {tot} architecture-map figure(s) over {len(made)} "
              f"categor{'y' if len(made) == 1 else 'ies'} under {figs_dir}",
              file=sys.stderr)
    return made


def main():
    args = parse_args()
    iso, genes = load_isoforms(args.gtf, args.tes_tol, args.inside_tol)
    if args.verbose:
        print(f"[classify] isoform models: {len(iso)} over {len(genes)} genes", file=sys.stderr)

    mods = set(args.mod_filter or [])

    with open(args.diff_tsv) as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        diff_rows = list(rd)
    if args.verbose:
        print(f"[classify] diff sites read: {len(diff_rows)}", file=sys.stderr)

    out_cols = [
        'gene_name', 'mod_code', 'chrom', 'start0', 'end0', 'strand',
        'class_key', 'bucket', 'event', 'direction', 'structural_delta_nt',
        'n_tx_tested', 'effect_max_abs_frac_diff', 'p_adj_bh',
        'hi_ZN', 'hi_arch', 'hi_frac', 'lo_ZN', 'lo_arch', 'lo_frac', 'anchor_ZN',
        'status_hi', 'status_lo', 'jd_hi', 'jd_lo',
        'stoich_tier', 'hi_stoich_level',
    ]
    rows = []
    fig_records = []
    counts = Counter()
    n_considered = n_no_model = 0
    for r in diff_rows:
        mod = r.get('mod_code', '')
        if mods and mod not in mods:
            continue
        try:
            padj = float(r.get('p_adj_bh', 'nan'))
            eff = float(r.get('effect_max_abs_frac_diff', 'nan'))
        except ValueError:
            continue
        if not (padj <= args.fdr and eff >= args.min_effect):
            continue
        # Keep pos = 0-based bedMethyl start for the OUTPUT columns / figures (unchanged
        # semantics). Use cpos = start0+1 (== end0), the base's 1-based coordinate, ONLY for
        # comparison against the 1-based GTF exon coords in classify()/status_in()/junctions --
        # else every site is tested 1 bp upstream of its true base and boundary/EJC calls are
        # systematically off by one.
        try:
            gene = r['gene_name']; pos = int(r['start0']); cpos = pos + 1; strand = r['strand']
        except (ValueError, TypeError, KeyError):
            continue  # one malformed row must not abort the whole classification stage
        n_considered += 1

        try:
            per_tx = json.loads(r.get('per_transcript_json', '[]') or '[]')
        except (ValueError, TypeError):
            per_tx = []   # malformed JSON on a row must not abort the whole stage
        cov_tx = [t for t in per_tx if int(t.get('Ncov', 0)) >= args.min_cov
                  and (gene, str(t['ZN'])) in iso]
        if len(cov_tx) < 2:
            n_no_model += 1
            ck = 'UNEXPLAINABLE__NO_MODEL'
            counts[ck] += 1
            rows.append([gene, mod, r['chrom'], pos, r.get('end0', pos + 1), strand,
                         ck, 'UNEXPLAINABLE', 'NO_MODEL', '', '',
                         r.get('n_tx_tested', ''), f"{eff:.4f}", f"{padj:.3e}",
                         '', '', '', '', '', '', '',
                         '', '', '', '', '', ''])
            continue
        hi = max(cov_tx, key=lambda t: t['frac'])
        lo = min(cov_tx, key=lambda t: t['frac'])
        hiZN = str(hi['ZN']); loZN = str(lo['ZN'])
        anchorZN = anchor_of(gene, genes[gene], iso)
        bucket, event, direction, info = classify_tree(
            gene, cpos, hiZN, loZN, iso, tes_tol=args.tes_tol, ejc_nt=args.ejc_nt)
        stier = stoich_tier(hi['frac'] - lo['frac'])
        hlvl = hi_stoich_level(hi['frac'])
        class_key = f"{bucket}__{event}" + (f"__{direction}" if direction else "")
        counts[class_key] += 1
        fig_records.append({
            'class_key': class_key, 'bucket': bucket, 'event': event,
            'gene': gene, 'mod': mod, 'chrom': r['chrom'],
            'start0': pos, 'end0': int(r.get('end0', pos + 1)), 'strand': strand,
            'effect': eff, 'padj': padj, 'per_tx': per_tx,
            'hiZN': hiZN, 'loZN': loZN, 'anchorZN': anchorZN or '',
            'hi_frac': hi['frac'], 'lo_frac': lo['frac'],
        })
        jd_hi = info.get('jd_hi', ''); jd_lo = info.get('jd_lo', '')
        rows.append([
            gene, mod, r['chrom'], pos, r.get('end0', pos + 1), strand,
            class_key, bucket, event, direction, info.get('delta_nt', ''),
            r.get('n_tx_tested', ''), f"{eff:.4f}", f"{padj:.3e}",
            hiZN, iso[(gene, hiZN)]['arch'], f"{hi['frac']:.4f}",
            loZN, iso[(gene, loZN)]['arch'], f"{lo['frac']:.4f}",
            anchorZN or '',
            info.get('status_hi', ''), info.get('status_lo', ''),
            jd_hi if jd_hi != '' and jd_hi != 10**9 else '',
            jd_lo if jd_lo != '' and jd_lo != 10**9 else '',
            stier, hlvl,
        ])

    rows.sort(key=lambda x: (x[6], x[0], x[3]))
    import os
    os.makedirs(os.path.dirname(args.out_tsv) or '.', exist_ok=True)
    with open(args.out_tsv, 'w', newline='') as out:
        w = csv.writer(out, delimiter='\t')
        w.writerow(out_cols)
        w.writerows(rows)

    # --- coverage-independent PRIVATE-site scan (separate table; not gated by the differential test)
    if args.private_out_tsv and args.zn_long and os.path.exists(args.zn_long):
        priv_rows, priv_cols = scan_private_sites(
            args.zn_long, iso, genes,
            min_frac=args.private_min_frac, min_cov=args.private_min_cov)
        os.makedirs(os.path.dirname(args.private_out_tsv) or '.', exist_ok=True)
        with open(args.private_out_tsv, 'w', newline='') as out:
            w = csv.writer(out, delimiter='\t')
            w.writerow(priv_cols)
            w.writerows(priv_rows)
        from collections import Counter as _C
        pc = _C(r[7] for r in priv_rows)
        print(f"[ok] wrote {args.private_out_tsv}: {len(priv_rows)} PRIVATE sites "
              f"(coverage-independent scan)  " + " ".join(f"{k}={v}" for k, v in pc.items()))

    if args.figs_dir and args.zn_long:
        render_category_figures(fig_records, args.zn_long, args.figs_dir,
                                args.figs_per_category, verbose=args.verbose)
    if args.arch_figs_dir:
        render_arch_figures(fig_records, iso, genes, args.arch_figs_dir,
                            args.figs_per_category, verbose=args.verbose)

    tot = sum(counts.values())
    print(f"[ok] wrote {args.out_tsv}: {tot} classified sites "
          f"(mod={'/'.join(sorted(mods)) or 'ALL'}, FDR<={args.fdr}, effect>={args.min_effect})")
    if tot:
        # counts is keyed by the primary class_key (structural_category__stoich_direction);
        # print largest first.
        for key, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            if n:
                print(f"    {key:<44} {n:>6}  ({100*n/tot:.1f}%)")


if __name__ == '__main__':
    main()
