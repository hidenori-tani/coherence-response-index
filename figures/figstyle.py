"""Shared matplotlib style used to draw the figures (journal print widths, TrueType Arial, minimum font sizes, Okabe-Ito colours). Copied from the author's lab-wide helper with comments and docstrings removed; the code is otherwise unchanged."""
from __future__ import annotations
import os
import re
import matplotlib
import matplotlib.pyplot as plt
MM_PER_INCH = 25.4
OKABE_ITO = {'black': '#000000', 'orange': '#E69F00', 'sky_blue': '#56B4E9', 'bluish_green': '#009E73', 'yellow': '#F0E442', 'blue': '#0072B2', 'vermillion': '#D55E00', 'reddish_purple': '#CC79A7'}
OKABE_ITO_CYCLE = [OKABE_ITO['blue'], OKABE_ITO['vermillion'], OKABE_ITO['bluish_green'], OKABE_ITO['orange'], OKABE_ITO['sky_blue'], OKABE_ITO['reddish_purple'], OKABE_ITO['black'], OKABE_ITO['yellow']]
CVD_SAFE_6 = {'red': '#EE6677', 'green': '#228833', 'cyan': '#66CCEE', 'purple': '#AA3377', 'indigo': '#332288', 'orange': '#E69F00'}
CVD_SAFE_6_LIST = list(CVD_SAFE_6.values())
SEQ3 = ['#F2F2F2', '#7FB2D9', '#08306B']
SEQ3_INK = ['#222222', '#0B2E4F', '#FFFFFF']

def ink_on(hex_color: str) -> str:
    h = hex_color.lstrip('#')
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))
    return '#FFFFFF' if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.45 else '#111111'
WIDTHS_MM = {'nature': {'single': 89.0, 'double': 183.0}, 'elsevier': {'single': 90.0, 'onehalf': 140.0, 'double': 190.0}, 'cellpress': {'single': 85.0, 'onehalf': 114.0, 'double': 174.0}, 'plos': {'single': 66.8, 'double': 190.5}, 'springer': {'single': 84.0, 'double': 174.0}, 'jbs': {'single': 85.0, 'double': 170.0}, 'acs': {'single': 84.7, 'double': 177.8}, 'tandf': {'single': 89.6, 'double': 182.8}, 'wiley': {'single': 80.0, 'double': 180.0}, 'clinchem': {'single': 85.0, 'double': 176.0}, 'cclm': {'single': 85.0, 'double': 176.0}, 'generic': {'single': 85.0, 'double': 174.0}}
BASE_FONT_PT = {'nature': 7.0, 'elsevier': 8.0, 'cellpress': 7.5, 'plos': 9.0, 'springer': 9.0, 'acs': 7.0, 'wiley': 8.0, 'clinchem': 8.5, 'cclm': 8.5, 'generic': 8.0, 'jbs': 8.8}
JOURNAL_MIN_PT = {'nature': 5.0, 'elsevier': 7.0, 'cellpress': 6.0, 'plos': 8.0, 'springer': 7.0, 'acs': 4.5, 'tandf': 7.0, 'generic': 7.0, 'clinchem': 8.0, 'cclm': 6.0, 'jbs': 8.0}

def safe_min_pt(journal: str='generic') -> float:
    return round(JOURNAL_MIN_PT.get(journal, 7.0) * 1.115, 1)

def pt_scale(journal: str='generic') -> dict:
    lo = safe_min_pt(journal)
    return {'tiny': lo, 'small': lo + 0.4, 'body': lo + 0.8, 'label': lo + 1.4, 'head': lo + 2.2, 'big': lo + 4.2}

def mm(value_mm: float) -> float:
    return value_mm / MM_PER_INCH

def _pick_sans() -> list:
    available = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
    order = ['Arial', 'Helvetica', 'Helvetica Neue', 'Liberation Sans', 'Nimbus Sans']
    picked = [name for name in order if name in available]
    return picked + ['DejaVu Sans']

def use_journal_style(journal: str='generic', base_pt: float | None=None) -> None:
    pt = base_pt if base_pt else BASE_FONT_PT.get(journal, 8.0)
    sans = _pick_sans()
    head = sans[0]
    plt.rcParams.update({'mathtext.fontset': 'custom', 'mathtext.rm': head, 'mathtext.it': f'{head}:italic', 'mathtext.bf': f'{head}:bold', 'mathtext.sf': head, 'mathtext.tt': head, 'mathtext.cal': f'{head}:italic', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'font.family': 'sans-serif', 'font.sans-serif': sans, 'font.size': pt, 'axes.titlesize': pt, 'axes.labelsize': pt, 'xtick.labelsize': max(pt - 0.5, safe_min_pt(journal)), 'ytick.labelsize': max(pt - 0.5, safe_min_pt(journal)), 'legend.fontsize': max(pt - 0.5, safe_min_pt(journal)), 'figure.titlesize': pt + 1, 'axes.linewidth': 0.6, 'grid.linewidth': 0.5, 'lines.linewidth': 1.0, 'lines.markersize': 3.0, 'patch.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.5, 'ytick.minor.width': 0.5, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.prop_cycle': plt.cycler(color=OKABE_ITO_CYCLE), 'legend.frameon': False, 'figure.dpi': 150, 'savefig.dpi': 600, 'savefig.bbox': 'standard', 'savefig.transparent': False, 'figure.facecolor': 'white', 'savefig.facecolor': 'white'})

def new_figure(journal: str='generic', column: str='single', height_mm: float=60.0, **kwargs):
    width_mm = WIDTHS_MM.get(journal, WIDTHS_MM['generic'])[column]
    return plt.subplots(figsize=(mm(width_mm), mm(height_mm)), **kwargs)
EM_STAMP_Y_MM = (1.49, 8.55)
EM_LABEL_MM = (1.76, 1.49, 30.4, 5.02)

def warn_if_em_stamp_collides(pdf_path: str) -> bool:
    try:
        import fitz
    except ImportError:
        return False
    to_pt = lambda v: v / 25.4 * 72
    page = fitz.open(pdf_path)[0]
    w, h = (page.rect.width, page.rect.height)
    stamp = fitz.Rect(w / 2, to_pt(EM_STAMP_Y_MM[0]), w, to_pt(EM_STAMP_Y_MM[1]))
    label = fitz.Rect(*(to_pt(v) for v in EM_LABEL_MM))

    def touches(r):
        if r.height <= 0 or r.width <= 0:
            r = fitz.Rect(r.x0 - 0.05, r.y0 - 0.05, r.x1 + 0.05, r.y1 + 0.05)
        return r.intersects(stamp) or r.intersects(label)
    hit = []
    for b in page.get_text('dict')['blocks']:
        for line in b.get('lines', []):
            for s in line.get('spans', []):
                if s.get('text', '').strip() and touches(fitz.Rect(s['bbox'])):
                    hit.append(s['text'].strip()[:24])
    for it in page.get_drawings():
        r = it['rect']
        if r.width >= w * 0.99 and r.height >= h * 0.99:
            continue
        if touches(r):
            hit.append('(図形)')
    if hit:
        print(f"  ⚠️ {os.path.basename(pdf_path)}: EM のスタンプ域に {len(hit)} 件（{'／'.join((repr(t) for t in hit[:3]))}）。Elsevier/Cell Press へ出すなら、上端 y {EM_STAMP_Y_MM[1]} mm・ページ幅の中央から右を空けること（FIGURE_GUIDE.md §6.4）")
    return bool(hit)

def warn_if_text_offpage(fig) -> list:
    fig.canvas.draw()
    W, H = (fig.bbox.x1, fig.bbox.y1)
    ignore = set()
    for ax in fig.get_axes():
        for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            lo, hi = (min(lim), max(lim))
            for tick in list(axis.get_major_ticks()) + list(axis.get_minor_ticks()):
                if not lo - 1e-09 <= tick.get_loc() <= hi + 1e-09:
                    ignore.add(id(tick.label1))
                    ignore.add(id(tick.label2))
    lost = []
    for txt in fig.findobj(match=lambda o: hasattr(o, 'get_text')):
        if id(txt) in ignore:
            continue
        s = (txt.get_text() or '').strip()
        if not s or not txt.get_visible():
            continue
        try:
            bb = txt.get_window_extent()
        except Exception:
            continue
        if bb.x1 < 0 or bb.y1 < 0 or bb.x0 > W or (bb.y0 > H):
            lost.append(s[:28])
    if lost:
        print(f"  🚨 紙の外に落ちた文字 {len(lost)} 件（{'／'.join((repr(t) for t in lost[:4]))}）＝出力ファイルには存在しない。余白か図の高さを増やすこと")
    return lost

def warn_if_clipped_by_axes(fig, tol_pt: float=0.05) -> list:
    fig.canvas.draw()
    px = fig.dpi / 72.0
    tol = tol_pt * px
    lost = []
    for ax in fig.get_axes():
        cands = list(ax.patches)
        for t in ax.texts:
            ap = getattr(t, 'arrow_patch', None)
            if ap is not None:
                cands.append(ap)
        for art in cands:
            if art is ax.patch or not art.get_visible() or (not art.get_clip_on()):
                continue
            cb = art.get_clip_box() or ax.bbox
            try:
                bb = art.get_window_extent()
            except Exception:
                continue
            pad = 0.0
            try:
                ec = art.get_edgecolor()
                if art.get_linewidth() > 0 and (len(ec) < 4 or ec[3] > 0):
                    pad = art.get_linewidth() * px / 2.0
            except Exception:
                pass
            over = {'左': cb.x0 - (bb.x0 - pad), '右': bb.x1 + pad - cb.x1, '下': cb.y0 - (bb.y0 - pad), '上': bb.y1 + pad - cb.y1}
            hit = {k: v for k, v in over.items() if v > tol}
            if hit:
                where = '・'.join((f'{k} {v / px / 72 * 25.4:.2f}mm' for k, v in hit.items()))
                lost.append(f'{type(art).__name__}({where})')
    if lost:
        print(f"  🚨 軸の枠で切られた図形 {len(lost)} 件（{'／'.join(lost[:4])}）＝出力では枠線や角が欠ける。座標を軸の内側に収めること")
    return lost

def figure_label(fig, text: str, width_mm: float, height_mm: float, left_mm: float=2.0, top_mm: float=6.6, fontsize: float=12.0) -> None:
    fig.text(left_mm / width_mm, 1.0 - top_mm / height_mm, text, fontsize=fontsize, fontweight='bold', ha='left', va='top', color='#111111')
_CMAP_CACHE: dict = {}

def _font_cmap(fp) -> set:
    path = matplotlib.font_manager.findfont(fp)
    if path not in _CMAP_CACHE:
        try:
            from fontTools.ttLib import TTFont
            f = TTFont(path, fontNumber=0, lazy=True)
            cps = set()
            for t in f['cmap'].tables:
                cps |= set(t.cmap.keys())
            _CMAP_CACHE[path] = cps
        except Exception:
            _CMAP_CACHE[path] = None
    return _CMAP_CACHE[path]

def warn_if_missing_glyph(fig) -> list:
    from matplotlib.text import Text
    missing = []
    for txt in fig.findobj(match=lambda o: isinstance(o, Text)):
        s = txt.get_text() or ''
        if not s.strip() or not txt.get_visible():
            continue
        fp = txt.get_fontproperties()
        cmap = _font_cmap(fp)
        if cmap is None:
            continue
        for ch in set(s):
            if ch in '\n\t ':
                continue
            if ord(ch) not in cmap:
                missing.append((ch, s.strip()[:26]))
    if missing:
        uniq = sorted({c for c, _ in missing})
        print(f"  🚨 書体に無い文字 {len(missing)} 箇所（{' '.join(uniq)}）＝PDF には代替字形（豆腐）が入る。例: {missing[0][1]!r}")
    return missing

def warn_if_foreign_font(pdf_path: str) -> list:
    try:
        import fitz
    except Exception:
        return []
    allowed = {re.sub('[^a-z]', '', n.lower()) for n in _pick_sans()}
    foreign = []
    with fitz.open(pdf_path) as doc:
        for page in doc:
            for f in page.get_fonts():
                base = re.sub('^[A-Z]{6}\\+', '', f[3])
                key = re.sub('[^a-z]', '', base.lower())
                if not any((key.startswith(a) for a in allowed)):
                    foreign.append(base)
    foreign = sorted(set(foreign))
    if foreign:
        print(f'  🚨 {os.path.basename(pdf_path)}: 指定外の書体が埋め込まれている {foreign}＝mathtext の記号（$\\sim$ / $\\approx$ 等）が原因のことが多い。Arial にある文字を直接書く')
    return foreign

def savefig(fig, stem: str, dpi: int=600, formats=('pdf', 'png')) -> list:
    os.makedirs(os.path.dirname(stem) or '.', exist_ok=True)
    warn_if_text_offpage(fig)
    warn_if_missing_glyph(fig)
    written = []
    for ext in formats:
        path = f'{stem}.{ext}'
        fig.savefig(path, format=ext, dpi=dpi)
        written.append(path)
        if ext == 'pdf':
            warn_if_em_stamp_collides(path)
            warn_if_foreign_font(path)
    return written

def cvd_safe_pair(name_a: str='blue', name_b: str='vermillion') -> tuple:
    return (OKABE_ITO[name_a], OKABE_ITO[name_b])

def _selftest(outdir: str='/tmp/figstyle_selftest') -> str:
    import numpy as np
    W, H = (WIDTHS_MM['elsevier']['single'], 58.0)
    use_journal_style('elsevier')
    fig = plt.figure(figsize=(mm(W), mm(H)))
    left, right, top, bottom = (13.0, 2.0, 10.5, 11.0)
    ax = fig.add_axes([left / W, bottom / H, (W - left - right) / W, (H - top - bottom) / H])
    x = np.linspace(0, 10, 100)
    for i, (label, shift) in enumerate([('control', 0.0), ('treated', 0.6)]):
        ax.plot(x, np.sin(x) + shift, label=label, linestyle=['-', '--'][i])
    ax.set_xlabel('time (h)')
    ax.set_ylabel('relative level')
    ax.legend(loc='upper right')
    stem = os.path.join(outdir, 'selftest')
    paths = savefig(fig, stem)
    plt.close(fig)
    return paths[0]

def _selftest_against_gate(pdf_path: str) -> int:
    import subprocess
    import sys as _sys
    here = os.path.dirname(os.path.abspath(__file__))
    r = subprocess.run([_sys.executable, os.path.join(here, 'check_figure.py'), pdf_path, '--journal', 'elsevier', '--column', 'single'], capture_output=True, text=True)
    m = re.search('合計: FAIL (\\d+) 件 / WARN (\\d+) 件', r.stdout)
    if not m:
        print('  ⚠️ check_figure.py の出力を読み取れなかった（判定していない）')
        return 1
    n_fail, n_warn = (int(m.group(1)), int(m.group(2)))
    print(f'  見本図を判定器に当てた結果: FAIL {n_fail} 件 / WARN {n_warn} 件')
    if n_fail:
        for line in r.stdout.splitlines():
            if '[FAIL]' in line:
                print('   ', line.strip())
    return 1 if n_fail else 0

def selftest_palettes() -> int:
    import itertools
    import sys as _sys
    import numpy as np
    _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from check_figure import simulate_cvd, srgb_to_lab, delta_e76, relative_luminance

    def _rgb(h):
        h = h.lstrip('#')
        return [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]

    def worst(hexes):
        cols = np.array([_rgb(h) for h in hexes])
        w = 1000000000.0
        for kind in ('protanopia', 'deuteranopia', 'tritanopia'):
            lab = srgb_to_lab(simulate_cvd(cols, kind))
            for i, j in itertools.combinations(range(len(cols)), 2):
                w = min(w, delta_e76(lab[i], lab[j]))
        return w

    def gap(hexes):
        lums = sorted((relative_luminance(np.array(_rgb(h))) for h in hexes))
        return min((b - a for a, b in zip(lums, lums[1:])))
    fails, checks = ([], 0)
    checks += 1
    if worst(CVD_SAFE_6_LIST) < 12:
        fails.append(f'CVD_SAFE_6 の最小ΔE {worst(CVD_SAFE_6_LIST):.1f} < 12')
    checks += 1
    if gap(CVD_SAFE_6_LIST) < 0.05:
        fails.append(f'CVD_SAFE_6 の灰階差 {gap(CVD_SAFE_6_LIST):.3f} < 0.05')
    checks += 1
    if worst(SEQ3) < 12:
        fails.append(f'SEQ3 の最小ΔE {worst(SEQ3):.1f} < 12')
    checks += 1
    okabe6 = [OKABE_ITO[k] for k in ('orange', 'sky_blue', 'bluish_green', 'yellow', 'blue', 'vermillion')]
    if worst(okabe6) >= 12:
        fails.append(f'陰性対照が落ちない＝色覚の判定が働いていない（Okabe-Ito 6色 minΔE={worst(okabe6):.1f}）')
    checks += 1
    if worst(okabe6[:5]) < 8:
        fails.append(f'Okabe-Ito 5色が FAIL 域＝判定が厳しすぎる（minΔE={worst(okabe6[:5]):.1f}）')
    checks += 1
    if worst(okabe6[:5]) >= 12:
        fails.append(f'Okabe-Ito 5色が PASS 域に入った＝実測(11.6)と違う。判定器が変わった可能性があるので CVD_SAFE_6 の根拠を測り直すこと')
    checks += 1
    for j in ('elsevier', 'nature', 'springer', 'generic'):
        sc = pt_scale(j)
        if min(sc.values()) < safe_min_pt(j):
            fails.append(f'{j}: pt_scale の最小 {min(sc.values())} < {safe_min_pt(j)}')
    checks += 1
    if JOURNAL_MIN_PT['elsevier'] >= safe_min_pt('elsevier'):
        fails.append('safe_min_pt が下限と同じ＝WARN 域を避けられていない')
    checks += 1
    if ink_on('#332288') != '#FFFFFF' or ink_on('#F0E442') != '#111111':
        fails.append('ink_on の白黒判定が反転している')
    checks += 1
    if 'Arial' not in {f.name for f in matplotlib.font_manager.fontManager.ttflist}:
        fails.append('Arial が見つからない＝DejaVu に落ちる')
    print(f'  CVD_SAFE_6 minΔE={worst(CVD_SAFE_6_LIST):5.1f} 灰階差={gap(CVD_SAFE_6_LIST):.3f}')
    print(f'  SEQ3       minΔE={worst(SEQ3):5.1f} 灰階差={gap(SEQ3):.3f}')
    print(f'  陰性対照 Okabe-Ito 6色 minΔE={worst(okabe6):5.1f}（12 未満で正常）')
    print(f'  対照 Okabe-Ito 5色     minΔE={worst(okabe6[:5]):5.1f}（8 以上 12 未満で正常＝WARN 域）')
    for f in fails:
        print('  FAIL', f)
    checks += 1
    fig = plt.figure(figsize=(mm(60), mm(40)))
    ax = fig.add_axes([0.02, 0.02, 0.96, 0.96])
    ax.set_ylabel('a very long axis label that cannot fit in this margin at all')
    if not warn_if_text_offpage(fig):
        fails.append('陽性対照が拾われない＝warn_if_text_offpage が何も見ていない')
    plt.close(fig)
    checks += 1
    fig = plt.figure(figsize=(mm(80), mm(50)))
    ax = fig.add_axes([0.18, 0.18, 0.78, 0.72])
    ax.plot([0, 10], [0, 1])
    ax.set_xlim(0, 10)
    if warn_if_text_offpage(fig):
        fails.append('表示範囲外の目盛ラベルを誤検出している')
    plt.close(fig)
    use_journal_style('elsevier')
    checks += 1
    fig = plt.figure(figsize=(mm(60), mm(30)))
    fig.text(0.5, 0.5, '★ star', fontsize=10)
    if not warn_if_missing_glyph(fig):
        fails.append('陽性対照が拾われない＝warn_if_missing_glyph が何も見ていない')
    plt.close(fig)
    checks += 1
    fig = plt.figure(figsize=(mm(60), mm(30)))
    fig.text(0.5, 0.5, '• → ≥ – ok', fontsize=10)
    if warn_if_missing_glyph(fig):
        fails.append('在る文字を『無い』と言っている＝誤検出')
    plt.close(fig)
    checks += 1
    from matplotlib.patches import Rectangle as _Rect
    fig = plt.figure(figsize=(mm(60), mm(30)))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.add_patch(_Rect((90, 2), 11, 6, fc='none', ec='0.3', lw=0.9))
    if not warn_if_clipped_by_axes(fig):
        fails.append('陽性対照が拾われない＝warn_if_clipped_by_axes が何も見ていない')
    plt.close(fig)
    checks += 1
    fig = plt.figure(figsize=(mm(60), mm(30)))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 10)
    ax.add_patch(_Rect((10, 2), 20, 6, fc='none', ec='0.3', lw=0.9))
    ax.barh([5], [40], left=0, height=1.0, color='0.6')
    ax.plot([-50, 150], [8, 8], lw=0.8)
    if warn_if_clipped_by_axes(fig):
        fails.append('切れていないものを『切れた』と言っている＝誤検出')
    plt.close(fig)
    print(f'配色・文字サイズの自己試験: {checks - len(fails)}/{checks}')
    return 1 if fails else 0
if __name__ == '__main__':
    rc = selftest_palettes()
    pdf = _selftest()
    print('wrote:', pdf)
    rc |= _selftest_against_gate(pdf)
    print(f"figstyle selftest: {('PASS' if rc == 0 else 'FAIL')}")
    raise SystemExit(rc)
