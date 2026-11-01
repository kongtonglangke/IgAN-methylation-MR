# -*- coding: utf-8 -*-
"""G4: ATAC peak 的细胞类群归属与最近距离

由 ChIP-Atlas assembled AllCell BED 的 name 字段解析 Cell group / cell line，
区分血液来源，并计算与目标 CpG 的最近峰距离。
"""
import collections
import csv
import io
import os
import re
import urllib.parse

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")

REGIONS = [
    ("DEFA1_8p23", "cg00705004", 6970292),
    ("HORMAD2_22q12", "cg00539347", 30196307),
]

# 常见血液/免疫相关细胞系的谱系归类
LINEAGE = {
    "B": ["gm12878", "raji", "daudi", "nalm6", "namalwa", "bjab", "ocily", "su-dhl",
          "rpmi", "jurkat-b", "b-cell", "b cell", "lymphoblast"],
    "T": ["jurkat", "cem", "molt", "hsb2", "t-cell", "t cell", "cd4", "cd8"],
    "Myeloid": ["thp-1", "thp1", "hl-60", "hl60", "u937", "k562", "monocyte",
                "macrophage", "neutrophil", "granulocyte", "mast"],
    "HSC/Progenitor": ["cd34", "hspc", "hematopoietic"],
}


def parse_name(name):
    s = urllib.parse.unquote(name)
    d = {}
    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
        elif kv.startswith("ID="):
            d["ID"] = kv[3:]
    return d, s.lower()


def lineage_of(low):
    for lin, keys in LINEAGE.items():
        for k in keys:
            if k in low:
                return lin
    return "other"


def main():
    rows = list(csv.DictReader(io.open(os.path.join(RES, "g4_atac_peaks.tsv"),
                                        encoding="utf-8"), delimiter="\t"))
    out = []
    for region, cpg, pos in REGIONS:
        sub = [r for r in rows if r["region"] == region]
        # 去重：同一区间保留任一行作代表
        seen = {}
        for r in sub:
            key = (r["start"], r["end"])
            if key not in seen:
                meta, low = parse_name(r["name"])
                seen[key] = {
                    "start": int(r["start"]), "end": int(r["end"]),
                    "group": meta.get("Cell group", "NA"),
                    "line": meta.get("cell line", meta.get("cell type", "NA")),
                    "lin": lineage_of(low),
                    "inside": r["cpg_inside"] == "YES",
                }
        uniq = list(seen.values())
        blood = [u for u in uniq if u["group"].lower().startswith("blood")]
        inside = [u for u in uniq if u["inside"]]
        inside_blood = [u for u in inside if u["group"].lower().startswith("blood")]

        def nearest(pool):
            if not pool:
                return None
            return min(min(abs(pos - u["start"]), abs(pos - u["end"])) for u in pool)

        lin_blood = collections.Counter(u["lin"] for u in blood)
        lin_inside = collections.Counter(u["lin"] for u in inside)

        out.append((region, cpg, pos, len(uniq), len(blood), len(inside),
                    len(inside_blood), nearest(uniq), nearest(blood),
                    lin_blood.get("B", 0), lin_blood.get("T", 0),
                    lin_blood.get("Myeloid", 0),
                    lin_inside.get("B", 0), lin_inside.get("T", 0),
                    lin_inside.get("Myeloid", 0)))

        print("== %s (%s @ %d)" % (region, cpg, pos))
        print("   唯一区间 %d | Blood 类 %d (B %d / T %d / Myeloid %d)" %
              (len(uniq), len(blood), lin_blood.get("B", 0), lin_blood.get("T", 0),
               lin_blood.get("Myeloid", 0)))
        print("   覆盖 CpG: 全部 %d | Blood %d" % (len(inside), len(inside_blood)))
        print("   覆盖 CpG 的谱系: %s" % dict(lin_inside))
        d_all, d_blood = nearest(uniq), nearest(blood)
        print("   最近峰距离: 全部 %s bp | Blood %s bp" %
              (d_all, d_blood if d_blood is not None else "无"))
        print()

    outp = os.path.join(RES, "g4_atac_summary.tsv")
    with io.open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(["region", "cpg", "cpg_pos", "n_unique_peak", "n_blood_peak",
                            "n_peak_covering_cpg", "n_blood_peak_covering_cpg",
                            "min_dist_all_bp", "min_dist_blood_bp",
                            "blood_B", "blood_T", "blood_myeloid",
                            "inside_B", "inside_T", "inside_myeloid"]) + "\n")
        for r in out:
            fh.write("\t".join(str(v) for v in r) + "\n")
    print("写出", outp)


if __name__ == "__main__":
    main()
