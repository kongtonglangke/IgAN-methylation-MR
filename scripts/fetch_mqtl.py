import json, urllib.request, ssl, os, time
from concurrent.futures import ThreadPoolExecutor
ctx = ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE
WIN = 500000

genes=[]
with open("data/gene_coords_grch37.tsv") as f:
    next(f)
    for line in f:
        p=line.rstrip("\n").split("\t")
        if p[1]=="NA": continue
        genes.append((p[0],p[1],int(p[2]),int(p[3])))

def fetch(g):
    name,chr_,s,e = g
    lo=max(1,s-WIN); hi=e+WIN
    # X/Y 跳过（GoDMC 主要常染色体）
    if chr_ in ("X","Y","MT"): return (name,None,"skip")
    url="http://api.godmc.org.uk/v0.1/assoc_meta/range/cpg/%s:%d-%d"%(chr_,lo,hi)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=180, context=ctx) as r:
                d=json.loads(r.read().decode())
            rows=d.get("assoc_meta",[])
            return (name, rows, "ok")
        except Exception as ex:
            if attempt==2: return (name,None,"ERR:"+str(ex)[:80])
            time.sleep(2)

allrows=[]
meta=[]
with ThreadPoolExecutor(max_workers=6) as ex:
    for name,rows,st in ex.map(fetch, genes):
        if rows is None:
            meta.append((name,0,st)); print("FAIL",name,st,flush=True); continue
        for r in rows: r["_gene"]=name
        allrows += rows
        meta.append((name,len(rows),st))
        print("OK %-12s rows=%d"%(name,len(rows)),flush=True)

os.makedirs("data",exist_ok=True)
with open("data/mqtl_raw.json","w") as f: json.dump(allrows,f)
cpgs=sorted({r["cpg"] for r in allrows})
print("TOTAL rows=%d distinct_cpg=%d"%(len(allrows),len(cpgs)))
with open("data/cpg_list.txt","w") as f: f.write("\n".join(cpgs))
