import urllib.request, ssl, json, re, time
from concurrent.futures import ThreadPoolExecutor
ctx=ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE

def bucket(n):
    if 90001001<=n<=90002000: return "GCST90001001-GCST90002000"
    if 90002001<=n<=90003000: return "GCST90002001-GCST90003000"
    return None

def get(n):
    g="GCST%d"%n
    b=bucket(n)
    if not b: return (g,"")
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/%s/%s/%s_buildGRCh37.tsv.gz-meta.yaml"%(b,g,g)
    for a in range(2):
        try:
            with urllib.request.urlopen(url,timeout=30,context=ctx) as r:
                t=r.read().decode("utf-8","ignore")
            m=re.search(r"trait_description:\n\s*-\s*(.+)",t)
            return (g, m.group(1).strip() if m else "?")
        except Exception as e:
            if a==1: return (g,"")
            time.sleep(1)

res=[]
with ThreadPoolExecutor(max_workers=12) as ex:
    for g,t in ex.map(get, range(90001391,90002122)):
        if t: res.append((g,t))
print("resolved=%d"%len(res),flush=True)
out=open("data/orru_traits.tsv","w")
out.write("gwas_id\ttrait\n")
for g,t in res: out.write("%s\t%s\n"%(g,t))
out.close()
bc=[(g,t) for g,t in res if re.search(r"\bB cell\b",t) and not re.search(r"T cell|NK|monocyte|neutrophil|platelet|RBC|basophil|eosinophil",t)]
print("=== B cell traits: %d ==="%len(bc))
for g,t in bc: print(g,t)
