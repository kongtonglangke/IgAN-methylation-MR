import json, urllib.request, ssl, time, os
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE
genes = ["TNFSF13","TNFRSF13B","CFB","CFH","CFHR1","CFHR3","C1GALT1","ITGAM","ITGAX",
         "DEFA1","HORMAD2","LIF","OSM","VAV3","PXK","SPRY2","ST6GAL1","AICDA",
         "HDAC1","HDAC2","HDAC3","ITGA4","ITGB7","CCR9","MADCAM1","AHR","IDO1",
         "TDO2","NR1H4","GPBAR1","FFAR2","FFAR3","SLC16A1","SLC5A8","TNFSF13B","TNFRSF17"]
url = "https://grch37.rest.ensembl.org/lookup/symbol/homo_sapiens"
body = json.dumps({"symbols": genes}).encode()
req = urllib.request.Request(url, data=body, headers={"Content-Type":"application/json","Accept":"application/json"})
try:
    with urllib.request.urlopen(req, timeout=60, context=ctx) as r:
        data = json.loads(r.read().decode())
    rows=[]
    for g,v in data.items():
        if v is None:
            rows.append((g,"NA",0,0,0)); continue
        rows.append((g, v["seq_region_name"], v["start"], v["end"], v.get("strand",1)))
    os.makedirs("data", exist_ok=True)
    with open("data/gene_coords_grch37.tsv","w",encoding="utf-8") as f:
        f.write("gene\tchr\tstart\tend\tstrand\n")
        for x in rows:
            f.write("\t".join(str(i) for i in x)+"\n")
    ok=sum(1 for x in rows if x[1]!="NA")
    print("OK resolved=%d/%d"%(ok,len(rows)))
    for x in rows: print("\t".join(str(i) for i in x))
except Exception as e:
    print("ERR", type(e).__name__, e)
