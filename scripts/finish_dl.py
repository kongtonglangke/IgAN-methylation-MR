"""Append the missing tail of specific parts, then merge and verify gzip integrity."""
import gzip, os, ssl, subprocess, sys, time, urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

url = sys.argv[1]
out = sys.argv[2]
total = int(sys.argv[3])
nthread = 8
step = total // nthread + 1

parts = []
for i in range(nthread):
    s = i * step
    e = min(total - 1, s + step - 1)
    if s > e:
        break
    parts.append((s, e, out + ".part%d" % i))


def fill(s, e, p):
    want = e - s + 1
    have = os.path.getsize(p) if os.path.exists(p) else 0
    while have < want:
        rs = s + have
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (rs, e)})
            with urllib.request.urlopen(req, timeout=300, context=CTX) as r, open(p, "ab") as f:
                while True:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    if have + len(c) > want:
                        f.write(c[: want - have])
                        have = want
                        break
                    f.write(c)
                    have += len(c)
        except Exception as ex:
            print("  retry %s %s" % (p, str(ex)[:50]), flush=True)
            time.sleep(8)
        if os.path.getsize(p) >= want:
            break
    return os.path.getsize(p) == want


allok = True
for s, e, p in parts:
    want = e - s + 1
    have = os.path.getsize(p) if os.path.exists(p) else 0
    if have == want:
        print("ok   %s" % os.path.basename(p), flush=True)
        continue
    print("fill %s have=%d want=%d" % (os.path.basename(p), have, want), flush=True)
    if not fill(s, e, p):
        print("FAILED %s" % p, flush=True)
        allok = False
if not allok:
    sys.exit(1)

tmp = out + ".merging"
with open(tmp, "wb") as fo:
    for s, e, p in parts:
        with open(p, "rb") as fi:
            while True:
                c = fi.read(1 << 20)
                if not c:
                    break
                fo.write(c)
print("merged=%d expected=%d" % (os.path.getsize(tmp), total), flush=True)
if os.path.getsize(tmp) != total:
    sys.exit(1)
os.replace(tmp, out)

rc = subprocess.run(["gzip", "-t", out], capture_output=True)
print("GZIP_OK=%s" % (rc.returncode == 0), flush=True)
if rc.returncode != 0:
    print(rc.stderr.decode("utf-8", "ignore")[:200], flush=True)
    sys.exit(1)
with gzip.open(out, "rt", errors="ignore") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        print("HEADER:", line.rstrip()[:300], flush=True)
        break
for s, e, p in parts:
    os.remove(p)
print("DONE", flush=True)
