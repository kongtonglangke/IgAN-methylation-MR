"""Re-download only the parts whose size does not match, then merge and verify gzip."""
import gzip, os, ssl, sys, time, urllib.request

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


def dl(s, e, p):
    if os.path.exists(p):
        os.remove(p)
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (s, e)})
            with urllib.request.urlopen(req, timeout=300, context=CTX) as r, open(p, "wb") as f:
                while True:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    f.write(c)
            if os.path.getsize(p) == e - s + 1:
                return True
        except Exception as ex:
            print("  retry %d %s" % (attempt + 1, str(ex)[:50]), flush=True)
        time.sleep(5 * (attempt + 1))
    return False


for s, e, p in parts:
    want = e - s + 1
    have = os.path.getsize(p) if os.path.exists(p) else 0
    if have == want:
        print("ok   %s %d" % (p, have), flush=True)
        continue
    print("fix  %s have=%d want=%d" % (p, have, want), flush=True)
    if not dl(s, e, p):
        print("FAILED %s" % p, flush=True)
        sys.exit(1)

with open(out, "wb") as fo:
    for s, e, p in parts:
        with open(p, "rb") as fi:
            while True:
                c = fi.read(1 << 20)
                if not c:
                    break
                fo.write(c)
print("merged=%d expected=%d" % (os.path.getsize(out), total), flush=True)
# verify gzip integrity + peek header
try:
    op = gzip.open(out, "rt", errors="ignore")
    for i, line in enumerate(op):
        if line.startswith("#"):
            continue
        print("HEADER:", line.rstrip()[:300], flush=True)
        break
    op.close()
    import subprocess
    print("VERIFY:", subprocess.run(["gzip", "-t", out], capture_output=True).returncode == 0, flush=True)
except Exception as ex:
    print("GZIP BROKEN", str(ex)[:100], flush=True)
    sys.exit(1)
for s, e, p in parts:
    os.remove(p)
print("OK", flush=True)
