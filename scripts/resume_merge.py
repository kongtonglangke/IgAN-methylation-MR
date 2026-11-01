"""Resume the final partial chunk and merge all parts into the output file."""
import os, ssl, sys, time, urllib.request

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

# resume any incomplete part
for s, e, p in parts:
    want = e - s + 1
    have = os.path.getsize(p) if os.path.exists(p) else 0
    if have >= want:
        print("complete %s (%d/%d)" % (p, have, want), flush=True)
        continue
    print("resume %s from %d to %d" % (p, have, want), flush=True)
    rs = s + have
    for attempt in range(6):
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (rs, e)})
            with urllib.request.urlopen(req, timeout=300, context=CTX) as r, open(p, "ab") as f:
                while True:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    f.write(c)
            break
        except Exception as ex:
            print("  retry %d %s" % (attempt + 1, str(ex)[:50]), flush=True)
            time.sleep(5 * (attempt + 1))

# verify and merge
ok = True
for s, e, p in parts:
    want = e - s + 1
    have = os.path.getsize(p) if os.path.exists(p) else 0
    if have != want:
        print("SIZE MISMATCH %s %d != %d" % (p, have, want), flush=True)
        ok = False
if not ok:
    sys.exit(1)

with open(out, "wb") as fo:
    for s, e, p in parts:
        with open(p, "rb") as fi:
            while True:
                c = fi.read(1 << 20)
                if not c:
                    break
                fo.write(c)
print("merged size=%d expected=%d" % (os.path.getsize(out), total), flush=True)
if os.path.getsize(out) == total:
    for s, e, p in parts:
        os.remove(p)
    print("OK", flush=True)
else:
    print("MISMATCH", flush=True)
    sys.exit(1)
