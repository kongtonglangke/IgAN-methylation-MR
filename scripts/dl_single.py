"""Single-stream resumable downloader. Safe to re-run: always resumes from the
current file size and never writes concurrently. Verifies gzip on completion.
Usage: python dl_single.py <url> <outfile> [expected_total]
"""
import gzip, os, ssl, subprocess, sys, time, urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE

url = sys.argv[1]
out = sys.argv[2]
total = int(sys.argv[3]) if len(sys.argv) > 3 else None

if total is None:
    r = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60, context=CTX)
    total = int(r.headers["Content-Length"])

have = os.path.getsize(out) if os.path.exists(out) else 0
print("total=%d have=%d" % (total, have), flush=True)

if have >= total:
    print("already complete", flush=True)
else:
    stall = 0
    while have < total:
        before = have
        try:
            req = urllib.request.Request(url, headers={"Range": "bytes=%d-%d" % (have, total - 1)})
            with urllib.request.urlopen(req, timeout=300, context=CTX) as r, open(out, "ab") as f:
                while have < total:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    f.write(c)
                    have += len(c)
        except Exception as ex:
            print("  retry at %d (%s)" % (have, str(ex)[:50]), flush=True)
            time.sleep(5)
        if have == before:
            stall += 1
            if stall > 8:
                print("STALLED at %d" % have, flush=True)
                sys.exit(1)
        else:
            stall = 0
        print("  progress %d/%d (%.1f%%)" % (have, total, 100.0 * have / total), flush=True)

print("downloaded %d" % os.path.getsize(out), flush=True)
rc = subprocess.run(["gzip", "-t", out], capture_output=True)
ok = rc.returncode == 0
print("GZIP_OK=%s" % ok, flush=True)
if not ok:
    print(rc.stderr.decode("utf-8", "ignore")[:150], flush=True)
    sys.exit(1)
with gzip.open(out, "rt", errors="ignore") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        print("HEADER:", line.rstrip()[:250], flush=True)
        break
print("DONE", flush=True)
