"""Multi-threaded range downloader for EBI GWAS summary statistics.
Falls back to single stream if server does not honour Range requests.
Usage: python dl.py <url> <outfile> [n_threads]
"""
import os, ssl, sys, threading, time, urllib.request

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE


def head_size(url):
    req = urllib.request.Request(url, method="HEAD")
    r = urllib.request.urlopen(req, timeout=60, context=CTX)
    n = r.headers.get("Content-Length")
    acc = r.headers.get("Accept-Ranges", "")
    return (int(n) if n else 0), acc


def part(url, start, end, path, idx, results):
    """Download [start,end], resuming from whatever is already on disk."""
    want = end - start + 1
    for attempt in range(6):
        have = os.path.getsize(path) if os.path.exists(path) else 0
        if have >= want:
            results[idx] = have
            return
        try:
            req = urllib.request.Request(
                url, headers={"Range": "bytes=%d-%d" % (start + have, end)})
            with urllib.request.urlopen(req, timeout=300, context=CTX) as r, open(path, "ab") as f:
                while True:
                    c = r.read(1 << 20)
                    if not c:
                        break
                    f.write(c)
        except Exception as ex:
            print("  part %d attempt %d: %s" % (idx, attempt + 1, type(ex).__name__), flush=True)
            time.sleep(3)
    have = os.path.getsize(path) if os.path.exists(path) else 0
    results[idx] = have if have >= want else -1
    if results[idx] < 0:
        print("  part %d FAILED have=%d want=%d" % (idx, have, want), flush=True)


if __name__ == "__main__":
    url = sys.argv[1]
    out = sys.argv[2]
    nthread = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    if os.path.exists(out) and os.path.getsize(out) > 0:
        print("exists, skip: %s (%d bytes)" % (out, os.path.getsize(out)), flush=True)
        sys.exit(0)
    total, acc = head_size(url)
    print("size=%d bytes  accept-ranges=%r" % (total, acc), flush=True)
    t0 = time.time()
    if total == 0 or "bytes" not in acc.lower():
        r = urllib.request.urlopen(url, timeout=600, context=CTX)
        with open(out, "wb") as f:
            while True:
                c = r.read(1 << 20)
                if not c:
                    break
                f.write(c)
        print("single-stream done in %.1fs" % (time.time() - t0), flush=True)
        sys.exit(0)
    step = total // nthread + 1
    parts = []
    for i in range(nthread):
        s = i * step
        e = min(total - 1, s + step - 1)
        if s > e:
            break
        parts.append((s, e, out + ".part%d" % i))
    res = [0] * len(parts)
    ths = []
    for i, (s, e, p) in enumerate(parts):
        t = threading.Thread(target=part, args=(url, s, e, p, i, res))
        t.start()
        ths.append(t)
    for t in ths:
        t.join()
    if any(x < 0 for x in res):
        print("some parts failed, aborting", flush=True)
        sys.exit(1)
    with open(out, "wb") as fo:
        for _, _, p in parts:
            with open(p, "rb") as fi:
                while True:
                    c = fi.read(1 << 20)
                    if not c:
                        break
                    fo.write(c)
            os.remove(p)
    got = os.path.getsize(out)
    if got != total:
        print("SIZE MISMATCH got=%d want=%d" % (got, total), flush=True)
        sys.exit(1)
    dt = time.time() - t0
    print("done %.1f MB in %.1fs -> %.2f MB/s" % (got / 1e6, dt, got / 1e6 / dt), flush=True)
