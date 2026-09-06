#!/usr/bin/env python3
"""Archive completed-run solver logs, verify every byte, then reclaim inodes.

Call only after confirming no jobs are writing within this run root.
"""
import argparse, csv, hashlib, io, sys, tarfile
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--run-root',required=True);p.add_argument('--archive',required=True)
args=p.parse_args();root=Path(args.run_root).resolve();archive=Path(args.archive)
if archive.exists():raise FileExistsError(archive)
# Unpublished atomic-save remnants contain no authoritative allocation.
pending=list(root.glob('**/.pending-*'))
for path in pending:
    if path.is_file() and not path.is_symlink():path.unlink()
print('Removed',len(pending),'unpublished temporary files',flush=True)
logs=sorted(root.glob('**/solver-logs/*.log'))
if not logs:
    print('No solver logs: all draws may be infeasible or undefined',flush=True)
    sys.exit(0)
manifest=[]
with tarfile.open(archive,'w:gz',compresslevel=1) as tar:
    for i,path in enumerate(logs):
        if path.is_symlink() or root not in path.resolve().parents:raise RuntimeError('Unexpected path')
        data=path.read_bytes();name=str(path.relative_to(root))
        item=tarfile.TarInfo(name);item.size=len(data);tar.addfile(item,io.BytesIO(data))
        manifest.append(dict(path=name,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data)))
        if (i+1)%10000==0:print('Archived',i+1,'logs',flush=True)
    text=io.StringIO();writer=csv.DictWriter(text,fieldnames=['path','sha256','bytes']);writer.writeheader();writer.writerows(manifest)
    data=text.getvalue().encode();item=tarfile.TarInfo('MANIFEST.csv');item.size=len(data);tar.addfile(item,io.BytesIO(data))
expected={r['path']:r for r in manifest};verified=set()
with tarfile.open(archive,'r:gz') as tar:
    for item in tar:
        if item.name=='MANIFEST.csv':continue
        row=expected[item.name];data=tar.extractfile(item).read()
        if len(data)!=row['bytes'] or hashlib.sha256(data).hexdigest()!=row['sha256']:raise RuntimeError('Archive verification failed')
        verified.add(item.name)
if verified!=set(expected):raise RuntimeError('Archive inventory mismatch')
for path in logs:
    row=expected[str(path.relative_to(root))]
    if hashlib.sha256(path.read_bytes()).hexdigest()!=row['sha256']:raise RuntimeError('Source changed while archiving')
    path.unlink()
print('Verified and archived',len(logs),'solver logs:',archive,flush=True)
