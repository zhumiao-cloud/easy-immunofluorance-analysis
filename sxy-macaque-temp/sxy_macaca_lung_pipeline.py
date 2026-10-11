#!/usr/bin/env python3
from pathlib import Path
import gzip, hashlib, os, shutil, subprocess, time
from concurrent.futures import ThreadPoolExecutor, as_completed

H=Path.home()
SRC=H/"PROJECTS/zhumiao/danxibao2/刘光慧食蟹猴肺"
OUT=H/"PROJECTS/zhumiao/danxibao2/数据清理整理/刘光慧食蟹猴肺"
WORK=H/"PROJECTS/zhumiao/.sxy_macaca_lung_work"
REF=H/"reference/Macaca_fascicularis_5_0_Ensembl100"
CR=H/"cellranger-10.0.0/bin/cellranger"
VAL=WORK/"sxy_macaca_lung_validate.R"
SAMPLES=[("Lung_YM1",924),("Lung_YM2",928),("Lung_YM3",932),("Lung_YM4",936),
("Lung_YF1",940),("Lung_YF2",944),("Lung_YF3",948),("Lung_YF4",952),
("Lung_OM2",956),("Lung_OM3",960),("Lung_OM4",964),
("Lung_OF1",968),("Lung_OF2",972),("Lung_OF3",976),("Lung_OF4",980)]
MEX=("matrix.mtx.gz","features.tsv.gz","barcodes.tsv.gz")

def log(*x): print(time.strftime("%F %T"),*x,flush=True)
def need(x,msg):
    if not x: raise RuntimeError(msg)
def cmd(a,**kw):
    log("RUN",a[0],a[1] if len(a)>1 else "")
    subprocess.run([str(x) for x in a],check=True,**kw)

def fh(p):
    with gzip.open(p,"rb") as f: r=[f.readline().rstrip(b"\r\n") for _ in range(4)]
    need(r[0].startswith(b"@") and r[2].startswith(b"+") and len(r[1])==len(r[3]) and len(r[1])>0,"bad FASTQ "+str(p))
    return r[0].split()[0].removesuffix(b"/1").removesuffix(b"/2"),len(r[1])

def runs(sample,start):
    for lane,n in enumerate(range(start,start+4),1):
        r=f"CRR{130000+n}"; d=SRC/r
        yield lane,r,d/(r+"_f1.fq.gz"),d/(r+"_r2.fq.gz")

def mounted_copy(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True)
    expected=src.stat().st_size
    need(expected>100,"mounted FASTQ has invalid size "+str(src))
    if dst.is_file() and dst.stat().st_size==expected:
        try:
            fh(dst)
            return
        except Exception:
            dst.unlink()
    tmp=dst.with_name(dst.name+".part")
    if tmp.exists() and tmp.stat().st_size>expected:
        tmp.unlink()
    for attempt in range(1,9):
        have=tmp.stat().st_size if tmp.exists() else 0
        if have==expected:
            break
        log("MOUNT_COPY",src.name,"offset",have,"of",expected,"attempt",attempt)
        p=subprocess.run([
            "timeout","900","dd",f"if={src}",f"of={tmp}","bs=8M",
            "iflag=fullblock,skip_bytes",f"skip={have}",
            "oflag=seek_bytes",f"seek={have}","conv=notrunc","status=none"
        ])
        if p.returncode==0 and tmp.exists() and tmp.stat().st_size==expected:
            break
        log("MOUNT_COPY_RETRY",src.name,"rc",p.returncode,
            "bytes",tmp.stat().st_size if tmp.exists() else 0)
        time.sleep(min(10*attempt,60))
    need(tmp.is_file() and tmp.stat().st_size==expected,
         "mounted copy incomplete "+str(src))
    os.replace(tmp,dst)
    fh(dst)

def stage(sample,start):
    d=WORK/"local_fastqs"/sample; d.mkdir(parents=True,exist_ok=True)
    lines=["sample_id\trun_id\tread1\tread2"]
    for lane,r,r1,r2 in runs(sample,start):
        q1=d/f"{sample}_S1_L{lane:03d}_R1_001.fastq.gz"
        q2=d/f"{sample}_S1_L{lane:03d}_R2_001.fastq.gz"
        mounted_copy(r1,q1)
        mounted_copy(r2,q2)
        i1,l1=fh(q1); i2,l2=fh(q2); need(i1==i2,"R1/R2 mismatch "+r)
        need(l1>=16 and l2>=25,"unexpected read length "+r)
        log("SOURCE_OK",sample,r,l1,l2)
        lines.append(f"{sample}\t{r}\t{r1}\t{r2}")
    return d,lines

def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

def put(src,dst):
    dst.parent.mkdir(parents=True,exist_ok=True); hs=sha(src)
    if dst.exists():
        need(sha(dst)==hs,"existing target differs "+str(dst)); return
    t=dst.with_name(dst.name+".uploading")
    need(not t.exists(),"partial target exists "+str(t))
    with open(src,"rb") as a,open(t,"xb") as b: shutil.copyfileobj(a,b,8*1024*1024)
    need(sha(t)==hs,"upload hash mismatch "+str(dst))
    t.rename(dst); need(sha(dst)==hs,"target readback hash mismatch "+str(dst))

def waitref():
    for i in range(720):
        if (REF/"reference.json").is_file() and (REF/"fasta/genome.fa").is_file(): return
        if i%10==0: log("WAIT_REFERENCE",i)
        time.sleep(30)
    raise RuntimeError("reference not ready")

def one(sample,start):
    final=OUT/sample
    if (final/"VALIDATION_PASS.txt").is_file():
        need(all((final/"10x"/x).is_file() for x in MEX),"marked sample incomplete")
        log("SKIP_COMPLETE",sample); return
    fq,manifest=stage(sample,start)
    counts=WORK/"counts"; counts.mkdir(parents=True,exist_ok=True); out=counts/sample/"outs"
    if not (out/"filtered_feature_bc_matrix/matrix.mtx.gz").is_file():
        if (counts/sample).exists(): log("RESUME_EXISTING_CELLRANGER_CHECKPOINT",sample)
        cmd([CR,"count","--id",sample,"--transcriptome",REF,"--fastqs",fq,"--sample",sample,
             "--chemistry","auto","--include-introns","true","--create-bam=false","--localcores","24","--localmem","240"],cwd=counts)
    for a,b in (("filtered_feature_bc_matrix","10x"),("raw_feature_bc_matrix","raw_10x")):
        for n in MEX: put(out/a/n,final/b/n)
    v=WORK/"validated"/sample; v.mkdir(parents=True,exist_ok=True)
    cmd(["Rscript","--vanilla",VAL,final/"10x",sample,v])
    put(v/"seurat_unfiltered_after_cellranger.rds",final/(sample+"_seurat.rds"))
    mf=v/"input_runs.tsv"; mf.write_text("\n".join(manifest)+"\n"); put(mf,final/"input_runs.tsv")
    mk=v/"VALIDATION_PASS.txt"; mk.write_text("MOUNTED_ONLY_SOURCE=TRUE\nREAD10X=PASS\nCREATESEURATOBJECT=PASS\n"); put(mk,final/"VALIDATION_PASS.txt")
    shutil.rmtree(fq,ignore_errors=True)
    log("SAMPLE_COMPLETE",sample)

def main():
    need(SRC.is_dir() and OUT.parent.is_dir(),"danxibao2 mount unavailable")
    need(CR.is_file() and VAL.is_file(),"required executable/script missing")
    waitref()
    # 以3路并发处理不同样本，每路24核、240GB，合计72核、720GB。
    # 只有全部样本通过完整10x和R校验才能报告全部完成。
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures={executor.submit(one,s,n):s for s,n in SAMPLES}
        failures=[]
        for future in as_completed(futures):
            name=futures[future]
            try:
                future.result()
            except Exception as exc:
                log("SAMPLE_FAILED",name,type(exc).__name__,str(exc))
                failures.append(name)
    need(not failures,"failed samples: "+",".join(failures))
    log("ALL_15_SAMPLES_COMPLETE")

if __name__=="__main__": main()
