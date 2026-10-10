#!/usr/bin/env python3
"""Validate and time exact symbolic TFIM propagation (standard library only)."""
import argparse
import csv
import cmath
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
FIELDS = ['engine','qubits','steps','mode','trial','build_s','eval_s','numeric_s','terms','observable','max_error']

def capture(command, env=None):
    return subprocess.check_output(command, cwd=ROOT, env=env, text=True, timeout=180)

def invoke(engine, env):
    command = [str(ROOT/'target/release/examples/symbolic_benchmark')] if engine=='ppvm' else ['julia','--project=benchmarks/symbolic','-t1','benchmarks/symbolic/surrogate.jl']
    return capture(command, os.environ | env | {'JULIA_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','RAYON_NUM_THREADS':'1'})

def coefficients(text):
    result={}
    for line in text.splitlines():
        if line.startswith('C,'):
            _, sample, word, value=line.split(',')
            result[(int(sample),word)]=float(value)
    if not result:
        raise ValueError('runner emitted no verification coefficients')
    return result

def dense_expectation(n, depth, mode, sample):
    """Independent state-vector check in forward Schrödinger gate order."""
    count=2 if mode=='shared' else depth*(2*n-1)
    values=[0.07+((i*17+sample*11)%53)/100 for i in range(count)]
    state=[0j]*(1<<n);state[0]=1+0j
    gate_index=0
    for _ in range(depth):
        for q in range(n):
            theta=values[0 if mode=='shared' else gate_index];gate_index+=1
            c=math.cos(theta/2);t=-1j*math.sin(theta/2)
            for k in range(1<<n):
                if not k&(1<<q):
                    partner=k^(1<<q);a,b=state[k],state[partner]
                    state[k],state[partner]=c*a+t*b,t*a+c*b
        for q in range(n-1):
            theta=values[1 if mode=='shared' else gate_index];gate_index+=1
            for k in range(1<<n):
                sign=1 if ((k>>q)&1)==((k>>(q+1))&1) else -1
                state[k]*=cmath.exp(-1j*theta*sign/2)
    return sum(abs(v)**2*(1 if k&1==0 else -1) for k,v in enumerate(state))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--qubits',default='3,4'); p.add_argument('--steps',default='1,2,3,4,5,6')
    p.add_argument('--trials',type=int,default=7); p.add_argument('--batch',type=int,default=2000)
    p.add_argument('--out',type=Path,default=ROOT/'benchmarks/symbolic/results')
    args=p.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    raw=[]; summary=[]; checks=[]
    for n in map(int,args.qubits.split(',')):
        for depth in map(int,args.steps.split(',')):
            for mode in ['independent','shared']:
                env={'QUBITS':str(n),'STEPS':str(depth),'MODE':mode,'TRIALS':str(args.trials),'BATCH':str(args.batch)}
                print(f'Checking and timing n={n}, depth={depth}, {mode}',flush=True)
                a,b=[coefficients(invoke(engine,env|{'VERIFY':'1'})) for engine in ['ppvm','paulipropagation-jl']]
                error=max(abs(a.get(k,0)-b.get(k,0)) for k in a.keys()|b.keys())
                if not math.isfinite(error) or error>1e-10:
                    raise ValueError(f'cross-engine coefficient mismatch: {error}')
                dense_error=max(abs(sum(value for (idx,w),value in a.items() if idx==sample and all(c in 'IZ' for c in w))-dense_expectation(n,depth,mode,sample)) for sample in range(5))
                if dense_error>1e-10:
                    raise ValueError(f'dense-state mismatch: {dense_error}')
                checks.append(dict(qubits=n,steps=depth,mode=mode,max_coefficient_error=error,max_dense_error=dense_error))
                for engine in ['ppvm','paulipropagation-jl']:
                    rows=list(csv.DictReader(invoke(engine,env).splitlines(),fieldnames=FIELDS))
                    if len(rows)!=args.trials:
                        raise ValueError('wrong trial count')
                    for row in rows:
                        if row['engine']!=engine or float(row['max_error'])>1e-10:
                            raise ValueError('invalid benchmark result')
                        for name in ['build_s','eval_s','numeric_s']:
                            if not math.isfinite(float(row[name])) or float(row[name])<=0:
                                raise ValueError('invalid timing')
                    raw.extend(rows)
                    item={key:rows[0][key] for key in ['engine','qubits','steps','mode','terms','observable']}
                    item.update({key:statistics.median(float(r[key]) for r in rows) for key in ['build_s','eval_s','numeric_s']})
                    saving=item['numeric_s']-item['eval_s']
                    item['break_even_evaluations']=math.ceil(item['build_s']/saving) if saving>0 else ''
                    summary.append(item)
                # Save after each case, so interrupted runs keep their measured data.
                for name,rows in [('samples.csv',raw),('summary.csv',summary),('validation.csv',checks)]:
                    with (args.out/name).open('w',newline='') as f:
                        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    metadata={'source_sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ['examples/symbolic_benchmark.rs','benchmarks/symbolic/surrogate.jl','benchmarks/symbolic/run.py']},'utc':datetime.now(timezone.utc).isoformat(),'platform':platform.platform(),'machine':platform.machine(),
              'ppvm_revision':capture(['git','rev-parse','HEAD']).strip(),'rust':capture(['rustc','--version']).strip(),
              'julia':capture(['julia','--version']).strip(),'paulipropagation_revision':'b2125c7d2dc61626dbd50f44a47064f290dc5ce9','julia_package':capture(['julia','--project=benchmarks/symbolic','-e','using PauliPropagation; println(pkgversion(PauliPropagation))']).strip(),
              'trials':args.trials,'batch':args.batch,'threads':1,'arguments':vars(args)|{'out':str(args.out)}}
    if sys.platform=='darwin':
        metadata['hardware']=capture(['sysctl','-n','machdep.cpu.brand_string']).strip()
    (args.out/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(f'Validated {len(checks)} cases. Results: {args.out}',flush=True)

if __name__=='__main__':main()
