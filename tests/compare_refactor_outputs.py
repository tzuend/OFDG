"""Compare full hexadecimal coefficient dumps from two builds."""
import argparse,json,math
from pathlib import Path

def compare(before,after):
    summary={'coefficients':0,'changed_hex_values':0,'max_absolute_difference':0.,'max_scaled_difference':0.,'files':[]}
    files=sorted(before.glob('*-rank*.txt'))
    if not files:raise RuntimeError('No baseline outputs found')
    for old in files:
        new=after/old.name;previous_case='';a=old.read_text().splitlines();b=new.read_text().splitlines()
        if len(a)!=len(b):raise RuntimeError(f'Output shape changed: {old.name}')
        for line,(x,y) in enumerate(zip(a,b),1):
            if x.startswith('case '):
                if x!=y:raise RuntimeError(f'Case changed: {old.name}:{line}')
                previous_case=x;continue
            xv=x.split();yv=y.split()
            if xv[:2]!=yv[:2] or len(xv)!=len(yv):raise RuntimeError('Vector shape changed')
            for i,(vx,vy) in enumerate(zip(xv[2:],yv[2:])):
                p,q=float.fromhex(vx),float.fromhex(vy)
                if not math.isfinite(p) or not math.isfinite(q):raise RuntimeError('Nonfinite output')
                absolute=abs(p-q);scaled=absolute/max(1,abs(p),abs(q))
                summary['coefficients']+=1;summary['changed_hex_values']+=(vx!=vy)
                if absolute>summary['max_absolute_difference']:
                    summary['max_absolute_difference']=absolute
                    summary['worst_absolute_case']={'file':old.name,'case':previous_case,'vector':xv[0],'coefficient':i,'before':p,'after':q}
                summary['max_scaled_difference']=max(summary['max_scaled_difference'],scaled)
        summary['files'].append({'name':old.name,'byte_identical':old.read_bytes()==new.read_bytes()})
    summary['passed']=summary['max_scaled_difference']<=1e-12
    return summary
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('before',type=Path);parser.add_argument('after',type=Path);args=parser.parse_args()
    s=compare(args.before,args.after);print(json.dumps(s,indent=2));raise SystemExit(0 if s['passed'] else 1)
