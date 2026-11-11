#!/usr/bin/env python3
"""Exact finite code audit for A38. Python 3 standard library only.
Task: decode an n-bit unsigned integer pair. Precision n is shared.
Tags 00/01/10 select equal, one-sided and generic; 11 is unused.
These conditional code lengths are not Python file sizes or physical laws.
"""
from fractions import Fraction
from collections import defaultdict
from pathlib import Path
import hashlib,json

def encode(tag,a,b,n):
    if tag not in ('00','01','10') or n<1 or not (0<=a<2**n and 0<=b<2**n):raise ValueError('Invalid input')
    if tag=='00' and a!=b:raise ValueError('Equal format requires a=b')
    if tag=='01' and b!=0:raise ValueError('One-sided format requires b=0')
    return tag+format(a,f'0{n}b')+(format(b,f'0{n}b') if tag=='10' else '')

def decode(bits,n):
    if n<1 or any(b not in '01' for b in bits):raise ValueError('Invalid code')
    tag=bits[:2]
    if tag not in ('00','01','10') or len(bits)!=2+n*(2 if tag=='10' else 1):raise ValueError('Invalid length/tag')
    a=int(bits[2:2+n],2)
    return (a,a) if tag=='00' else (a,0) if tag=='01' else (a,int(bits[2+n:],2))

def audit(n):
    codes=[];weights=defaultdict(Fraction);classes=defaultdict(Fraction)
    for tag in ('00','01','10'):
        pairs=((a,a) for a in range(2**n)) if tag=='00' else ((a,0) for a in range(2**n)) if tag=='01' else ((a,b) for a in range(2**n) for b in range(2**n))
        for a,b in pairs:
            c=encode(tag,a,b,n);assert decode(c,n)==(a,b);codes.append(c);w=Fraction(1,2**len(c));classes[tag]+=w;weights[a,b]+=w
    ordered=sorted(codes)
    assert len(set(codes))==len(codes)
    assert all(not b.startswith(a) for a,b in zip(ordered,ordered[1:]))
    assert all(w==Fraction(1,4) for w in classes.values())
    assert sum(weights.values())==Fraction(3,4)
    assert weights[0,0]==Fraction(2,2**(n+2))+Fraction(1,2**(2*n+2))
    return dict(n=n,codewords=len(codes),distinct_outputs=len(weights),format_weights={k:str(v) for k,v in classes.items()},output_weight_sum=str(sum(weights.values())),zero_pair_weight=str(weights[0,0]))

def main():
    n=16;a=137;short=encode('00',a,a,n);long=encode('10',a,a,n)
    assert decode(short,n)==decode(long,n)==(137,137)
    raw=Path(__file__).read_bytes()
    print(json.dumps(dict(scope='Finite conditional coding illustration, not a physical theory comparison or Kolmogorov minimum.',shared_inputs=['precision n','the decoder and tag definitions'],exhaustive_checks=[audit(n) for n in range(1,6)],matched_example=dict(output=[137,137],n=16,equal_code=short,generic_code=long,equal_bits=len(short),generic_bits=len(long),individual_weight_ratio=2**(len(long)-len(short))),implementation=dict(utf8_bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),charged='Shared implementation; not included in the conditional tag/payload lengths.')),indent=2))
if __name__=='__main__':main()
