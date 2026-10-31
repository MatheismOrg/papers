"""Optional source refresh for audit. Downloads and parses literal numeric data only.
Normal run.py reproduction needs only the already supplied parameter CSV.
"""
from pathlib import Path
from urllib.request import urlopen
import ast,csv,json,hashlib
R=Path(__file__).resolve().parent
URL='https://zenodo.org/records/15299349/files/ClusterInfo.py?download=1'
def main():
    data=urlopen(URL,timeout=60).read()
    for n in ast.parse(data.decode()).body:
        if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='cluster_data' for t in n.targets) and isinstance(n.value,ast.List):
            raw=ast.literal_eval(n.value);break
    rows=[]
    for v in raw:
        c=list(map(float,v[2:10]))
        if v[0]=='MACS1720':c[5]=0.
        rows.append(dict(zip(['system','ne0','r0','alpha','re0','beta0','ne1','re1','beta1'],[v[0]]+c)))
    with (R/'clash_gas_parameters.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (R/'source_provenance.json').write_text(json.dumps({'URL':URL,'sha256':hashlib.sha256(data).hexdigest(),'authors':'B. Famaey, L. Pizzuti, I. D. Saltas; gas profiles credited to V. Salzano','paper':'https://doi.org/10.1103/dccw-srks','dataset_DOI':'10.5281/zenodo.15299349','extraction':'AST literal parameter rows, with MACS1720 single-beta treatment explicitly required by source; no third-party code executed','parameter_csv_sha256':hashlib.sha256((R/'clash_gas_parameters.csv').read_bytes()).hexdigest()},indent=2))
if __name__=='__main__':main()
