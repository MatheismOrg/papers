from pathlib import Path
import subprocess
D=Path(__file__).resolve().parent;root=D.parent
header=D/'pdf_header.tex'
header.write_text('\\usepackage{sectsty}\n\\allsectionsfont{\\raggedright}\n\\usepackage{float}\n\\floatplacement{figure}{H}\n')
cmd=['pandoc',str(root/'PAPER/A42_T19.md'),'-o',str(root/'PAPER/A42_T19.pdf'),'--pdf-engine=xelatex',
 '--resource-path='+str(D),'-H',str(header),'-V','mainfont=DejaVu Serif','-V','sansfont=DejaVu Sans','-V','monofont=DejaVu Sans Mono',
 '-V','geometry:margin=0.8in','-V','fontsize=10pt','-V','colorlinks=true']
with (D/'pdf_build.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
