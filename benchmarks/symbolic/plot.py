#!/usr/bin/env python3
"""Plot the symbolic benchmark in the paper's QuEra palette."""
import argparse
import csv
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--csv',type=Path,default=Path(__file__).parent/'results/summary.csv')
p.add_argument('--qubits',type=int,default=4)
p.add_argument('--out',type=Path,default=Path(__file__).parent/'results/symbolic_propagation.pdf')
args=p.parse_args()
rows=list(csv.DictReader(args.csv.open()))
plt.rcParams.update({'font.family':'sans-serif','font.size':8,'axes.linewidth':0.5,'lines.linewidth':1.2,'lines.markersize':3,'xtick.direction':'in','ytick.direction':'in','pdf.fonttype':42})
colors={'ppvm':'#8A40F7','paulipropagation-jl':'#FF7C24'}
labels={'ppvm':'PPVM','paulipropagation-jl':'PauliPropagation.jl'}
fig,axes=plt.subplots(2,2,figsize=(6.8,3.8),sharex=True,sharey="row")
for col,mode in enumerate(['independent','shared']):
    for engine in colors:
        series=sorted([r for r in rows if int(r['qubits'])==args.qubits and r['mode']==mode and r['engine']==engine],key=lambda r:int(r['steps']))
        x=[int(r['steps']) for r in series]
        axes[0,col].plot(x,[float(r['build_s'])*1e3 for r in series],'-o',color=colors[engine],label=labels[engine])
        axes[1,col].plot(x,[float(r['eval_s'])*1e6 for r in series],'-o',color=colors[engine])
        axes[1,col].plot(x,[float(r['numeric_s'])*1e6 for r in series],'--',color=colors[engine])
    axes[0,col].set_title('Independent angles' if mode=='independent' else 'Shared angles')
    axes[1,col].set_xlabel('Trotter steps')
    for row in range(2):
        ax=axes[row,col];ax.set_yscale('log');ax.set_xticks(x)
        ax.spines[['top','right']].set_visible(False)
        ax.text(-0.14,1.02,f'({"abcd"[row*2+col]})',transform=ax.transAxes,fontweight='bold')
axes[0,0].set_ylabel('Construction time (ms)')
axes[1,0].set_ylabel('Time per parameter set (µs)')
axes[0,0].legend(frameon=False,fontsize=7)
fig.legend(handles=[Line2D([0],[0],color='#2B313E',label='Symbolic evaluation'),Line2D([0],[0],color='#2B313E',linestyle='--',label='Numeric propagation')],loc='lower center',ncol=2,frameon=False,fontsize=7)
fig.tight_layout(pad=1.0,rect=(0,0.055,1,1))
args.out.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(args.out,bbox_inches='tight')
fig.savefig(args.out.with_suffix('.png'),dpi=180,bbox_inches='tight')
