from RNAstructure import RNAstructure
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import SeqIO
import sys

plt.ioff()

def ParseArg():
  p=argparse.ArgumentParser( description = "Plot duplex fold energy between miRNA and target RNA" )
  p.add_argument("--miRNA",type=str,required=True,help="fasta file of miRNA")
  p.add_argument("--RNA",type=str,required=True,help="fasta file of target RNA")
  p.add_argument("-o","--output",type=str,help="specify output file")
  p.add_argument("-A","--Ago",type=str,help="if specified, positions will be plotted")
  p.add_argument("-R","--RNAstructureExe",type=str,default="/home/yu68/Software/RNAstructure/exe/",help="folder of RNAstrucutre suite excutable")
  if len(sys.argv)==1:
    print >>sys.stderr, p.print_help()
    sys.exit(0)
  return p.parse_args()
  
def Main():
  energies=[]
  args=ParseArg()
  seq1_handle=open(args.miRNA,"rU")
  seq1=SeqIO.parse(seq1_handle,"fasta").next()
  #print seq1.seq
  Seq1=str(seq1.seq)
#  print type(Seq1)
  
  seq2_handle=open(args.RNA,"rU")
  seq2=SeqIO.parse(seq2_handle,"fasta").next()
  Seq2=str(seq2.seq)
  #print seq2.seq
  
  Output=open("energy.txt","w")
  RNA_prog = RNAstructure(exe_path=args.RNAstructureExe)
  print(len(Seq2)-len(Seq1))
  for i in range(len(Seq2)-len(Seq1)):
    energy=RNA_prog.DuplexFold(Seq1,Seq2[i:i+len(Seq1)]) 
    energies.append(energy)
    if i==701:
        print i, Seq2[i:i+len(Seq1)]
    print >>Output, '\t'.join([str(i),str(energy)])

  if(args.Ago!=""):
    ago=open(args.Ago,'r')
    left=[]
    width=[]
    while True:
      line=ago.readline()
      if line=="":break
      line=line.strip().split('\t')
      left.append(int(line[0]))
      width.append(int(line[1])-int(line[0]))
  Output2=open("local_minimum_seq.txt","w")
#  energies=np.array(energies)

  for j in range(len(energies)):  
    if float(energies[j])<-20:
      if (j==len(energies)-1):
        if(j<energies[j]):
          print >>Output2, '\t'.join([str(j),str(energies[j]),Seq2[j:j+len(Seq1)]])
      if (float(energies[j])<float(energies[j-1]) and float(energies[j])<float(energies[j+1])):
        print >>Output2, '\t'.join([str(j),str(energies[j]),Seq2[j:j+len(Seq1)]])

  energies=np.array(energies)
  plt.figure(figsize=(6,4))
  plt.plot(energies,label="Wtap")
  plt.xlim(1,len(Seq2))
 # plt.plot((148,148),(-30,0), color='r')
  plt.bar(left,[0.1]*len(left),width,facecolor='r',edgecolor="r")  
#  plt.xlabel("Transcript position")
  plt.ylabel("Energy")
  plt.savefig(args.output)
  plt.close() 

if __name__=="__main__":
  Main()
     
