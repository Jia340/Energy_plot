import argparse
import string
import sys,os
from RNAstructure import RNAstructure
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import SeqIO
plt.ioff()
#import Duplex_Energy

def ParseArg():
  p=argparse.ArgumentParser( description = "Plot duplex fold energy distribution between miRNA and target RNA" )
  p.add_argument("--miRNA",type=str,required=True,help="bed6 file of miRNAs.")
  p.add_argument("-T","--target",type=str,required=True,help="file of target RNA.Contains gene name.")
  p.add_argument("--exons",type=str,required=True,help="annotation file downloaded from UCSC table browser. Contains chr, start, end, name, strand, exons. Recommand to be filtered before.")
#  p.add_argument("-o","--output",type=str,help="specify output file")
  p.add_argument("-A","--Ago",type=str,help="bed3 file+name. If specified, positions will be plotted. Recommand to be filtered before.")
  p.add_argument("--smRNA",type=str,help="bed3 file+name. If specified,small RNA positions will be plotted")
  p.add_argument("-H","--HiC", type=str, help="bed3 file+name. Extract from the strong interaction file of RNA Hi-C. If specified, interaction sites will be plotted")
  p.add_argument('-g',"--genomeFa",type=str,default="/home/yu68/Software/bowtie-0.12.7/indexes/mm9.fa",help="genomic sequence,need to be fadix-ed")
  p.add_argument("-R","--RNAstructureExe",type=str,default="/home/yu68/Software/RNAstructure/exe/",help="folder of RNAstrucutre suite excutable")
  p.add_argument('-s','--samtool_path',dest='spath', type=str,metavar='samtool_path',help="path for the samtool program",default='samtools')
  if len(sys.argv)==1:
    print >>sys.stderr, p.print_help()
    sys.exit(0)
  return p.parse_args()


def Plot_energy(Seq1,Seq2,miRNA,output,RNAstructureExe,Ago=[],smRNA=[], HiC=[]):
  energies=[]
  e_p=[]
 
  RNA_prog = RNAstructure(exe_path=RNAstructureExe)
  c=1

  print(len(Seq2)-len(Seq1))
  for i in range(0,len(Seq2)-len(Seq1)):
    if i%100==0: print >>sys.stderr, "%d\r"%(i),
    energy=RNA_prog.DuplexFold(Seq1,Seq2[i:i+len(Seq1)])
    energies.append(energy)
    e_p.append(i+int(len(Seq1)/2))


  energies=np.array(energies)
  plt.figure(figsize=(50,8))
  pl1=plt.plot(e_p,energies,color="#4169E1")
  plt.xlim=(e_p[0],e_p[-1])

  if Ago!=[]:
    Ago_left=[x[0] for x in Ago]
    Ago_width=[x[1]-x[0]+1 for x in Ago]

  pl2=plt.bar(Ago_left,[0.5]*len(Ago_left),Ago_width,label="Ago HITS-CLIP",facecolor="#FA8072",edgecolor="#FA8072")    
  
    
  if smRNA!=[]:
    smRNA_left=[x[0] for x in smRNA]
    smRNA_width=[x[1]-x[0]+1 for x in smRNA]
    pl3=plt.bar(smRNA_left,[0.5]*len(smRNA_left),smRNA_width,bottom=[2]*len(smRNA_left),label="small RNA-Seq",facecolor="#006400",edgecolor="#006400")
#  Output2=open("local_minimum_seq.txt","w")
  
  if HiC!=[]:
    HiC_left=[x[0] for x in HiC]
    HiC_width=[x[1]-x[0]+1 for x in HiC]
    pl4=plt.bar(HiC_left,[0.5]*len(HiC_left),HiC_width,bottom=[4]*len(HiC_left),label="RNA Hi-C",facecolor="#FFD700",edgecolor="#FFD700")
  
  if(miRNA[0]!=-1 and miRNA[1]!=-1):  
    plt.bar(miRNA[0],0.5,miRNA[1]-miRNA[0]+1,facecolor="none",edgecolor="#000000")
    plt.bar(miRNA[0],0.5,miRNA[1]-miRNA[0]+1,bottom=2,facecolor="none",edgecolor="#000000")
 
  plt.ylabel("Energy")
  '''
  Ago_patch = mpatches.Patch(color='#FA8072', label='Ago HITS-Clip')
  smRNA_patch= mpatches.Patch(color='#006400', label='smallRNA-Seq')
  HiC_patch= mpatches.Patch(color='#FFD700',label='RNA Hi-C')
  plt.legend(handles=[Ago_patch,smRNA_patch,HiC_patch])
  '''
  plt.legend(loc="lower right")
  plt.savefig(output+".pdf")
  plt.close()    
  return

def Overlap(c1,s1,e1,c2,s2,e2):
  if c1!=c2:
    return 0
  if c1==c2:
    if (e1<s2 or e2<s1):
      return 0
    return 1

def Contains(c1,s1,e1,c2,s2,e2):
  if c1!=c2:
    return 0
  if c1==c2:
    if(s1<s2 and e1>e2) or (s2<s1 and e2>e1):
      return 1
    return 0 

def Convert_Pos(exon_s,exon_e,s,e,trans_l,strand):
  if strand=="+":
    n_s=s-exon_s+trans_l+1
    n_e=e-exon_s+trans_l+1
        
  if strand=="-":
    n_s=exon_e-e+trans_l+1
    n_e=exon_e-s+trans_l+1
  return([max(trans_l+1,n_s),min(n_e,trans_l+1+exon_e-exon_s)])    

rev_table=string.maketrans('ACGTacgt', 'TGCATGCA')
def revcomp(seq, rev_table):
  return seq.translate(rev_table)[::-1]

def fetchSeq(chro,start,end,strand,fasta,s_path):
  ''' s_path is the path of samtools  '''
  region = '%s:%d-%d'%(chro,start,end)
  seq="".join(os.popen(s_path+" faidx "+fasta+" "+region).read().split('\n')[1:])
  if strand=="-":
    seq = revcomp(seq,rev_table)
  return seq

def Main():
  args=ParseArg()

  ####index all the input files####
  Ago={}
  smRNA={}
  HiC={}
  Exons={}
  miRNA=[]
  target=[]
  Input1=open(args.target,"r")
  for line in Input1.readlines():
    line=line.strip()
    target.append(line)
    Exons[line]=[]
#  print target  

  Input2=open(args.miRNA,"r")
  for line in Input2.readlines():
    line=line.strip().split("\t")
    miRNA.append(line)
#  print miRNA
  if len(miRNA)!=len(target):
    print >>sys.stderr, "Error: miRNA and target doesn't match!" 
    exit()

  exon_file=open(args.exons,"r")
  for line in exon_file.readlines():
    line=line.strip().split("\t")
    ##only store those whose name is in the targets##
    if line[3] in Exons:
      Exons[line[3]].append(line)
#  print Exons

  if args.Ago!="":
    Ago_file=open(args.Ago,"r")
    while True:
      line=Ago_file.readline()
      if line=="":break
      line=line.strip().split("\t")
      if line[3] in target:
        if line[3] not in Ago:
          Ago[line[3]]=[line]
        else:
          Ago[line[3]].append(line)
#  print Ago        
 # print len(Ago["Wtap"])

  if args.smRNA!="":
    smRNA_file=open(args.smRNA,"r")
    while True:
      line=smRNA_file.readline()
      if line=="":break
      line=line.strip().split("\t")
      if line[3] in target:
        if line[3] not in smRNA:
          smRNA[line[3]]=[line]
        else:
          smRNA[line[3]].append(line)
  
  if args.HiC!="":
    HiC_file=open(args.HiC,"r")
    while True:
      line=HiC_file.readline()
      if line=="":break
      line=line.strip().split("\t")
      if line[3] in target:
        if line[3] not in HiC:
          HiC[line[3]]=[line]
        else:
          HiC[line[3]].append(line)
          
#  print len(smRNA["Wtap"])
  
#  print smRNA
  
  for i in range(len(target)):
#    print len(target)
    seq1=fetchSeq(miRNA[i][0],int(miRNA[i][1]),int(miRNA[i][2]),miRNA[i][5],args.genomeFa,args.spath)
    
 #   print seq1
    target_name=target[i]
    print target_name
    for e in Exons[target_name]:
      Ago_temp=[]
      smRNA_temp=[]
      HiC_temp=[]
      miRNA_s=-1
      miRNA_e=-1
      e_s=e[6].split(",")
      e_e=e[7].split(",")
      
#      print e_s
#      print e_e
      seq2=""
      exons_len=[int(x)-int(y) for x,y in zip(e_e[:-1],e_s[:-1])]
#      print exons_len
      if e[5]=="+":
        trans_len=[0]+[sum(exons_len[:k+1]) for k in range(len(exons_len))]
      elif e[5]=="-":
        trans_len=[0]+[sum(exons_len[-k:]) for k in range(1,len(exons_len))]
        trans_len=trans_len[::-1]
#      print trans_len
      
      for j in range(len(e_s)-1):           
        e_s[j]=int(e_s[j])+1
        seq2=seq2+fetchSeq(e[0],int(e_s[j]),int(e_e[j]),"+",args.genomeFa,args.spath)

        if Overlap(miRNA[i][0], int(miRNA[i][1]),int(miRNA[i][2]),e[0],int(e_s[j]),int(e_e[j]))==1:
           if(e[5]=="+"):
             miRNA_s=int(miRNA[i][1])-int(e_s[j])+trans_len[j]
             miRNA_e=int(miRNA[i][1])-int(e_e[j])+trans_len[j]
           elif(e[5]=="-"):
             miRNA_s=int(e_e[j])-int(miRNA[i][2])+trans_len[j]
             miRNA_e=int(e_e[j])-int(miRNA[i][1])+trans_len[j]
                     
        for a in Ago[target_name]:
          if Overlap(a[0],int(a[1]),int(a[2]),e[0], int(e_s[j]), int(e_e[j]))==1:
            Ago_temp.append(Convert_Pos(int(e_s[j]),int(e_e[j]),int(a[1]),int(a[2]),trans_len[j],e[5]))
            
        for s in smRNA[target_name]:
          if Overlap(s[0],int(s[1]),int(s[2]),e[0], int(e_s[j]), int(e_e[j]))==1:
            smRNA_temp.append(Convert_Pos(int(e_s[j]),int(e_e[j]),int(s[1]),int(s[2]),trans_len[j],e[5]))
           
        for h in HiC[target_name]:
          if Overlap(h[0],int(h[1]),int(h[2]),e[0], int(e_s[j]), int(e_e[j]))==1:
            HiC_temp.append(Convert_Pos(int(e_s[j]),int(e_e[j]),int(h[1]),int(h[2]),trans_len[j],e[5]))
                  

#      print Ago_temp
#      print smRNA_temp
#      print miRNA_s,miRNA_e
      if e[5]=="-":
        seq2=revcomp(seq2, rev_table)
#     print seq2      
      output=miRNA[i][3]+"_"+target_name+"-"+e[4]+"_"+str(i)     
#      print len(Ago_temp),len(smRNA_temp),len(HiC_temp) 
      Plot_energy(seq1,seq2,[miRNA_s,miRNA_e],output,args.RNAstructureExe,Ago_temp,smRNA_temp,HiC_temp)      
     
if __name__=='__main__':
  Main()
