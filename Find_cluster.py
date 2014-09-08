import sys,argparse
from operator import attrgetter
import copy
import UnionFind
from data_structure import *
import pp
from xplib import DBI

def ParseArg():
    p=argparse.ArgumentParser(description="Find clusters based on the input bed file. The reads should be on the same chromosome")
    p.add_argument("-i","--input",type=str,required=True,help="input file which is a bed file containing all reads")
    p.add_argument('-o','--output',type=str,help="specify output file")
#    p.add_argument("-P","--parallel",type=int,default=5,help="number of workers for parallel computing, default: 5")
#    p.add_argument("-a","--anno",type=str, help="Exon positions. Should be bed3 + strand")
    p.add_argument("-P","--parallel",type=int,default=5,help="number of workers for parallel computing, default: 5")
    if len(sys.argv)==1:
        print >>sys.stderr,p.print_help()
        sys.exit(0)
    return p.parse_args()

def cluster_regions(part):
    N=len(part)
    print N
    cluster_loc=copy.deepcopy(part)
    chro=part[0].chr
    uf_object=UnionFind.UF(N) # union find object
    for i in range(N):
        for j in range(i+1,N):
            if part[i].overlap(part[j]):
                print "overlap"
                uf_object.merge(i,j)
                #update cluster location (expand)
                cluster_loc[uf_object.find(j)].Update(min(part[i].start,part[j].start),max(part[i].end,part[j].end)) #position may become larger after merging
            if part[i]<part[j]:
                break
#            print >>sys.stderr, "j is: %d\r"%(j)
        if i%1000==0: print >> sys.stderr, "  Merging segment for clusters %s,(%d/%d)\r"%(chro,i,N),
    
    c_pool=[] # cluster pool
    for i in range(N):
        c=uf_object.find(i)
        #c_info="%d:%d"%(cluster_loc[c].start,cluster_loc[c].end)
        part[i].Cluster(chro+".%d"%(c))
        c_pool.append(c)
    

    cluster_pool={}
    for c in set(c_pool):
        count=c_pool.count(c)
        cluster_pool[chro+".%d"%(c)]=cluster_loc[c]
        cluster_pool[chro+".%d"%(c)].cluster=count
      
    return (cluster_pool)

'''
def convertPos(exon,exon_length,pos,strand):
    k=0
#    print strand    
    for e in exon:
        if(int(e[1])>=pos and int(e[0])<=pos):            
            if(k==0):
                if(strand=="+"):return(pos-int(e[0])+1)
                else:
      #              print e[0],e[1]
                    return(int(e[1])-pos+1)
            else:
                if(strand=="+"):return(sum(exon_length[0:k])+pos-int(e[0])+1)
                else:
       #             print e[0],e[1],sum(exon_length[0:k])
                    return(sum(exon_length[0:k])+int(e[1])-pos+1)
        else:
            if(strand=="+"):
                if(int(e[0])>pos):
                    return(sum(exon_length[0:k]))
                    break     
            else:
                if(int(e[1])<pos):
                    return(sum(exon_length[0:k]))
                    break
        k+=1  
'''
     
def Main():
    args=ParseArg()
    inp=open(args.input,'r')
#    inp2=open(args.anno,'r')
    output=open(args.output,'w')
    k=0
    part=[]
    chr_list=[]
    ncpus=args.parallel
   # exons=[]
   # exon_length=[]

    for line in inp.read().split('\n'):
        if line=='':continue
        line=line.strip().split('\t')
        p=annotated_bed('\t'.join(line[0:3]+[line[5]]),id=k)
        part.append(p)
        k+=1
        if p.chr not in chr_list: chr_list.append(p.chr)
    print >> sys.stderr, "Get total %d pairs"%(k) 
    part=sorted(part,key=attrgetter("start"))
    part=sorted(part,key=attrgetter("chr"))
    
   # chro=part[0].chr
    
  #  for line in inp2.readlines():
   #     line=line.strip().split('\t')
    #    if line[0]!=chro:
     #       print >>stderr, "Wrong chromosome!"
      #      exit()
      #  exons.append(line[1:3])
   
   # Strand=line[3]

   # if(Strand=="+"):
    #    exons=sorted(exons)
   # elif(Strand=="-"):
    #    exons=sorted(exons,reverse=True)
   # else:
    #    print >> stderr, "Wrong strand!"
    #    exit()
   # for e in exons:
    #    exon_length.append(int(e[1])-int(e[0])+1)
   
    ppservers=()
    job_server=pp.Server(ncpus,ppservers=ppservers)
    jobs=[]
#    print ncpus
    for chro in chr_list:
        part_temp=filter(lambda p: p.chr==chro, part)
        if len(part_temp)>0:
            jobs.append(job_server.submit(cluster_regions,(part_temp,),(annotated_bed,),("UnionFind","copy",)))
#    clusters=cluster_regions(part) 
    
    cluster_pool={}
    for job in jobs:
      cluster_pool.update(job())
    
    for Item in cluster_pool.values():
      print >>output, '\t'.join([Item.chr,str(Item.start),str(Item.end)])

    

#    for Item in clusters.values():
     #   print Item.start, Item.end
 #       s1=convertPos(exons,exon_length,int(Item.start),Strand)
  #      e1=convertPos(exons,exon_length,int(Item.end),Strand)
     #   print s1,e1
  #      if(Strand=="+"):
   #         print >>output, '\t'.join([str(s1),str(e1)])   
   #     elif(Strand=="-"):
    #        print >>output, '\t'.join([str(e1),str(s1)])
            
    
if __name__=="__main__":
  Main()
    
  
        
