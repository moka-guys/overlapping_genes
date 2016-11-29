'''
Created on 23 Nov 2016

@author: aled
'''
import pandas as pd

class overlapping_genes():
    def __init__(self):
        # variables used to read and write bed files 
        self.folder='/home/aled/Documents/161125_coverage_reports/bedfiles/'
        #input bedfile
        self.bed_file=self.folder+"Pan647data.bed"
        # input refseq file
        self.refseq_file=self.folder+"Pan647dataRefSeqFormat.txt"
        
        # empty variable for panel number
        self.panname = ""
        
        # empty dict to hold the bed file  
        self.dict = {}
        
        # dict for start and stop of each gene 
        self.gene_start_stop_dict={}
        
        # genes with a start < stop of another gene
        self.overlap_genes=[]
        
        
    def loop_through_bed(self):
        # take the bedfile name and extract the Pan number 
        split_pan=self.bed_file.split("/")
        last_element=split_pan[-1]
        self.panname=last_element.replace("data.bed","")
        
        # open the bedfile 
        bed=open(self.bed_file,'r')
        #ignore header
        for line in bed:
            if line.startswith("#"):
                pass
            else:
                # split
                splitline=line.split("\t")
                chr=splitline[0]
                start=splitline[1]
                stop=splitline[2]
                Accession=splitline[3]
                GeneSymbol=splitline[-1]
                # if genesymbol in dict add a tuple containing chr start and stop to the list
                if GeneSymbol in self.dict:
                    self.dict[GeneSymbol].append((chr,start,stop))
                else:
                    # if gene symbol not already in the list add tuple to empty list 
                    self.dict[GeneSymbol]=[(chr,start,stop)]
        
    def define_gene_boundaries(self):
        # loop through the list of genes 
        for GeneSymbol in self.dict:
            # set empty vars
            genestart=0
            genechr=0
            genestop=0
            # loop through the tuples for the gene
            for i in self.dict[GeneSymbol]:
                # capture the chromosome if first time, else check if it's the same chromosome
                if genechr == 0:
                    genechr = i[0]
                else:
                    assert genechr==i[0]
                #capture the start of the gene
                if genestart == 0:
                    genestart= i[1]
                # if the start is less than the captured start then replace
                elif genestart > i[1]:
                    genestart = i[1]
                # if current gene stop is larger than the captured stop replace
                if genestop < i[2]:
                    genestop=i[2]
            #minus chr
            genechr=genechr.replace('chr','')
            
            #intergerise the chromosome to ensure sorting is int then text (not alphanumeric)
            if genechr in ('M','X','Y'):
                pass
            else:
                genechr=int(genechr)
            # add to dict where key is genesymbol and value is a tuple of chr, start and stop
            self.gene_start_stop_dict[GeneSymbol]=((genechr,int(genestart),int(genestop)))
        #print self.gene_start_stop_dict
    
    def identify_overlaps(self):
        # put the dictionary into a dataframe 
        df=pd.DataFrame.from_dict(self.gene_start_stop_dict, orient='index')
        #name columns
        df.columns = ['chr','start','stop']
        #sort columns on chr and then start, replacing the df
        df.sort_values(['chr','start'],inplace=True)
        
        # for each chromosome
        for i in set(df['chr']):
            # create a mini dataframe (filter on chromosome) 
            a=df[df['chr'] == i]
            #set stop
            stop=0
            # for each gene on chromosome 
            for index, row in a.iterrows():
                # for first gene on chr set stop == gene stop
                if stop==0:
                    stop=row['stop']
                else:
                    # for subsequent genes. If the start of the row is less than the stop of the prev gene add the index (genename) to overlap genes list
                    if row['start'] <= stop:
                        self.overlap_genes.append(index.rstrip())
                    else:
                        # if no overlap change stop to the last base of this gene
                        stop=row['stop']
        # report overlapping genes
        print len(set(self.overlap_genes))
        print self.overlap_genes
    
    def write_seperate_bedfiles(self):
        # open output bedfiles, A+B
        bedfile1=open(self.folder+self.panname+"A.bed",'w')
        bedfile2=open(self.folder+self.panname+"B.bed",'w')
        
        #open original bedfile and loop through
        bed=open(self.bed_file,'r')
        for line in bed:
            #write headers to both files
            if line.startswith("#"):
                bedfile1.write(line)
                bedfile2.write(line)
            else:
                # capture gene symbol from bedfile
                splitline=line.split("\t")
                GeneSymbol=splitline[-1].rstrip()
                #if is a gene which overlaps another gene write to bedfile B
                if GeneSymbol in self.overlap_genes:
                    bedfile2.write(line.replace('chr',''))
                else:
                    #otherwise write to bedfile A
                    bedfile1.write(line.replace('chr',''))
                
        bedfile1.close()
        bedfile2.close()
        bed.close()
    
    def split_refseqfile(self):
        # repeat for refseq file
        # open refseq file
        refseq=open(self.refseq_file,'r')
        #open output refseq files
        refseqfile1=open(self.folder+self.panname+"A_RefSeqFormat.txt",'w')
        refseqfile2=open(self.folder+self.panname+"B_RefSeqFormat.txt",'w')
        
        for line in refseq:
            # write headers to both output files
            if line.startswith("#"):
                refseqfile1.write(line)
                refseqfile2.write(line)
            else:
                #capture genesymbol
                splitline=line.split("\t")
                genename=splitline[1]
                #if a gene which overlaps write to refseq file B
                if genename in self.overlap_genes:
                    refseqfile2.write(line)
                else:
                    #otherwise write to refseq file A
                    refseqfile1.write(line)
        
        #close files
        refseq.close()
        refseqfile1.close()
        refseqfile2.close()        
    
    def create_sambamba_bed(self):
        '''convert to a bedfile with all columns needed for sambamba'''
        #open input bedfile
        bed=open(self.bed_file,'r')
        # open output sambamba file
        sambamba=open(self.folder+self.panname+"sambamba.bed",'w')
        #loop through bedfile
        for line in bed: 
            strand = "X"
            #write headers
            if line.startswith("#2016"):
                sambamba.write(line)
            elif line.startswith("#Chr"):
                sambamba.write("#chr\tstart\tstop\tname\tscore\tstrand\ttranscript\tgene symbol\n")
            else:
                # capture the required info from bedfile 
                splitline=line.split("\t")
                chr=str(splitline[0].replace("chr",""))
                start=str(splitline[1])
                stop=str(splitline[2])
                gene_acc=splitline[-1].rstrip()
                splitgene_acc=gene_acc.split(";")
                genesymbol=splitgene_acc[0].rstrip()
                transcripts=splitgene_acc[1]
                
                #create name and score values
                F3=chr+"-"+start+"-"+stop 
                F4 = "0"
                #open refseq file to capture strand
                refseq=open(self.refseq_file,'r')
                #loop through refseq file
                for line2 in refseq:
                    #skip if strand is already captured or if a header
                    if strand !="X" or line2.startswith("#"):
                        pass
                    else:
                        #capture the gene symbol
                        splitline2=line2.split("\t")
                        gene_acc2=splitline2[1].rstrip()
                        splitgene_acc2=gene_acc2.split(";")
                        genesymbol2=splitgene_acc2[0].rstrip()
                        transcripts2=splitgene_acc2[1]
                        # if gene symbol in refseq file matches that of bedfile capture the strand
                        if str(genesymbol2)== str(genesymbol):
                            strand=splitline2[3].rstrip()
                # important to close refseq file so it loops through again for next line in bedfile
                refseq.close()
                
                #write sambamba bedfile line
                tab="\t"
                sambamba.write(chr+tab+start+tab+stop+tab+F3+tab+F4+tab+strand+tab+transcripts+tab+genesymbol+"\n")
                
                    
if __name__ == '__main__':
    a=overlapping_genes()
    a.loop_through_bed()
    a.define_gene_boundaries()
    a.identify_overlaps()
    a.write_seperate_bedfiles()
    a.split_refseqfile()
    #a.create_sambamba_bed()