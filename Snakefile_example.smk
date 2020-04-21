import os
import re
import pandas as pd
import math
import numpy as np
from snakemake.shell import shell

#directory with input files, and where output files will be generated
DIR="/g/steinmetz/vijayram"


def genomeList(dir):
    """
    (str)-> (list of str)
    gets list of genomes to sequenced from folder
    """
    listfiles=[]
    for file in os.listdir(dir):
        if file.endswith(".fasta"):
            filename=file.split('.')[0]
            listfiles.append(filename)
    return listfiles
def checkGenome(filename):
    """
    (str)-> (bool)
    checks genome to check if its circular or linear
    """
    file=str(filename) + ".fasta"
    # print (file)
    f=open(file, "r")
    if f.mode=="r":
        seqs=f.read()
    else:
        raise ValueError("File did not open")
        shell("File did not open")
    flag=bool(re.search('chr',seqs))
    return flag


def returnList(dir):
    """
    (str)-> (list of str)
    returns list of genomes after splitting into circular and linear

    """
    listfiles=[]
    for file in os.listdir(dir):
        if file.endswith(".fasta"):
            filename=file.split('.')[0]
            f1=dir +"/"+filename + ".fasta"
            f=open(f1, "r")
            seqs=f.read()
            flag=bool(re.search('chr',seqs))
            if flag:
                listfiles.append(filename+"_lin")
                listfiles.append(filename+"_circ")
            else:
                listfiles.append(filename+"_circ")
    return listfiles 


def genomeSize(file):
    """(str) -> (int)
    returns reference genome size
    """
    os.system('samtools faidx "'+ file +'"')
    faidx_file=file+'.fai'
    df=pd.read_csv(faidx_file, sep='\t', header=None)
    value=df.sum()[1]
    print("The genome size is"+ str(value))
    return value

def read_fasta(fname):
    """ (str) -> (list of tuples)
    reads a fasta file and returns tuples of sequence names and sequence strings.

    """
    f=open(fname)
    lines=f.readlines()
    line=[x.strip('\n') for x in lines]

    #print (line)
    sequences=[]
    for i in range(len(line)):
        if line[i].startswith('>'):
            tup1=line[i]
            tup2='-'
            j=i+1
            while j<len(line) and line[j].startswith('>')== False:
                tup2+=line[j]
                #print (tup2)
                j+=1
            tup=(tup1.strip('>'),tup2.strip('-'))
#             sequences[0].append(tup[0])
#             sequences[1].append(tup[1])
            sequences.append(tup)

    return sequences # a list of (sequence_name, sequence) tuples

def getGFF(string):
    """(str)->(str)
    gets the path of the annotation file corresponding to the wildcard
    """
    splitup=string.split('_')
    filename=splitup[0]
    path="/g/steinmetz/project/IESY/genomes/annotations/scramble/gff/"+filename+".gff"
    # print (path)
    return path

def getTelomeres(output,type):
    """(str) -> (dataframe)
    get table of telomeres in each chromosome with wildcard
    """
    df=pd.read_csv(getGFF(output),sep='\t', header=None)
#     print (df)
    df2=df.loc[df[2]=='telomere']
    df3=df2.iloc[:,[0,2,3,4]]
    # print (df3)
    sequences=read_fasta('/g/steinmetz/vijayram/split/'+output+'.fasta')
    # print(sequences)
    dfseq=pd. DataFrame(sequences,columns=["Name","Sequence"])
    dfseq = dfseq.replace(to_replace='None', value=np.nan).dropna()
    # print (dfseq)
    dftelomeres=pd.DataFrame(index=np.array(dfseq.iloc[:,dfseq.columns.get_loc('Name')]),columns=['Tel1','Tel2'])
    string=''
    for i in range(len(dfseq['Name'])):
        sequence_name=dfseq.iloc[i,dfseq.columns.get_loc('Name')]
        sequence=dfseq.iloc[i,dfseq.columns.get_loc('Sequence')]        
        # print (sequence_name)
        df4=df3.loc[df3[0]==sequence_name]
        # print (sequence[int(df4.iloc[0,2]):int(df4.iloc[0,3])])
        # print(sequence_name)
        # print(len(sequence))
        dftelomeres.loc[sequence_name,'Tel1']=sequence[(int(df4.iloc[0,2])-1):(int(df4.iloc[0,3])-1)]
        string=string+' '+sequence[(int(df4.iloc[0,2])-1):(int(df4.iloc[0,3])-1)]
        if len(df4[0])<2:
            dftelomeres.loc[sequence_name,'Tel2']=math.nan
        else:
            # print(int(df4.iloc[1,3]))
            dftelomeres.loc[sequence_name,'Tel2']=sequence[(int(df4.iloc[1,2])-1):(int(df4.iloc[1,3])-1)]
            string=string+' '+sequence[(int(df4.iloc[1,2])-1):(int(df4.iloc[1,3])-1)]

        # df5=pd.DataFrame({'Tel1':sequence[int(df4.iloc[0,2]):int(df4.iloc[0,3])],'Tel2':sequence[int(df4.iloc[0,2]):int(df4.iloc[0,3])]}, index=[sequence_name])
        # # df5=df5.transpose()
        # print (df5)
        # dftelomeres.append(df5,ignore_index=True)
    if type=='s':
        return (string)
    if type=='df':
        return (dftelomeres)

def getCoverage(file, coverage):
    """
    (str, int)-> (int)
    calculates number of reads for genome based on coverage and genome size, assuming average read length
    """
    size=genomeSize(file)
    mean_read_len=8000
    number_reads=round((coverage*size)/mean_read_len)
    return number_reads

SAMPLES=genomeList(DIR+"/genome")#requisite genome fasta files in /genome directory in DIR
STRAINS=returnList(DIR+"/genome")
print (STRAINS)
COVERAGE=['5', '10', '15', '30', '100', '200']#different coverages to sequence at

# SAMPLES='JS599_test2'
# STRAINS='JS599_test2_circ'
# COVERAGE=10

rule all:
    input:
        expand(DIR+"/output1/split/{sample}_circ.fasta", sample=SAMPLES),
        expand(DIR+"/output1/split/{strain}.fasta", strain=STRAINS),
        # expand([DIR+"/errors/training.maf",DIR+"/errors/model_profile"]),
        # expand([DIR+"/output/{output}_sim.fa", DIR+"/output/{output}_sim.log"], output=OUTPUTS)
        # expand([DIR+"/assembly/{output}/assembly.fasta", DIR+"/assembly/{output}/assembly_graph.gfa", DIR+"/assembly/{output}/assembly_info.txt"],output=OUTPUTS),
        expand(DIR+"/output1/assembly/{strain}_cov{coverage}/repeat.pdf", coverage=COVERAGE, strain=STRAINS),
        expand(DIR+"/output1/quast/{strain}_cov{coverage}/report.pdf", coverage=COVERAGE, strain=STRAINS),
        expand(DIR+"/tapestry/{strain}_cov{coverage}/{strain}_cov{coverage}.tapestry_report.html", coverage=COVERAGE, strain=STRAINS)
    threads:1

#splits genome into circular and linear portion
rule circlinear:
    input:
        DIR+"/genome/{sample}.fasta"
    output:
        circular=DIR+"/output1/split/{sample}_circ.fasta",
        linear= DIR+"/output1/split/{sample}_lin.fasta"
    threads:1
    run:
        flag=checkGenome(DIR+"/genome/"+wildcards.sample) #checks to see if genome is linear or circular
        shell("faidx --invert-match --regex '[ERCC|chr]' {input} > {output.circular}")
        if flag==True:
            print('This file contains linear sequences')
            shell("faidx --invert-match --regex '[ERCC|JS]' {input} > {output.linear}")
        else:
            shell("touch {output.linear}")

# rule nanosim_training:
#     input:
#         ref=DIR+"/training/reference.fasta",
#         reads=DIR+"/training/reads.fasta"
#     output:
#         DIR+"/errors/training.maf",
#         DIR+"/errors/model_profile"
#     params:
#         prefix=DIR+"/errors",
#         threads=16 
#     shell:
#         "nanosim-h-train "
#         "-i {input.reads} "
#         "{input.ref} "
#         "{params.prefix}"

#simulates nanosim reads
rule nanosim_h:
    input:
        DIR+"/output1/split/{strain}.fasta" 
    output:
        reads=DIR+"/output1/nanosim/{strain}_cov{coverage}.fa",
    params:
        training= "yeast",#error profile used for simulation
        # max_length= ,
        min_length=50,
    threads:24
    run:
        coverage=int(wildcards.coverage)# coverage of reads
        # print (wildcards.strain)
        flag=checkGenome(DIR+"/output1/split/"+wildcards.strain)#checks type of reference to simulate reads accordingly
        if flag==True:
            type=""
            print('This genome is linear')
        else:
            type="--circular"
            print('This genome is circular')
        args ={
        "DIR":DIR,
        "type":type,#whether circular or linear, inferred from flag
        "reference":input,#reference genome from which reads are simulated
        "training":params.training,
        "prefix":(output.reads).split('.')[0],
        "num":getCoverage(str(input), coverage),#number of reads approximately calculated
        # "max_len":params.max_length,
        "min_len":params.min_length
        }   
        command=("nanosim-h {type} "
            "-n {num} "
            "-p {training} "
            "-o {prefix} "
            # "--max-len {max_len} "
            "--min-len {min_len} "
            # "-m -i "
            "{reference} "
            )
        command = command.format(**args)
        print(command)
        shell(command)

#assembles reads
rule flye:
    input:
        reads=DIR+"/output1/nanosim/{strain}_cov{coverage}.fa" #nanosim reads
    output:
        contigs=DIR+"/output1/assembly/{strain}_cov{coverage}/assembly.fasta",# assembled contigs
        graph=DIR+"/output1/assembly/{strain}_cov{coverage}/assembly_graph.gfa",# assembled graph
        text=DIR+"/output1/assembly/{strain}_cov{coverage}/assembly_info.txt",# assembly stat
        repeat_graph=DIR+"/output1/assembly/{strain}_cov{coverage}/20-repeat/graph_before_rr.gv"
    params:
        dir=DIR,
        reference=DIR+"/output1/split/{strain}"+".fasta",#reference genome
        threads=16,
        outdir=DIR+"/output1/assembly/{strain}_cov{coverage}"#output directory
    threads:16
    run:
        args={
        "DIR":params.dir,
        "reads":input.reads,
        "size":genomeSize(params.reference),#size of reference genome (usually approx size parameter in flye)
        "outdir":params.outdir,#output directory
        "threads":params.threads,
        }
        command=("/g/steinmetz/vijayram/anaconda3/envs/flye/bin/flye "
            "--nano-raw {reads} "
            "-g {size} "
            "-o {outdir} "
            "-t {threads} "
            "--no-trestle" #skipping trestle step where unbridged repears are resolved based on repeat heterogeneities
            )
        command = command.format(**args)
        print(command)
        shell(command)

#converts repeat graph to pdf
rule repeat_graph:
    input: 
        DIR+"/output1/assembly/{strain}_cov{coverage}/20-repeat/graph_before_rr.gv"
    output:
        DIR+"/output1/assembly/{strain}_cov{coverage}/repeat.pdf"
    shell:
        "dot -Tpdf {input} -o {output}" #getting pdf of repeat graph image

#evaluates assembly - misassemblies, NGA/NA values, GC content graphs
rule quast:
    input:
        contigs=DIR+"/output1/assembly/{strain}_cov{coverage}/assembly.fasta" #contigs generated by flye assembly
    output:
        text=DIR+"/output1/quast/{strain}_cov{coverage}/report.txt",
        tsv=DIR+"/output1/quast/{strain}_cov{coverage}/report.tsv",
        tex=DIR+"/output1/quast/{strain}_cov{coverage}/report.tex",
        icarus=DIR+"/output1/quast/{strain}_cov{coverage}/icarus.html",
        pdf=DIR+"/output1/quast/{strain}_cov{coverage}/report.pdf",
        html=DIR+"/output1/quast/{strain}_cov{coverage}/report.html"
    params:
        reference=DIR+"/output1/split/{strain}"+".fasta", #reference genome for comparison
        outdir=DIR+"/output1/quast/{strain}_cov{coverage}", #out directory for quast
        threads=6,
        annotation=lambda  wildcards: getGFF(wildcards.strain) #genome annotation from separate folder
    conda:DIR+"/envs/quast.yml" #environment containing quast
    threads:6
    shell:
        "quast "
        "-o {params.outdir} "
        "-r {params.reference} "
        "-t {params.threads} "
        "-g {params.annotation} "
        # "{type} "
        "{input.contigs}"

#converts fasta reads to fastq for tapestry
rule fastq:
    input:
        reads=DIR+"/output1/nanosim/{strain}_cov{coverage}.fa" #fasta file of reads
    output:
        gzip=DIR+"/output1/fq_reads/{strain}_cov{coverage}.fq.gz" #gzip file of fastq reads
    run:
        fastq=(output.gzip).strip('.gz')
        shell('/g/steinmetz/vijayram/anaconda3/envs/asimu2/bin/reformat.sh in={input.reads} out='+fastq+' qfake=20 ow=t')
        shell('gzip '+fastq)

#evaluates flye assembly - contig-wise evaluation, GC%, read depth
rule tapestry:
    input:
        assembly=DIR+'/output1/assembly/{strain}_cov{coverage}/assembly.fasta', #fasta file of flye assembled contig-wise
        reads=DIR+"/output1/fq_reads/{strain}_cov{coverage}.fq.gz" #gzipped fastq files
    output:
        report=DIR+"/tapestry/{strain}_cov{coverage}/{strain}_cov{coverage}.tapestry_report.html"
    params:
        dir="{strain}_cov{coverage}" #out directory within tapestry folder
    run:
        flag=checkGenome(DIR+"/output1/split/"+wildcards.strain) #checking if genome is linear or circular
        if flag==True: 
            tel=getTelomeres(wildcards.strain,'s') #telomeres for linear genome
            print('This genome is linear')
        else:
            tel="" #no telomeres for circular genome
            print('This genome is circular')
        # if os.path.isdir(DIR+'/tapestry_output')==False:
        #     os.mkdir(DIR+'/tapestry_output')
        args={
        "DIR":DIR,
        "assembly":input.assembly,
        "reads":input.reads,
        "outdir":params.dir,
        "telomeres":tel
        }
        command=("cd {DIR}/tapestry \n" #running in tapestry folder itself (needs to be run that way)
            "./weave "
            "-a {assembly} "
            "-r {reads} "
            "-o {outdir} "
            # "-t {telomeres}"
            )
        command = command.format(**args)
        # print(command)
        os.system(command)
        cwd=os.getcwd()
        print (cwd)

