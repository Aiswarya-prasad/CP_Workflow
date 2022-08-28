UPDATE: The Bracken, Kraken approach is not suitable for nanopore data. This repository may be updated at a later point in time with more appropriate tools.

This is an experimental snakemake pipeline that can be used to carry out pre-processing and taxonomic classification of nanopore sequencing data (from fast5 to figures) without assembly and using limited resources.

> This is not (yet) intended to be a ready-to-use pipeline. I learnt Snakemake from scratch through the course of this project and would have found a repo like this immensely useful. That is the main motivation I had to create this repository. Most of the code can be easily understood and reused.<br/>

The pipeline combines a set of existing tools for which there are multiple alternatives available in and outside of GitHub. These particular tools were chosen such that they were easy to set up and use and could run on limited resources. This pipeline was written and tested in snakemake version 5.5.4.

## Rulegraph (DAG) for this pipeline

<p align="center">
  <img src="rulegraph.png" width="500" />
</p>


# Tools used

Below is a list of tools used by this pipeline with a link to their page/repo where information about installation and usage can be found.

- [Guppy](https://denbi-nanopore-training-course.readthedocs.io/en/latest/basecalling/basecalling.html) (you need to be part of the nanopore community for access to official docs and [download](https://community.nanoporetech.com/downloads))
- [NanoPlot](https://github.com/wdecoster/NanoPlot)
- [qcat](https://github.com/nanoporetech/qcat)
- [NanoFilt](https://github.com/wdecoster/nanofilt/)
- [Kraken 2](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual)
- [Centrifuge](https://ccb.jhu.edu/software/centrifuge/)

If you will find this repo useful if improved (in terms of documentation, error handling, features etc.), please let me know. Also, if you are a veteran and somehow come across this repo I will be thrilled to hear your sugesstions and more so if you have the time to help me polish it. Finally, if you find any part of this repo useful for your work, do let me know. I will be very happy to hear!

# Setup and use</br>

Firstly, all these tools need to be installed and working. Some of this can be installed through conda. Snakemake can handle this and I will soon add the conda environment set up for this. For now, they need to be installed and set up manually.

In the config.yaml file, make sure to configure the following:</br>

- **RAWDIR**: Path to directory containing links to fast5 directory (or) path to fast5 directory to be used by Guppy</br>
- **ROOT**: Path to project directory to be used by shell commands and scripts if necessary</br>
- **runnames**: List of run names eg. \['Run0', 'Run1', 'Run2', 'Run3', 'Run4'\]</br>
- **barcode_kit**: 'NBD103/NBD104' parameter used by qcat</br>
samples: List of sample IDs eg. \['01', '03', '06', '07', '10', '11', '12', '13', '14', '15', '17', '18', '19', '20', '21', '22'\]</br>
- **sample_dict**: a nested dictionary mapping runs to barcodes and sample IDs. ie. the first level of keys are run names, within each dict which is a value for the key run name, is a dictonary where keys are barcode numbers and values are sample ID associated with each barcode in that particular run. *This is useful for demultiplexing each run and separating out the files based on sample ID rather than barcode as barcodes may not be unique unless coupled with run names leading to long messy names. leave this out and edit Snakefile accordingly if this is not desired*
eg. (written in yaml format will be read by snakemake as a python style nested dictionary)</br>
```yaml
sample_dict:
    'Run0':
        '04': '01'
    'Run1':
        '02': '03'
        '03': '13'
    'Run2':
        '04': '06'
        '05': '07'
        '06': '10'
        '07': '14'
    'Run3':
        '08': '11'
        '09': '12'
        '10': '15'
        '11': '17'
        '12': '18'
    'Run4':
        '01': '19'
        '02': '20'
        '03': '21'
        '04': '22'
```
- **kraken_db**: '$MINIKRAKEN_DB' environment varaible if it is configured in ~/.bashrc. If not, write the path to the kraken database here.</br>
- **centrifuge_db**: path to the centrifuge database</br>


Make sure that these programs are installed and configured properly (eg. added to the path variable in ~/.bashrc as is appropriate).

--------------------------------------------------------------------------------------------------------------------------


## I/O file format <br/>
At the moment, most of the information about output and input formats can be found in the code or as **comments in the code**. Rules handling runs read and write fastq files which can be zipped by including zipped files in the target rule. All the rules handling samples read and write fastq files as zipped files (.fastq.gz). They can be unzipped by including fastq files in the target rule and run fastq files can be zipped and stored elsewhere. rule zip and unzip do take care of this and can be edited as needed (see comments in the code).

## Note about NanoFilt<br/>
input to rule running kraken2 and centrifuge can be edited depending on which fastq files (filtered Q > 7, Q > 10 or unfiltered) are to be used for downstream analysis. Edit this rule (filterSamples) as needed. Refere to the NanoFilt repository for more details.

## directory structure of output after classification<br/>
classified<br/>
└── <sample#><br/>
    ├── bracken<br/>
    ├   ├── genus_report<br/>
    ├   └── species_report<br/>
    ├── centrifuge<br/>
    ├   ├── report<br/>
    ├   └── result<br/>
    └── kraken2_Minidb<br/>
        ├── report<br/>
        └── result<br/>


# There were duplicate reads so flye failed,
# from 03, I removed
# @7d5a02c6-8bf0-4ed1-8b7e-8990b3eb81d4 runid=cc18a835305a98d8cc7f2ff70b00bef64778926f read=13 ch=203 start_time=2019-10-25T15:59:43Z flow_cell_id=FAL22147 protocol_group_id=Exp1_25Oct sample_id=sample3_13 barcode=unclassified
# CTCCACGCTCTTCCCATTCCGAACAGAGAAGTTAGCTCTACCACGCCGATGGTACTACAATGCAATGCGAGAGATAGGTAGCCCTCTTTCAATTTCAAGCCTCGATTACAATTTGTAGTCAGGAGGCTTTTTATTTTTATCTGCTCAGATGTTTTTATCAGCTTGAATGACGTTGATTACGGATTTTTATTATCTTTGCAGCATGATTTGAATGGAAGTAGAATAAGGCAATTGAAGATTTTTTGGTGATAGCAGCAATTGTTTATTGCTGTCGCCTCTCTTGTTGTGTCTCATTTTCTTGTTCGTGATCTAGCAGAGAAACCTGGCCGTCTATGGCGGTTTGGGCTGAGCTATGCGTGTACGCTCAATAATGCCGATGAGAATACTGACGACTTCTGGTAACGAAGGTTATCAGTAATAATTCCATTCCGGTTATCGTAATGGATTCTGAAAACCAATATCAGACTTTTCGTAATGTTGATGGGAGGAAAGATTATGCAGATCCCTTGCTTGTAGCTACCGTCGGTCAGCGGCTCTTGAAACAAGGCAAGAATATTAAAATTGAACTAGATGATTCAACCGCTGACTATATTCGGGTTTGCTATGATGAATCATGATACGACGTCTGTCAGCTTATCCTTATATTCAGTTGGGGGTAGTCATGTACGTGGTTGTGCCATCTTTGCTTTTAGCACTTCCAAAGCCGGAACAGAATAAGGTGTGGGTAAGGATTTGTCTAAAGCTTGCTCATCGGTTGGTGCTCCAATTTCTTCCAGTCTGATGGCTTGGATTGAAATTCTGAAGATGAATTATCGGATGACGAACTCATTCCTGAAATGAACAAGGATATCCAGCGTTTTGTACTGATGGCAGACCGTTTCTCAAAAATAGGTGCCCTGCCAGAACCTGTTCCACCTGAAGGATTTGAACCGAAGTGATGGAATCATATTATAGACCACTGGATCGTCGAACATCAAGAATGTTAAAATGGTGAGAGGAAATACACGAACAGGATATTTATAGATTTAAGATGAATGCACTCTCTCTTTTGAATGGGTGATAGAAAATCTATTCAGAAGAATACTTGTTTAAATTATTTATAGAGGCAAATGGTGGACTGAATTGCGCTCTATAATGGAAGAAATACGAGTATTAAAGCTATTATAGAAGTTTCAGATACGGGAAATGGTGTAATAGCTTTATCGTCAGCTTCCATAAGCGCATAATCTATCACCATTGCTCCATAACGCTCAATGGCATTACTCGACATATCGCATTCAGAAAGAACATTCGCCAGCGGCTAACGTCTGACCAGCCATTTCTCACCG
# +
# :*'4%%$$,3:;<;;<;=3>63($$$+$%'.(*$''(+$#%&#$$#$+7999EA;*&8:776240%$$&##&#$%%%&))&*)$&(&'.(%((55=.--BE9A43)</,,,%,61-/$4#&+)%0&889<74..,*.399?<6?2:3(D=AHENED.*:45**))).+503*6;FA98,%&*&AEMLLC@C@@=>>E?<')064128E>?;/(../.*7/77=?A</,+))%%))7DBCKKE887.37,$$//.)/5$$%;DDM448>CF93().+//-1552-.-/2:=>9::EJDBE6/%4/'8;9<E,(=;853*/$$&&$$$#%#'"#%)&&$%,5B=5550332,)+,&',-$$$+&'3996>;?@<=EFEE::60,'&%%$)'$%&$$&$&%&''%#&'+077=147)'%065.8-=A@::9*/11/)%06$$$$1996;?;45586+%$%)&%&0''577FE?E65;810-*(')%%$*.+%1*0./1-/.213%(%+'&'%#')'&'##&&./%.@61,*%'((+:<AA:,*/95;,977:3%(&233>:4/<(?AA?BF@@>:IHH9B<(&&(+/4++34($8+<).@CDF?@B?>?<><:6//.1A=88*.668:AF73@401146++0+54FEG*@%9@?)))%('((/0)'%(%&*03/;@%0&-32,$#(;/**+-'(#(-;E9$;HD(+06374>>097:;I9.5..)),(%&&-+$$%:51<9IA6&./18;<>0,<0/68&'$$*+%@B@<?.$*3&:(44*%/(520,0033/12'3>1:?AB;/.>7))(('1)%4*46%'1.07;@@;995,8(',,+-7;;>E>/A.))((:7057863,$%%&&***(,5:>@?=;A@312::2*,63532<9@>D?;=403102134;<'9%&%#%'%&($01)*(*'%#%(%"'%+(%%)'&%*+1.&((%$%$$&*&)$.(%))%&&10//%&&&+'%+,33?D><4211&)%(-8;<+4)',0133B?.+++.90.,$&*('--'*--10*'#$+/&):><>GJ/23;><6/17/,)))8;8/+)($&$$$53((&%'%%($'),*)&#&##$'%%%1#((.:@@=>-()/4$$%%%%+&(%%+#$,+,(8.@?E=B?,&&#)$&''$&*-+-5(%%,333/--)-3-7+,/*+%%'-'-)*)%&%)-'&&+**+*'*+$*&$#$&&+.306,+%&'&%&&++&%,1/-2('+(&''-$$)&%%$)&$$#$%%$#$#&(('$#$'(.0#$#%&01'((%$$$%$#'($$&%$&&&#$####'%#*&#)%)%'&&)(('(',,$
#
# from 13, I removed
#
# @0894d2f4-3d3e-4045-9735-41a0bf345728 runid=cc18a835305a98d8cc7f2ff70b00bef64778926f read=11 ch=505 start_time=2019-10-25T15:59:39Z flow_cell_id=FAL22147 protocol_group_id=Exp1_25Oct sample_id=sample3_13 barcode=unclassified
# TCAACTCCAGTGACACATCATACGCACGAATCAGGCCTGCGCAAGACTGCGATTTTCACTTTTTTTACGCTACAATCGTATAAATTTACGGTGAATGGAGAAAAAGCGCTTATGCCTTTCTCGCCAGCGTATTATTTTTAACCAATTTGATTTGCTTGATTAAAAGAAACAGAATAGCGGGCAGCTTGATCAGTTGGCACGTGACTATCACGGCAATCACTATCACCTCTCCGACGTGCAGCCAGTCATCGCTGTCTGAATGCAAAAAATACGAGCAATCCGTTTCGGTCAGAAGATATAGGTACACAGTCGCGTTAACGGTGGGATCGGCAGCGGTAAAGCACGGTCGCCGACGCCTTCGCGCAACTCGGGGTTAACGTTGTTGATGCCGATGTTATCGCCCGCCGGTGGTCGAACACGGCACCCCCAGCTACTCAGGCCATCATGCGCGCCACTTTGGGTCGCAGATGATTGCCCCCGACGGCACGCTGAACCGCCGCCTGCTGCGGGAAAAAATTTTTTGCTCACGTGGAGATAAGCACGGCTTAACGCTGTTGCATCCGCGTGATCAGCGAAACGCGCCGCCAGATGCAGGCTGCTACCTCCCTTATCTCACTGGGTCGTACCTTTGCTGGTCGAAAACCGCCTGTCAGTCAGGCCGATCGCGTACTGGTCGTCGATGTGCCAAAGCGCGGATCGAACGCACCATGCTGCGTGACAAGGTCAGCCGCGAACATGCTGAACATATTCTTGCCGCTCAGGCCACGCGCGAGCGGCCGCCGCGGACGATCATTATTGAAAACACTGGCACGCCGGATGCAGTGGCATCGGATGTTGCCCTTACCACGAAAAAGTATTTAATGCACAGCATCGCAGGCCGCCTCACAGGAAAACTCGTAATGCATACCCCCGTTCTGTTTGAACATCCTCTTAATGAGAAAATGCGTGCACGGCTGCGCATCGAATTTCTGATCAGCAAATGGCTTTCCATCCGATGATTGCCACTCATGCCGATGCGCTCCATTTTTCCGTAACCTTAAACGATCTGTTGGATGTGCTGGAGGCGTGGCGAAGTGCGTACCGACTGGTTAAAGAACTGGGCGCCAGCAGCGCAAACTCCAGTCACGGGCTGAAGTACCGGGCGTCGATCGGGAGGCGCATTAACGAACTGCGCCAGCAGCTGAAGCAATCCTCCAGCACCTGATGGCCGCGCCGCCCGGCTATTGTCCAGTTTCTGCGTGAAAGATCGCCTGATTGCGCTGGTGCGTCAGCGCCACAGATGTCAGGGCGGCTGCTGCGGCGCTGATTTACCGACCCTGCATATCTGGCTCCATATGCCACACAGGCGCATCGGGACGAGCAGGTGGCCAGCTGGCTCGCCAGCCTCGATCCGCTGATCCAGTCCTGAACCTGATCCTCGACCTGATCCGCCAACTCTGCGCTGTTTCGCAAACAGACCAGCCTCAGGCTTTACCCGAGGATAATGGCGAAAGATGCCGATCTGCTACGTCTGCGTCTCGATCTGGCTCACCAGCTGTATCCGCAGATATCCGGCCATAAGAGCCACTGGCTATCGTTTCCTGCCGCTGGACGTGGTAATAGGATCGTGCGAGCCTTTGATTTTGAACTGGCCTGCTGCTAAGGGAGATAACATGTCTGAAGAGACCATCGTCAACTACCCGACCTGCGGCAAAAACGTGGTGTGGGGCGAAGAGCAGAGTCCCTTTCGTGCATTCTGCAGCAAGCGCTGTCAGCGCTGATCGATCTGGGAGAATGGGCCGCAGAAGAAACGCATCCCCAGTGCGGGCGATCTCTCCGACAGCGATGACCAGACGGTAGCAGCCTGATTTACCGGCGCTAACCCGATAATTATCGGGTTAGCGTTATAGTAACCTCAGCAGGCCGATAAGTTTGTTGATAACCAGTTCATTCGCCAGAGGAAACGCCTCTGGATCAGTCACGACCGCCTCGCCGTCCCGGTTGCCCTTCTTTACCCCATGGCTCGCCTTGCCAACTTTCACAAGGAAAAACCACAGCGTTATATGGCGGTCAGGAACTGATGCCAGCTTATCAAATAGCGTTGCGAACCTTACCGTAATTCCAACCTCTTCCTGCAGCTCGCGAATTAGCGCCTGCTCCGGCGTCTCATCAGATTCAATCTTCCGCCGGGAACTCCAGCTTATTCGCCATGGGCATCAGCGCTGGTAATAAAAATTTCGCCTTGCGGATTGCGAATTATCCGACGACGTTGCAGTTTTTCGTTGGATCAAGCTCCCCGCCAAAAAGGCGCAGAACTCTGCGCCTTTGTTTTTATCAGGCCTTAGCTCAGACGGCCATCGGCACTGTTATATTTTTACCGGAACCGCACGGACGGGATCGTTACGCCCCACTTTACGCTCGTAAGTTTTGTGCCGCCGGATCCTCGGCAGCGGCCTCGTCATCGCTCTGATGGCTGAGGCTGCTGCATCTCTTGAGCGCTCGGCTTCTTCACGACGCTGTTTGCTCCATCGCTGACTCTTCCGGCATACGCACCTGGACTTTGCTCAGAGTGCTAATCACTTCATACTTCAGTGGCTCAACATCACCCAGCAAACATGAGAACGATTCGCTTGTACTCCTGCTTCGGATCTTTCTGCGCATAGCCACGCAGATGGATGCCGAGCGCGAATAGTCCATCGCTGCGGGTGCTCTTTCACAGGAGTCGGGGTCTGCAGCGTCACGCCTTTCTCAAAGTGACCATCATCTCGGCACGACCACTTCTTCTTTACGCTGGTAAGTTTCTACCGCGCTCTGCAGAATACGTTCACGCAGGGTCTCTTCATGCAGCTCATTCTTTATCCAGCCACTCTTTGATCGGCAGCTCGAGGTCGAAGTCGTTTTTCAGACGCTCCTGCAGGCGCTAATATCCCACATCTCTTCCAGCGACTGCGGCGGAATATGCGCATCGATAGTCGCTTTAAAGGCAATCTTCGCGGATGCTGTTGATGGTTTTCTTGACATCAGACACGTCCAGCAGTTCGTTACGCTGAGTGTAAGATCGCGCGGCGCTGGTCGTTAACCACATCATCATATTCCAGCAGCTGCTTACGAATATCAAAGTTACGGCTTTCCACTTTACGCTGGGCGTTGGCGATAGCCTTAGTCACCCATGGGTGCTCGATGGCTTCACCCGGTTTCATCCCCGGTTGCGCATCATGCCGGACACACGATCGGAGAGCGAAGGCTGCGCATCAGAGCATCTTCCATCGACAGGTGAAGCGGGAAGAACCGGCATCACCCTGACGACCCGCACGGCCACGCAGCTGGTTATCGATACGGCGTGATTCATGACGCTTTCAGTACCGATAATGTGTAGACCGCCTGCGGCCAGAACCCTCATGGCGCACCTGCCAGTCGCTTTGATTTTTCAATCTGTTCCGCCGTCGGATTTTCCAGCGCGGCGACTTCCGCCTGCCAGCTACCGCCGAGCATAATATCGGTACCACGACCCGCCATGTTGGTGGCAATGGTCACCGCTGCCGGATGCCCGCGGGCGACGATATCCGCTTCGCTGGCGTGGAACGGCGTTCAGTACGTTGTGCCCAATCTCGCTTTGGTCAGCTCGCGCAGGCTACTTCGGATTTTTCGATCGAATAGTCCCCACCAGCATAGAAGCTGGTAAATGGCGGTGCAGGTTTTAATATCTTCGATAATCGCCTGAATTTTTCCACATTTCGGTCATATAGACCAGATCCGCCCATGTCTTACGGATCATCAGACGGTTAGTCGGCACGACAACGGTATCAACTATAGATAGAAGGCTGAATTCGAAGGCTTCGGTATCCGCTGTCATAGGTCATCCCGCCAGTTTCTCGTACGAGCGGAAGTAGTTCTGGAAGGTGATAGGGGCCAGAGTGCCCAGTTTTCGTTCTGGATCTCACGCAGCTACGGCTTCCCACCCGCTGGTGCAGATATCGGACAGCTTGGCGGCCTTGCATCGTACGGCCGAGTGTGCTCATCGACGATGATAACTTCGCCATCTTTAACGATGTAGTCCAACGTCGCGGGTGAACAACGCATGGGCGCAGCGGCGGTCACATGGTGCATCAGCGCTGATGTTGGTCGAGTACAGCGACTCACCTTCTTCCATAATGCCTTCCCTGCCAGCAGCTCTTCGATGAGCACCAGACCGCGTTCGGTCAGGTTAACTGGCGGGCTTCTCATCAACAGAGAGTGGCCTTCACCGGTAAGGTATCCGAATCCTCTTTCCTGACGGATCA
# +
# ,(.)+++;&./-,))'58;+((59431++)+?+/4777:8@MVF8DDHFDB:CEC*4$&-012333,437/-=C>:<,&((-/13;7&##'*)A?@@/0*@D=52*<C==;9,,1(.0+*1<?BB@2,'($(,08>=;06-*'('(/;=9A?7&$$1A@KCGL<&17::=>B?5**,')++'$.&&&0(%>&8';<0-189D@74($%#%+0*'&)5<?C<*&1'+,91@-9=C10)((#'.6444>-);((7)&>6-$$+''$(=8?3000((1AJ;:A;B<9.'%(%8<88:;38.(&2+)-'967>?F:8B;3$$$$+/-2<62?.35(?>>>:D@3)'((5A@>6A>))&(,,5,.26994%,-=.16E<D:79<IK925?:99745$107:413'$30-32-,<;/'?/**/,&(17;NMHA8+/)();*D39=65&)2((%')$)12)3358A<9:@2689??BB=E=;ONM9@G...&&'-''(5:A>19?AAB:C?=,'+.@.B>;9-/BJED+5>BD:<A55-0??CE6..(%3*'''2067>B>+8DA;%%644+'*&(3122&&-2/.%7<<ASNKMF@GB()@BC6'-<6894<44<02.*02262/;6<.?@A554@2+1-72@,,0/%AA>FE<+%;$'&4,:A@E@1+1&%*,,.&,19:.;=E6:;<62420($*,-&&56?%%)779;:5:?:37'$*$-%(&(.<892,EE4$$$9<?====@73-.?@?=BDC9BBDH@E?AIJN<>F0>543(%)563-246:?88.1,)%#-/((*'-46%++-/5677:1IHA@0,#,,4?EC=4/A8;?78=@964###%04776BA-+((2+.19//8*#$'),-769$$$$+JJ:C56:.8E1@1AI?:;<59=<B@D=98GB@4/;JA@B=NLCA@?B@..7:/9@EF<<6&92-865('$%&&#%%+65//..*:0/0,6)<:B)$D996951-.-(457?A@,,B61&(+./+2CCBBHAA0986/53594834=;IKDC=B:>LD$'$$44(((*99A;BBND@;;2A?CAA-/2;?@F@>.+*2=951,-35,)'(':GH76$2,/6:B6<3,DG?983.--/1/*(,(;550/($&/AAC?::4D<KHK<G6:3C9.3:%(9%(0221+)+..*.29;LD@LPED61'3.)*(-+5;:7:AC0IA566HDAIE:;;4A>75:360*)&&$%&&,+#)*(3(;>??53@9>71$.;:>?BJ<A69;?7869,2CE=93H>B=92,''()%%$%%#(+002?<@B=?1--/%(%%1%%(<'2-4589'+213:;;@>B@<;22%(+**)('$'(17A?C@-&<6:E?10-%(),/*'&%%68:50)69A>FF895*(283,,('(*7:..B$>45922,(573662260>?:>;211&&/('(5+'=::?:=GBEDBHE&)++'5459@=;984'1%%(###(,,**,&-<>DG@337;223-0243+.2/:<:.*6&/.5114-058.32%278;C5:LCD?$;335)424++0''(##-*:;752%&,)',000&&'$*.0'$;1AAA=;5599=368;9#$:&.<B1,,,,::<,9.$&+.2,$/8?=87/..2789'0:76@;>;5764%10,299789:;76C?=17'07-;;:=>41GJ?9/+98=?A=*2<7421JC20&$$4%,8.15752>)))&%-10%&59:J11@A6)40:779<7<;,0-54+/0+,+&A32%*7=F@BD>,;;83207+,2;6=B;8*8*-5+$0&/***?%%8234/1.%%'''2,776391522/05/((08'"#$$$%'*&'/10+.,2=967D4@?=1))931(9:9?6@G@<;)0D//-)'),100/#$%$%33((,9(,-'(A@A32:=64203(DAF>9>32;&26??>ACH3<0&'%&&$(633332=??63.?*-2'&$(28<;C>*3$78>%%<6).<,=;:<;?@KKG0510.-))/==7<=;<;.666'///,76<?(AA1BED:>@*%%<>;DD0:;9<AB8617&),#/97>6))6:8.---76:=@4..?=CLL;?43:0.-.**++0BC;A&=:./1;>3(8#$)422/1=9;7;7:>94<9:<;=8;.CH<:3,<7<//B7;8??BD?=@@<0.8,--)+4<9<8EC=ABD>=:HCHJ634:=:99=(%&/573222857@7,(35BCCECBG<;*.,*+00333.:,0152/.1''&41678=?@E>;<%-/565%%$&$/0>G==&%%JI=,%1@=@BA=@<<==@0*&%$##%*000/73/8227@=??'A?*?14:21787((&%(),78A>@2:?A>@?>7=>:877@<00149@<%%)578:6@>52#$$$$134575.1'&4$%'(1<BGA::9-1'(&&+4((.->93-++8365264.8@;=C2:399B@=DC(-849:889%)1//-'$%&&,,-/355?=*$156935/02%'+055654$"%'+*/=32330-3,/<99@.(535-)*1A@CD?4)*<>3:89=9*2?630>-$1466/40*+6;66+,$$&96C<&:831%%&$-)+:9+'*$$%*:96:6@;30<BC935(+):'/565>=?==<;=@>48,/2/562/39<@?11=;=87((78%%-899%$+((7577-69-9B25,083-.9366080009(.2/B34:40>+*0.)$#$$&*5:6.9<5.C?::-'&'$%%=A79<5*&+('%%%(*$()0--9:78F@>?67;><>9GJH:7,?CB<59@<B:<@>439:=?C?67'*6421+4,0<9A8:/*),(.1.*&$'14:85:=AB5HC65BE??220:=D<<3:=1=??D?7?,,,477:JIMB:>10B>ACA@G??A+02*'#$((2:11@-,..>=@$BCJH52*;;:A1+2/''=;4**$$(&.66CD2,,+&*(%'+,85+)$$$&7=@?2.3,?743(()%%%$11*1876&%*<4-,7<&&,**+1''8<,;%;59@5:778?@>>(09811+(754(%,5:?8;FBF2,76>><*2;.:;==(2;+2349$>10/+.0475698+,-2/0*=D<=)23-,.6AHGG3G9==>IC>@>BFE?06LNEB:?:<;:-779:=>92:83/%%)0<5J8<6*'%)**44'36449;=>C<78=>=%'=E5&-107=;440@/33/-,527-'=359%(37-$$)1'%)%#"(-5.-'.22/6-1&/6:3;*$$%%%$%*)..3.;76F=BE@A?><>DD?-7JF27;<<A-30359G68699.74;@AA@.>F=99:7:*)8EIFE26.2>0+*$36.;*&4:8?=?HKBC@%%$&&&8(9??:6221%57:83?=71504485;?@?<)<LM<-$>20079:9BCOE<@A6--"##'#*(4641>-1-:AI<6<'</79<@>>5;87*.23*-9>;591.+@ABAB=:8:54:1.56:+,3/%5616EKD73SM@AD<=,.++0<HC?B=:+49949==GB?:7'$..)&+&1<8>3)-167:;F)456:9>=@>%'5+3/4,-/15G0E5:?GCCD8820/11+(+)'))2128(.78%'754'&%$%&(067HM;D;=88IKD9AC2-1)+)895@CJN98$E;1.&#(#)&$###)+.)&#$*/'-94?**AHHC842074,+,14>@E-.4'''%&%%%;;86(,$$#%(-.*/58%%.,34339,*:**6450)&&%((%(/0,>842/-9)9;--0,.&*&3&)..).><F678=D+&7/%"(.+/13("%$(/10/2/+/47*),+,$$9227?E>AB-=;411&&'.3;?@&&3)'AB;FAJA46C566((',$$*%(/57@IDKFM8/7C44=811-((;6<%)#$("#$$&77@E>::=9,+(+11+-6'%$%#####"$%$1+42'0221&#'&'5*2'&)&)7==?*''(())*0)*7E77<3%+48&7//02$1&#3;<9>=745<=5,*,.&')'*()5;@HLG@<@DFDPC.100*$&)%%AG:34=6-'&))%(.2**%,448?C@>4133-+%$&'**&##$.332054.=-2'''9A;;>DGC93=*$*51/3)&2&893<977;8;>02*,35:99795+/#,(.00'--+%2>DDD>9+++1899:0>/92B1,><=%+0@.581:6)$$%+4*.2./1..530A<&00/5#&-++'+,1426IBIA<;8(5:9$$*88@%%,*<GEB*7-3/0-//)2$
