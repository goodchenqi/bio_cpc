# -*- coding: utf-8 -*-
"""
Created on Tue Nov 13 10:58:24 2018

@author: Administrator
"""
import random
import time
import sys
import os,re
import multiprocessing
import logging
import math
import argparse
from collections import OrderedDict
import numpy as np
from config import Config

####Description####
usage = '''
Author : chenqi
Email  : chenqi@gooalgene.com
Date   : 2018-11-14 16:16:18
Version: v1.0
Description:
    
Example:
    python %s 
''' % (__file__[__file__.rfind(os.sep) + 1:])

class HelpFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def check_dir(dirs):
    if not os.path.exists(dirs):os.mkdir(dirs)

def check_software(software_path):
    if os.path.exists(software_path):
        logging.debug("Choose software:" + software_path + "!")
    else:
        output = os.popen('which ' + software_path)
        software_temp = output.read().strip()
        if os.path.exists(software_temp):
            software_path = software_temp
            logging.debug("Choose software:" + software_path + "!")
        else:
            logging.error("Can't locate the " + software_path + "!")
            return 0
    return software_path

def run_cmd(cmd):
    logging.info(cmd)
    flag = os.system(cmd)
    if flag != 0:
        logging.error("Command fail: " + cmd)
        sys.exit(2)
    return 0

def fmt_time(spend_time):
    spend_time = int(spend_time)
    day = 24 * 60 * 60
    hour = 60 * 60
    min = 60
    if spend_time < 60:
        return "%ds" % math.ceil(spend_time)
    elif spend_time > day:
        days = divmod(spend_time, day)
        return "%dd%s" % (int(days[0]), fmt_time(days[1]))
    elif spend_time > hour:
        hours = divmod(spend_time, hour)
        return '%dh%s' % (int(hours[0]), fmt_time(hours[1]))
    else:
        mins = divmod(spend_time, min)
        return "%dm%ds" % (int(mins[0]), math.ceil(mins[1]))

###########blastx.feat result#########
def extrace_blastx_features(file_table,file_ref):
    existed_fead=OrderedDict()
    new_existed_fead=OrderedDict()
    curr_query_id,curr_sbjct_id,curr_hit_seq_count,curr_hit_HSP_count,curr_hit_score,curr_hit_frame_scores='','','','',[],[[],[],[]]

    for line in open(file_table,'r'):
        if line.strip().startswith('$') or line.strip().startswith('#'):
            continue
        info=line.strip().split('\t')
        if not curr_query_id:
            curr_query_id=info[0]
            curr_sbjct_id=info[1]
            curr_hit_seq_count=1
            curr_hit_HSP_count=0
            curr_hit_score=[]
            curr_hit_frame_scores=[[],[],[]]

        if curr_query_id!=info[0]:
            for i in range(len(curr_hit_frame_scores)):
                if len(curr_hit_frame_scores[i])==0:
                    curr_hit_frame_scores[i]=[0]
            curr_hit_score=np.mean(curr_hit_score)
            # curr_frame_score=(np.std([len(curr_hit_frame_scores[0]),len(curr_hit_frame_scores[1]),len(curr_hit_frame_scores[2])],ddof=1))**2
            curr_frame_score2=(np.std([np.mean(curr_hit_frame_scores[0]),np.mean(curr_hit_frame_scores[1]),np.mean(curr_hit_frame_scores[2])],ddof=1))**2
            if len(curr_hit_frame_scores[0])==0:
                print(np.mean(curr_hit_frame_scores[0]),np.mean(curr_hit_frame_scores[1]),np.mean(curr_hit_frame_scores[2]))
            if curr_query_id not in existed_fead:existed_fead[curr_query_id]='\t'.join(['999','1:%s'%str(curr_hit_seq_count),'2:%s'%str(curr_hit_score),'3:%s'%str(curr_frame_score2)])
            curr_query_id=info[0]
            curr_sbjct_id=info[1]
            curr_hit_seq_count=1
            curr_hit_HSP_count=0
            curr_hit_score=[]
            curr_hit_frame_scores=[[],[],[]]

        curr_hit_HSP_count+=1
        if curr_sbjct_id!=info[1]:
            curr_hit_seq_count+=1
            curr_sbjct_id=info[1]
        if float(info[10])>0:
            hit_score=-1*math.log(float(info[10]))/math.log(10)
        else:
            hit_score=250
        curr_hit_score.append(hit_score)
        curr_hit_frame_scores[int(info[6])%3].append(hit_score)

    if curr_query_id:
        curr_hit_score=np.mean(curr_hit_score)

        # curr_frame_score=(np.std([len(curr_hit_frame_scores[0]),len(curr_hit_frame_scores[1]),len(curr_hit_frame_scores[2])],ddof=1))**2
        curr_frame_score2=(np.std([np.mean(curr_hit_frame_scores[0]),np.mean(curr_hit_frame_scores[1]),np.mean(curr_hit_frame_scores[2])],ddof=1))**2
        if curr_query_id not in existed_fead:existed_fead[curr_query_id]='\t'.join(['999','1:%s'%str(curr_hit_seq_count),'2:%s'%str(curr_hit_score),'3:%s'%str(curr_frame_score2)])
    for line in open(file_ref,'r'):
        if line.startswith('>'):
            query_id=re.split(r'\s+',line.strip())[0].lstrip('>')
            if query_id not in existed_fead:
                new_existed_fead[query_id]='\t'.join(['999', '1:0', '2:0', '3:0'])
            else:
                new_existed_fead[query_id]=existed_fead[query_id]
    existed_fead.clear()
    return new_existed_fead

###########ff.feat result###########
def extract_framefinder_feats(file):
    ff_result=OrderedDict()
    for line in open(file,'r'):
        if line.startswith('>'):
            info=re.split(r'\s+',line.strip())
            query_id=info[0].lstrip('>')
            sta,end=re.split(r'[(|)|,]',info[2])[1:3]
            length=int(end)-int(sta)+1
            score=info[3].split('=')[1]
            used=info[4].split('=')[1]
            strict=1 if re.split(r'[{|}|,]',info[-2])[2]=='strict' else 0
            if query_id not in ff_result:ff_result[query_id]='\t'.join(['4:%s'%str(length),'5:%s'%score,'6:%s'%used,'7:%s'%str(strict)])+'\n'
    return ff_result

#############lsv_cbind###################
def lsv_cbind(arg_lsv1,arg_lsv2,outfile):
    w=open(outfile,'w')
    for key in arg_lsv1.keys():
        w.write(arg_lsv1[key]+'\t'+arg_lsv2[key])
    w.close()

def predict(arg_seq,file_predict,outfile):
    w=open(outfile,'w')
    query_ids,seq_len,curr_query_id,curr_query_seq_len=[],[],'',0
    for line in open(arg_seq,'r'):
        if line.strip().startswith('$') or line.strip().startswith('#'):continue
        if line.startswith('>'):
            if curr_query_id:
                query_ids.append(curr_query_id)
                seq_len.append(curr_query_seq_len)
            curr_query_id=re.split(r'\s+',line.strip())[0].lstrip('>')
            curr_query_seq_len=0
        else:
            curr_query_seq_len+=len(re.sub(r'[^a-zA-Z]*','',line.strip()))
    query_ids.append(curr_query_id)
    seq_len.append(curr_query_seq_len)
    count=0
    for line in open(file_predict,'r'):
        info=re.split(r'\s+',line.strip())
        if not re.findall(r'^-?\d+',info[0]):continue
        query_id='NA' if count>=len(query_ids) else query_ids[count]
        query_seq_len='NA' if count>=len(seq_len) else seq_len[count]
        class_type='noncoding'

        if float(info[1])>0:class_type='coding'
        # print('\t'.join([query_id,str(query_seq_len),class_type,info[1]])+'\n')
        w.write('\t'.join([query_id,str(query_seq_len),class_type,info[1]])+'\n')
        count+=1
    w.close()

#######must keep fasta file each sequence same length###########
def trans_file_get(merged_gtf,ref_fasta,gene_list,outdir):
    gtf_file=open(merged_gtf,'r')
    gene_file=open(gene_list,'r')
    gtf=gtf_file.readline()
    gene=gene_file.readline()
    dict_result,dict_id={},{}
    transcript_id=''
    w=open('%s/transcript.fa'%outdir,'w')
    for line in open(ref_fasta,'r'):
        line=line.strip()
        if line.startswith('>'):
            key=line.replace(">","")
            dict_result[key]=[]
            continue
        dict_result[key].append(line)
    each_seq=len(dict_result[key][0])
    for value in dict_result.values():
        if len(value[-1])<each_seq:
            for i in range(each_seq):
                if len(value[-1])==each_seq:
                    break
                value[-1]+=' '

    for line in open(gene_list,'r'):
        info=line.strip().split('\t')
        if info[1] not in dict_id:dict_id[info[1]]=0
    list_seq=[]

    for line in open(merged_gtf,'r'):
        info=line.strip().split('\t')
        match=re.match(r'(.*)transcript_id "(.*)"; exon(.*)',info[8])
        if transcript_id!=match.group(2):
            if list_seq:
                w.write('>%s\n'%transcript_id+''.join(list_seq)+'\n')
            transcript_id=match.group(2)
            list_seq=[]
        if transcript_id in dict_id:
            sta_index=int((int(info[3])-1)/each_seq)
            end_index=int((int(info[4])-1)/each_seq)+1
            sta_site=(int(info[3])-1)%each_seq if (int(info[3])-1)%each_seq>=0 else 0
            end_site=(int(info[4])-1)%each_seq-each_seq
            seq_transcript=''.join(dict_result[info[0]][sta_index:end_index])[sta_site:end_site]
            list_seq.append(seq_transcript)
    if list_seq:
        w.write('>%s\n'%transcript_id+''.join(list_seq)+'\n')
    w.close()
    list_seq=[]
    dict_result.clear()
    dict_id.clear()

def getopt():
    parser = argparse.ArgumentParser(
        formatter_class=HelpFormatter, description=usage)
    parser.add_argument(
        '-i',help='input fasta file[must keep each seq line same length]',dest='input',type=str)
    parser.add_argument(
        '-g',help='input gtf list file',type=str,dest='gtf_list')
    parser.add_argument(
        '-o',help='ouput dir',type=str,dest='outdir',default='.')
    parser.add_argument(
        '--fasta_format',help='if input fasta file not ncbi format,then use it[y|n]',default='n',type=str,dest='format')
    parser.add_argument(
        '--process',help='set the num of cuffmerge process',dest='process',type=int,default=10)
    parser.add_argument(
        '--cpu',help='set the num of blastx threads',dest='cpu',type=int,default=10)
    parser.add_argument(
        '--resume',help='resume file',dest='resume',type=str,default='novel_trans.resume')    
    args = parser.parse_args()
    if not args.input:
        print('fasta file must be given!!!')
        sys.exit(1)
    elif not args.gtf_list:
        print('gtf list file must be given!!!')
        sys.exit(1)
    return args

def noval_transcript_stat(trans_file,cod_file,outdir):
    total_novel_transcript_num,coding_transcript_num,noncoding_transcript_num,novel_isoform_num,novel_gene_num=0,0,0,0,0
    dict_cod={}
    for line in open(cod_file,'r'):
        if line.strip() not in dict_cod:dict_cod[line.strip()]=0
        coding_transcript_num+=1
    for line in open(trans_file,'r'):
        info=line.strip().split('\t')
        total_novel_transcript_num+=1
        if info[0]=='-' and info[1] in dict_cod:
            novel_gene_num+=1
        elif info[1] in dict_cod and info[0]!='-':
            novel_isoform_num+=1
    noncoding_transcript_num = total_novel_transcript_num - coding_transcript_num
    w=open('%s/NovelTranscriptSummary.xls'%outdir,'w')
    w.write('total_novel_transcript_num\tcoding_transcript_num\tnoncoding_transcript_num\tnovel_isoform_num\tnovel_gene_num\n%d\t%d\t%d\t%d\t%d'%(total_novel_transcript_num,coding_transcript_num,noncoding_transcript_num,novel_isoform_num,novel_gene_num))
    w.close()

def transe_fasta_format(file,outdir,samtools):
    for line in open(file,'r'):
        if line.startswith('>'):
            line=re.split(r'\s+',line.strip())[0].lstrip('>')
            cmd='%s faidx %s %s >>%s/novel_trans.fasta'%(samtools,file,line,outdir)
            run_cmd(cmd)
    return '%s/novel_trans.fasta'%(outdir)
    
def main():
    args=getopt()
    config=Config()
    #########check sortware#######################
    cuffmerge=check_software('cuffmerge') if check_software('cuffmerge') else config.cuffmerge
    cuffcompare=check_software('cuffcompare') if check_software('cuffcompare') else config.cuffcompare
    cufflinks=check_software('cufflinks') if check_software('cufflinks') else config.cufflinks
    gtf_to_sam=check_software('gtf_to_sam') if check_software('gtf_to_sam') else config.gtf_to_sam
    sort=check_software('sort') if check_software('sort') else config.sort
    blastx=check_software('blastx') if check_software('blastx') else config.blastx
    framefinder=check_software('framefinder') if check_software('framefinder') else config.framefinder
    svm_scale=check_software('svm-scale') if check_software('svm-scale') else config.svm_scale
    svm_predict2=check_software('svm-predict2') if check_software('svm-predict2') else config.svm_predict2
    samtools=check_software('samtools') if check_software('samtools') else config.samtools
    outdir=args.outdir
    check_dir(outdir)
    resume_num=0
    if os.path.exists('%s/%s'%(outdir,args.resume)):
        for line in open('%s/%s'%(outdir,args.resume)):
            resume_num=int(line.strip()) if int(line.strip())>=resume_num else resume_num
        w=open('%s/%s'%(outdir,args.resume),'w')
    else:
        w=open('%s/%s'%(outdir,args.resume),'w')
    args.input=transe_fasta_format(args.input,outdir,samtools) if args.format=='y' else args.input
    ##############merge stringtie####################
    count_resume=0
    cmd='%s -o %s -p %s -s %s %s'%(cuffmerge,args.outdir,args.process,args.input,args.gtf_list)
    count_resume+=1
    w.write('%d\n'%count_resume)
    if count_resume>=resume_num:run_cmd(cmd)
    count=0
    for line in open('%s/logs/run.log'%outdir):
        if not count:
            count+=1
            continue
        if count==1:check_dir(outdir+'/tmp')
        count_resume+=1
        w.write('%d\n'%count_resume)
        if count_resume>=resume_num:run_cmd(line.strip())
    cmd='awk \'{if($3=="j"||$3=="u"||$3=="o"||$3=="i") print $2"\t"$5}\' %s/tmp_meta_asm.merged.gtf.tmap > %s/gene_transcript.list'%(outdir,outdir)
    count_resume+=1
    w.write('%d\n'%count_resume)
    if count_resume>=resume_num:run_cmd(cmd)
    ###############transcript file get###########
    trans_file_get('%s/merged.gtf'%outdir,args.input,'%s/gene_transcript.list'%outdir,outdir)
    ##############predict#######################
    predict_outdir=outdir+'/cpc'
    check_dir(predict_outdir)
    #######step1:blastx######
    cmd='%s -query %s/transcript.fa -evalue 1e-10 -ungapped -threshold 14 -num_threads %s -comp_based_stats F -db %s -outfmt 6 > %s/blastx.table'%(blastx,outdir,args.cpu,config.nr_database,predict_outdir)
    count_resume+=1
    w.write('%d\n'%count_resume)
    if count_resume>=resume_num:run_cmd(cmd)
    ######step2:framefinder#####
    cmd='cat %s/transcript.fa |%s -r False -w %s /dev/stdin > %s/ff.fa1'%(outdir,framefinder,config.fram_model,predict_outdir)
    count_resume+=1
    w.write('%d\n'%count_resume)
    if count_resume>=resume_num:run_cmd(cmd)
    ######step3:features stat#######
    blastx_result=extrace_blastx_features('%s/blastx.table'%predict_outdir,'%s/transcript.fa'%outdir)
    ff_result=extract_framefinder_feats('%s/ff.fa1'%predict_outdir)
    lsv_cbind(blastx_result,ff_result,'%s/test.lsv'%predict_outdir)
    ######step4:svm-scale######
    cmd='%s -r %s %s/test.lsv > %s/test.lsv.scaled'%(svm_scale,config.lib_range,predict_outdir,predict_outdir)
    count_resume+=1
    w.write('%d\n'%count_resume)
    if count_resume>=resume_num:run_cmd(cmd)
    ######step5:svm-predict####
    cmd='%s %s/test.lsv.scaled %s %s/test.svm0.predict> %s/test.svm0.stdout 2> %s/test.svm0.stderr'%(svm_predict2,predict_outdir,config.lib_model,predict_outdir,predict_outdir,predict_outdir)
    count_resume+=1
    w.write('%d\n'%count_resume)
    if count_resume>=resume_num:run_cmd(cmd)
    #########result###########
    predict('%s/transcript.fa'%outdir,'%s/test.svm0.predict'%predict_outdir,'%s/result.txt'%outdir)
    ######################get new fasta file###############
    cmd='awk \'{if($3=="coding") print $1}\' %s/result.txt |while read i;do echo "%s faidx %s/transcript.fa ${i} >> %s/coding_transcript.fa" ;done |sh'%(outdir,samtools,outdir,outdir)
    count_resume+=1
    w.write('%d\n'%count_resume)
    if count_resume>=resume_num:run_cmd(cmd)
    cmd='cat %s/coding_transcript.fa %s>%s/%s'%(outdir,args.input,outdir,args.input+'.new.fa')
    count_resume+=1
    w.write('%d\n'%count_resume)
    if count_resume>=resume_num:run_cmd(cmd)
    cmd='awk \'{if($3=="coding") print $1}\' %s/result.txt > %s/coding_transcript.list'%(outdir,outdir)
    count_resume+=1
    w.write('%d\n'%count_resume)
    if count_resume>=resume_num:run_cmd(cmd)
    w.close()
    ###################final stat#################
    noval_transcript_stat('%s/gene_transcript.list'%outdir,'%s/coding_transcript.list'%outdir,outdir)

if __name__=='__main__':
    try:
        t1 = time.time()
        time1 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t1))
        print('Start at : ' + time1)

        main()
        t2 = time.time()
        time2 = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(t2))
        print('End at : ' + time2)
        t3 = t2 - t1
        print('Spend time: ' + fmt_time(t3))
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) See you!\n")
        sys.exit(0)
