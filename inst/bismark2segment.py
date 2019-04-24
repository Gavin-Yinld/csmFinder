
def split_bismark_file(file):
    from  itertools import islice
    f = open(file)
    file_list = {}
    for line in islice(f, 1, None):
        sub_line=line.split()
        chr = sub_line[2]
        if chr in file_list :
            file_list[chr].write(line)
        elif len(chr)<=5 :
            file_list[chr]=open(file +"."+chr+".bis_",'w')
            file_list[chr].write(line)
    for k,v in file_list.items():
        v.close()

def handle_bismark(path="./") :
    import os
    for root, dirs, files in os.walk(path):
        for file in files :
            if file.endswith(".bis_") :
                read = {}
                input = open(file)
                for line in input:
                    line=line.strip()
                    [read_id,meth_call,chr,loci,meth_state]=line.split('\t')
                    if read_id in read:
                        read[read_id][int(loci)] = [chr,loci,meth_state]
                    else :
                        read[read_id] = {}
                        read[read_id][int(loci)] = [chr,loci,meth_state]
                input.close()
                segment = {}
                for key,value in read.items() :
                    if len(value) >= 4 :
                        loci = sorted(value.keys())
                        for i in range(0,len(loci)-3) :
                            seg = value[loci[i]][0] + ':' + value[loci[i]][1] + '_' + value[loci[i+1]][1] + '_' + value[loci[i+2]][1] + '_' + value[loci[i+3]][1]
                            pattern = value[loci[i]][2]+value[loci[i+1]][2]+value[loci[i+2]][2]+value[loci[i+3]][2]
                            pattern = pattern.replace('z','0')
                            pattern = pattern.replace('Z','1')
                            if seg in segment :
                                if pattern in segment[seg] : 
                                    segment[seg][pattern] = segment[seg][pattern] + 1
                                else : segment[seg][pattern] = 1
                            else :
                                segment[seg] = {}
                                segment[seg][pattern] = 1
                del read
                out_file = open(file+".segment_",'w')
                for key,value in segment.items():
                    out_file.write(key+'\t')
                    for k2,v2 in value.items():
                        out_file.write(k2+':'+str(v2)+';')
                    out_file.write("\n")
                out_file.close()
                del segment

def process_segment(CpG_file,path="./"):
    import pandas as pd
    import os
    import re
    #cpg = pd.read_csv('C:/Users/Yin/Desktop/test.panda/mm10_CpG+',sep="\t", header = None)
    cpg = pd.read_csv(CpG_file,sep="\t", header = None)
    cpg.columns = ["chrom","loci"]
    cpg = cpg.sort_values(["chrom","loci"], ascending=[True, True])
    cpg["loci"] = cpg["loci"].apply(str)  
    chr_patt = re.compile(r'(chr[0-9]{1,2})')
    for root, dirs, files in os.walk(path):
        for filename in files :
            if filename.endswith(".bis_.segment_") :
                segment = pd.read_csv(filename,sep="\t", header = None)
                segment.columns = ["chrom_loci","pattern"]
                loci = segment["chrom_loci"].str.split("_",expand=True)[0]
                loci_mat = segment["chrom_loci"].str.split(":",expand=True)[1].str.split("_",expand=True)
                search = chr_patt.search(filename)
                if search : 
                    one_chr = search.group()
                    temp_cpg = cpg[cpg["chrom"]==one_chr]
                    cpg_index = temp_cpg["chrom"] + ':' +temp_cpg["loci"]
                    rev_cpg_index = pd.Series(list(set(loci) - (set(cpg_index))))
                    index = loci.isin(rev_cpg_index)
                    segment.ix[index,0] = one_chr + ':' + (loci_mat.ix[index,0].apply(int)-1).apply(str) +'_' + (loci_mat.ix[index,1].apply(int)-1).apply(str) +'_' + (loci_mat.ix[index,2].apply(int)-1).apply(str) +'_' + (loci_mat.ix[index,3].apply(int)-1).apply(str)               
                    temp_len = len(cpg_index)
                    temp_cpg.index = range(0,temp_len)
                    loci_1 = temp_cpg.ix[0:(temp_len-4),1]
                    loci_2 = temp_cpg.ix[1:(temp_len-3),1]
                    loci_3 = temp_cpg.ix[2:(temp_len-2),1]   
                    loci_4 = temp_cpg.ix[3:,1]
                    loci_1.index = range(0,temp_len-3)
                    loci_2.index = range(0,temp_len-3)
                    loci_3.index = range(0,temp_len-3)
                    loci_4.index = range(0,temp_len-3)
                    ref_4CG = one_chr+':'+loci_1+'_'+loci_2+'_'+loci_3+'_'+loci_4
                    segment = segment[segment["chrom_loci"].isin(ref_4CG)]
                    segment.to_csv(filename+"filter_", sep='\t', encoding='utf-8',header=False,index=False)


def merge_segment(path="./"):
    import os
    for root, dirs, files in os.walk(path):
        for cell_id in files:
            if cell_id.endswith("filter_") :
                full_file_path = root + '/' + cell_id
                input = open(full_file_path)
                segment = {}
                for line in input:
					line = line.strip()
					coordinate, pattern = line.split('\t')[0:2]
					if coordinate not in segment : segment[coordinate] = {}
					pattern_list = pattern.split(';')[0:-1]
					for one_patt in pattern_list :
						tmp=one_patt.split(':')
						patt=tmp[0]
						cov=int(tmp[1])
						if patt in segment[coordinate] : 
							segment[coordinate][patt] = segment[coordinate][patt] + cov
						else : segment[coordinate][patt] = cov
                out_file = root + '/' + cell_id + ".final_"
                out=open(out_file,"w")
                for k in sorted(segment.keys()):
					out.write(k + "\t")
					for k2 in sorted(segment[k].keys()) : 
						out.write(k2 + ':' + str(segment[k][k2]) + ';')
					out.write("\n")
                out.close()


def handle_bismark_single(file):
    from  itertools import islice
    f = open(file)
    read = {}
    for line in islice(f, 1, None):
        line=line.strip()
        [read_id,meth_call,chr,loci,meth_state]=line.split('\t')
        if read_id in read:
            read[read_id][int(loci)] = [chr,loci,meth_state]
        else :
            read[read_id] = {}
            read[read_id][int(loci)] = [chr,loci,meth_state]
    f.close()
    segment = {}
    for key,value in read.items() :
        if len(value) >= 4 :
            loci = sorted(value.keys())
            for i in range(0,len(loci)-3) :
                seg = value[loci[i]][0] + ':' + value[loci[i]][1] + '_' + value[loci[i+1]][1] + '_' + value[loci[i+2]][1] + '_' + value[loci[i+3]][1]
                pattern = value[loci[i]][2]+value[loci[i+1]][2]+value[loci[i+2]][2]+value[loci[i+3]][2]
                pattern = pattern.replace('z','0')
                pattern = pattern.replace('Z','1')
                if seg in segment :
                    if pattern in segment[seg] : 
                        segment[seg][pattern] = segment[seg][pattern] + 1
                    else : segment[seg][pattern] = 1
                else :
                    segment[seg] = {}
                    segment[seg][pattern] = 1
    del read
    out_file = open(file+".segment_",'w')
    for key,value in segment.items():
        out_file.write(key+'\t')
        for k2,v2 in value.items():
            out_file.write(k2+':'+str(v2)+';')
        out_file.write("\n")
    out_file.close()
    del segment

def process_segment_single(CpG_file,path="./"):
    import pandas as pd
    import os
    cpg = pd.read_csv(CpG_file,sep="\t", header = None)
    cpg.columns = ["chrom","loci"]
    cpg = cpg.sort_values(["chrom","loci"], ascending=[True, True])
    cpg["loci"] = cpg["loci"].apply(str)  
    ref_4CG = pd.Series([])
    chrom = cpg["chrom"].unique()
    for one_chr in chrom:
        temp_cpg = cpg[cpg["chrom"]==one_chr]
        temp_len = len(temp_cpg)
        temp_cpg.index = range(0,temp_len)
        loci_1 = temp_cpg.ix[0:(temp_len-4),1]
        loci_2 = temp_cpg.ix[1:(temp_len-3),1]
        loci_3 = temp_cpg.ix[2:(temp_len-2),1]   
        loci_4 = temp_cpg.ix[3:,1]
        loci_1.index = range(0,temp_len-3)
        loci_2.index = range(0,temp_len-3)
        loci_3.index = range(0,temp_len-3)
        loci_4.index = range(0,temp_len-3)
        ref_4CG = ref_4CG.append(one_chr+':'+loci_1+'_'+loci_2+'_'+loci_3+'_'+loci_4)
    for root, dirs, files in os.walk(path):
        for filename in files :
            if filename.endswith(".segment_") :
                segment = pd.read_csv(filename,sep="\t", header = None)
                segment.columns = ["chrom_loci","pattern"]
                loci = segment["chrom_loci"].str.split("_",expand=True)[0]
                loci_mat = segment["chrom_loci"].str.split(":",expand=True)[1].str.split("_",expand=True) 
                cpg_index = cpg["chrom"] + ':' +cpg["loci"]
                rev_cpg_index = pd.Series(list(set(loci) - (set(cpg_index))))
                index = loci.isin(rev_cpg_index)
                segment.ix[index,0] = segment["chrom_loci"].str.split(":",expand=True)[0] + ':' + (loci_mat.ix[index,0].apply(int)-1).apply(str) +'_' + (loci_mat.ix[index,1].apply(int)-1).apply(str) +'_' + (loci_mat.ix[index,2].apply(int)-1).apply(str) +'_' + (loci_mat.ix[index,3].apply(int)-1).apply(str)               
                segment = segment[segment["chrom_loci"].isin(ref_4CG)]
                segment.to_csv(filename+"filter_", sep='\t', encoding='utf-8',header=False,index=False)
