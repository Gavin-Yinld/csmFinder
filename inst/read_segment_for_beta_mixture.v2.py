def extract_seg_info(coordinate,pattern,cell_id):
    chr_id,loci =  coordinate.split(':')
    each_loci = loci.split('_')
    start = each_loci[0]
    end = each_loci[3]
    length = str(int(end)-int(start)+1)
    cell_info = ""
    cell_pose_info=["","","",""]
    all_seg_pattern = pattern.split(';')[0:-1]
    mc_count_all=0
    tc_count_all=0
    mc_per_cpg = [0,0,0,0]
    tc_per_cpg = 0
    for one_seg_pattern in all_seg_pattern :
        one_pattern, pattern_count = one_seg_pattern.split(':')
        cpg_number = one_pattern.count('1')
        mc_count_all = mc_count_all + cpg_number * int(pattern_count)
        tc_count_all = tc_count_all + 4*int(pattern_count)
        tc_per_cpg = int(pattern_count) + tc_per_cpg
        for i in range(0,4) :
            mc_per_cpg[i]=mc_per_cpg[i] + int(pattern_count) * int(one_pattern[i])
    for i in range(0,4):
        cell_pose_info[i] = '_'.join([each_loci[i],str(mc_per_cpg[i]),str(tc_per_cpg)])
    cell_info = cell_id +':' +'_'.join([str(mc_count_all),str(tc_count_all),'4'])
    cell_pose_info = cell_id + ':' + ':'.join(cell_pose_info)
    seg_info_one_cell = dict(chr_id=chr_id, start=start , end=end, length=length,cell_count=1 ,cell_list=cell_id ,cell_info=cell_info, cell_posi_info=cell_pose_info )
    return(seg_info_one_cell)

def read_segment_for_beta_mixture(path="./",out_path="./",split_by_chr=0):
    import os
    seg_info = {}
    for root, dirs, files in os.walk(path):
        for cell_id in files:
			if 	cell_id.endswith(".segment_filter_.final_") :
				full_file_path = root + '/' + cell_id
				input = open(full_file_path)
				cell_id = cell_id[0:-23]
				for line in input:
				   line = line.strip()
				   coordinate, pattern = line.split('\t')[0:2]
				   seg_info_one_cell = extract_seg_info(coordinate, pattern,cell_id)
				   if seg_info_one_cell["chr_id"] in seg_info:
					   if coordinate in seg_info[seg_info_one_cell["chr_id"]]:
						   seg_info[seg_info_one_cell["chr_id"]][coordinate]["cell_count"] = seg_info[seg_info_one_cell["chr_id"]][coordinate]["cell_count"] +1 
						   seg_info[seg_info_one_cell["chr_id"]][coordinate]["cell_list"] = seg_info[seg_info_one_cell["chr_id"]][coordinate]["cell_list"] +';'+ seg_info_one_cell["cell_list"]
						   seg_info[seg_info_one_cell["chr_id"]][coordinate]["cell_info"] = seg_info[seg_info_one_cell["chr_id"]][coordinate]["cell_info"] +';'+ seg_info_one_cell["cell_info"]
						   seg_info[seg_info_one_cell["chr_id"]][coordinate]["cell_posi_info"] = seg_info[seg_info_one_cell["chr_id"]][coordinate]["cell_posi_info"] +';'+ seg_info_one_cell["cell_posi_info"]
					   else :
						   seg_info[seg_info_one_cell["chr_id"]][coordinate] = seg_info_one_cell
				   else :
					   seg_info[seg_info_one_cell["chr_id"]] = dict(coordinate = seg_info_one_cell)
				input.close()
    if split_by_chr == 1 :
        out= {}
        for k,v in seg_info.items():
            if(k in out ) :
                for k2,v2 in v.items():
                    out[k].write('\t'.join([v2["chr_id"],v2["start"],v2["end"],v2["length"],str(v2["cell_count"]),v2["cell_list"],v2["cell_info"],v2["cell_posi_info"]])+'\n')
            else : 
                out[k] = open(out_path+"/beta.mixture.input."+k+".txt",'w')
                for k2,v2 in v.items():
                    out[k].write('\t'.join([v2["chr_id"],v2["start"],v2["end"],v2["length"],str(v2["cell_count"]),v2["cell_list"],v2["cell_info"],v2["cell_posi_info"]])+'\n')
        for k,v in out.items():
            v.close()        
    else :
        out=open(out_path+"/beta.mixture.input.txt",'w')
        for k,v in seg_info.items():
            for k2,v2 in v.items():
                out.write('\t'.join([v2["chr_id"],v2["start"],v2["end"],v2["length"],str(v2["cell_count"]),v2["cell_list"],v2["cell_info"],v2["cell_posi_info"]])+'\n')
        out.close()       





