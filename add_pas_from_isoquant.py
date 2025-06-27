import argparse
import sys
from datetime import datetime
import os
import gzip

def format_datetime():
    now = datetime.now()
    formatted_time = now.strftime('%Y-%m-%d %H:%M:%S')
    return formatted_time

def get_isoquant_fsm(input):
    path,out = input[0],input[1]
    with gzip.open(path,'rt') as file,gzip.open(out,'wt') as output:
        for line in file:
            if line[0] == '#':
                continue
            info = line.strip().split('\t')
            asstye = info[5]
            assevent = info[6]
            if asstye == 'unique':
                if 'fsm' in assevent.split(','):
                    output.write(line)

def get_isoquant_mem(input):
    path,out = input[0],input[1]
    with gzip.open(path,'rt') as file,gzip.open(out,'wt') as output:
        for line in file:
            if line[0] == '#':
                continue
            info = line.strip().split('\t')
            asstye = info[5]
            assevent = info[6]
            if asstye == 'unique':
                if 'mono_exon_match' in assevent.split(','):
                    output.write(line)

def filter_bed12(input_list):
    bed12,full_length,output = input_list[0],input_list[1],input_list[2]
    with gzip.open(bed12,'rt') as bed12, gzip.open(full_length,'rt') as fl, gzip.open(output,'wt') as out:
        fl_dict = {}
        for line in fl:
            read, isoform, gene = line.strip().split('\t')[0], line.strip().split('\t')[3], line.strip().split('\t')[4]
            fl_dict[read] = gene + '|' + isoform
        for line in bed12:
            info = line.strip().split('\t')
            read = info[3]           
            if read in fl_dict:
                info[4] = fl_dict[read]
                out_line = '\t'.join(info)
                out.write(out_line+'\n')

def get_junctions_bed12(info):
    chrstart = int(info[1])
    starts = [int(n) + chrstart + 1 for n in info[11].split(',')]
    sizes = [int(n) - 1 for n in info[10].split(',')]
    if len(starts) == 1:
        return 'mono_exon'
    junctions = []
    for b in range(len(starts) - 1):
        junctions.append(str(starts[b] + sizes[b]) + '-' + str(starts[b + 1]))
    out = ','.join(junctions)
    return out

def get_pas_site(tes_counts,width = 50): 
    peaks_dict = {}
    while len(tes_counts) > 0:
        value_max = max(list(tes_counts.values()))
        key_max = [key for key, val in tes_counts.items() if val == value_max][0]
        left = key_max
        right = key_max
        peaks_dict[((left,right),key_max)] = tes_counts.pop(key_max)

        key_list = list(tes_counts.keys())
        for key in key_list:
            tmp_store = peaks_dict[((left,right),key_max)]
            if key_max - width <= key <= key_max + width:
                peaks_dict.pop(((left,right),key_max))
                if key < key_max:
                    left = key
                else:
                    right = key
                peaks_dict[((left,right),key_max)] = tmp_store + tes_counts.pop(key)
    sorted_by_key = dict(sorted(peaks_dict.items()))   
    return(sorted_by_key)

def is_overlap(flank1,flank2,summit1,summit2,width=50):
    limit1 = (summit1-width,summit1+width)
    limit2 = (summit2-width,summit2+width)
    if summit1 > summit2:
        limit1,limit2 = limit2,limit1
        flank1,flank2 = flank2,flank1
    if limit1[1] > limit2[0]:
        return((True,(flank1[0],flank2[1])))
    else:
        return((False,'_'))   
def merge_overlap(peaks_dict,width=50):
    merge_peaks_dict = {}
    data = list(peaks_dict.keys())
    flank_list = []
    summit_dict = {}

    for flank, summit in data:
        flank_list.append(flank)
        summit_dict[flank] = summit

    while len(flank_list) > 0:
        flank = flank_list[0]

        if len(flank_list) == 1:
            merge_peaks_dict[flank,summit_dict[flank]] = peaks_dict[flank,summit_dict[flank]]
            break

        flank_list.pop(0)
        total_num = len(flank_list)
        n = 0

        for rest in flank_list:
            n += 1

            overlap, new_flank = is_overlap(flank,rest,summit_dict[flank],summit_dict[rest],width=width)
            if overlap:
                flank_list.remove(rest)
                if peaks_dict[flank,summit_dict[flank]] >= peaks_dict[rest,summit_dict[rest]]:
                    new_summit = summit_dict[flank]
                else:
                    new_summit = summit_dict[rest]

                peaks_dict[new_flank,new_summit] = peaks_dict[flank,summit_dict[flank]] + peaks_dict[rest,summit_dict[rest]]
                flank_list.insert(0,new_flank)
                summit_dict[new_flank] = new_summit
                break

            if n == total_num:
                merge_peaks_dict[flank,summit_dict[flank]] = peaks_dict[flank,summit_dict[flank]]
    return(merge_peaks_dict)

def get_pas(bed12,width):
    with gzip.open(bed12,'rt') as bed:
        pas_dict = {}
        for line in bed:
            info = line.strip().split('\t')
            start,end,read,gene,isoform,strand = info[1],info[2],info[3],info[4].split('|')[0],info[4].split('|')[1],info[5]
            if strand == '+':
                tes = end
            else:
                tes = int(start) + 1
            sjc = get_junctions_bed12(info)
            key = gene+'|'+strand
            if key not in pas_dict:
                pas_dict[key] = {}
            if sjc not in pas_dict[key]:
                pas_dict[key][sjc] = []
            pas_dict[key][sjc].append(tes)
    for gene in pas_dict:
        for sjc in pas_dict[gene]:
            tes_list = pas_dict[gene][sjc]
            tes_counts = {}
            for tes in tes_list:
                if int(tes) not in tes_counts:
                    tes_counts[int(tes)] = 0
                tes_counts[int(tes)] += 1
            pas_dict[gene][sjc] = get_pas_site(tes_counts,width = width)
            pas_dict[gene][sjc] = merge_overlap(pas_dict[gene][sjc],width = width)
    import copy
    pas_dict_filter = copy.deepcopy(pas_dict)
    gene_list = list(pas_dict_filter.keys())      
    for gene in gene_list:
        sjc_list = list(pas_dict_filter[gene].keys())
        for sjc in sjc_list:
            pas_list = list(pas_dict_filter[gene][sjc].keys())
            pas_sum = sum(list(pas_dict_filter[gene][sjc].values()))
            for pas in pas_list:
                if pas_dict_filter[gene][sjc][pas] < 30:
                    pas_dict_filter[gene][sjc].pop(pas)
                if len(pas_dict_filter[gene][sjc]) < 1:
                    pas_dict_filter[gene].pop(sjc)
                if len(pas_dict_filter[gene]) < 1:
                    pas_dict_filter.pop(gene)
    return pas_dict_filter

def get_junctions_bed12_v2(info):
    chrstart = int(info[1])
    starts = [int(n) + chrstart + 1 for n in info[11].split(',')[:-1]]
    sizes = [int(n) - 1 for n in info[10].split(',')[:-1]]
    if len(starts) == 1:
        return 'mono_exon'
    junctions = []
    for b in range(len(starts) - 1):
        junctions.append(str(starts[b] + sizes[b]) + '-' + str(starts[b + 1]))
    out = ','.join(junctions)
    return out

def main():

    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-1','--annotation_bed',required = True, help="bed file from isoquant output extended_annotation.gtf")
    parser.add_argument('-2','--reads_assignments',required = True, help="isoquant output reads_assignments.tsv.gz")
    parser.add_argument('-3','--corrected_reads',required = True, help="isoquant output corrected_reads.bed.gz")
    parser.add_argument('-w','--width',type=int,default=10,required = False, help="PAS width")
    parser.add_argument('-o','--output_path',required = True, help="output path")
    args = parser.parse_args(sys.argv[1:])

    time = format_datetime()
    print("----------------------------------------------------------------------------------------------------------------------------------------------------------")
    print(f'[{time}] python' ,' '.join(sys.argv), sep = '\t')

    time = format_datetime()
    print(f'[{time}] Progressing isoquant outputs...')
    if not os.path.exists(args.output_path+'/reads_assignments_fsm.tsv.gz'):
        get_isoquant_fsm([args.reads_assignments,args.output_path+'/reads_assignments_fsm.tsv.gz'])
    else:
        print('reads_assignments_fsm.tsv.gz existed')
    if not os.path.exists(args.output_path+'/reads_assignments_mem.tsv.gz'):
        get_isoquant_mem([args.reads_assignments,args.output_path+'/reads_assignments_mem.tsv.gz'])
    else:
        print('reads_assignments_mem.tsv.gz existed')
    if not os.path.exists(args.output_path+'/corrected_reads_fsm.bed.gz'):
        filter_bed12([args.corrected_reads,args.output_path+'/reads_assignments_fsm.tsv.gz',args.output_path+'/corrected_reads_fsm.bed.gz'])
    else:
        print('corrected_reads_fsm.bed.gz existed')
    if not os.path.exists(args.output_path+'/corrected_reads_mem.bed.gz'):
        filter_bed12([args.corrected_reads,args.output_path+'/reads_assignments_mem.tsv.gz',args.output_path+'/corrected_reads_mem.bed.gz'])
    else:
        print('corrected_reads_mem.bed.gz existed')

    time = format_datetime() 
    print(f'[{time}] Identifying PAS... width = {str(args.width)}')
    pas_dict_fsm = get_pas(args.output_path+'/corrected_reads_fsm.bed.gz',args.width)
    pas_dict_mem = get_pas(args.output_path+'/corrected_reads_mem.bed.gz',args.width)

    time = format_datetime()
    print(f'[{time}] Adding PAS information...')

    with open(args.annotation_bed,'r') as bed:
        annot_dict = {}
        for line in bed:
            info = line.strip().split('\t')
            gene,isoform,strand = info[3].split('|')[1],info[3].split('|')[0],info[5]
            info[3] = isoform
            info[4] = gene
            line = '\t'.join(info)+'\n'
            sjc = get_junctions_bed12_v2(info)
            key = gene+'|'+strand
            if key not in annot_dict:
                annot_dict[key] = {}
            if sjc not in annot_dict[key]:
                annot_dict[key][sjc] = {}
                annot_dict[key][sjc]['isoform'] = {}
            if isoform not in annot_dict[key][sjc]['isoform']:
                annot_dict[key][sjc]['isoform'][isoform] = line

    for gene in pas_dict_fsm:
        for sjc in pas_dict_fsm[gene]:
            if gene in annot_dict:
                if sjc in annot_dict[gene]:
                    annot_dict[gene][sjc]['pas'] = []  
                    for pas in pas_dict_fsm[gene][sjc]:
                        annot_dict[gene][sjc]['pas'].append(int(str(pas).strip(')').split(',')[-1]))
    for gene in pas_dict_mem:
        for sjc in pas_dict_mem[gene]:
            if gene in annot_dict:
                if sjc in annot_dict[gene]:
                    annot_dict[gene][sjc]['pas'] = []  
                    for pas in pas_dict_mem[gene][sjc]:
                        annot_dict[gene][sjc]['pas'].append(int(str(pas).strip(')').split(',')[-1]))
    merge_num = 0
    for gene in annot_dict:
        for sjc in annot_dict[gene]:      
            isoform_list = list(annot_dict[gene][sjc]['isoform'].keys())
            m = len(isoform_list) 
            tss_list = [] 
            for isoform in isoform_list:
                line = annot_dict[gene][sjc]['isoform'][isoform]
                info = line.strip().split('\t')
                info[8] = str(m)
                annot_dict[gene][sjc]['isoform'][isoform] = '\t'.join(info)
                strand,start,end = info[5],info[1],info[2]
                if strand == '+':
                    tss_list.append(int(start))
                else:
                    tss_list.append(int(end))
            if  m > 1:
                merge_num +=1
                if strand == '+':
                    select =  tss_list.index(min(tss_list))    
                else:
                    select =  tss_list.index(max(tss_list))
                isoform_list.pop(select)
                for isoform in isoform_list:
                    annot_dict[gene][sjc]['isoform'].pop(isoform)
    
    with open(args.output_path+'/added_pas.bed', 'w') as file:
        for gene in annot_dict:
            sjc_num = len(annot_dict[gene])
            sjc_no = 0
            for sjc in annot_dict[gene]:
                sjc_no += 1
                if 'pas' not in annot_dict[gene][sjc]:
                    for id in annot_dict[gene][sjc]['isoform']:
                        info = annot_dict[gene][sjc]['isoform'][id].strip().split('\t')
                        new_isoform = gene+'|'+str(sjc_num)+'.'+str(sjc_no)+'.'+'0.0'+'|'+id
                        info.append('#')
                        info[3] = new_isoform
                        line = '\t'.join(info)+'\n'
                        file.write(line)
                else:
                    for id in annot_dict[gene][sjc]['isoform']:
                        info = annot_dict[gene][sjc]['isoform'][id].strip().split('\t')
                        isoform = info[3]
                        m = info[8]
                        pas_no = 0
                        pas_num = len(annot_dict[gene][sjc]['pas'])
                        
                        start,end,strand = int(info[1]),int(info[2]),info[5]
                        if strand == '+':
                            tes = end
                        else:
                            tes = start

                        no_old_ref = True
                        write_old_ref = False
                        for pas in annot_dict[gene][sjc]['pas']:
                            if abs(int(pas) - int(tes)) <= 100:
                                if not write_old_ref:
                                    new_isoform = gene+'|'+str(sjc_num)+'.'+str(sjc_no)+'.'+str(pas_num)+'.'+'0'+'|'+id
                                    import copy
                                    info_copy = copy.deepcopy(info)
                                    info_copy.append(str(pas))
                                    info_copy[3] = new_isoform
                                    info_copy[-1] = '#'
                                    line = '\t'.join(info_copy)+'\n'
                                    file.write(line)
                                    write_old_ref = True
                                    no_old_ref = False
                            else:
                                pas_no += 1        
                                new_isoform = gene+'|'+str(sjc_num)+'.'+str(sjc_no)+'.'+str(pas_num)+'.'+str(pas_no)+'|'+id
                                import copy
                                info_copy = copy.deepcopy(info)
                                info_copy.append(str(pas))
                                info_copy[3] = new_isoform
                                line = '\t'.join(info_copy)+'\n'
                                file.write(line)
                                
                        if no_old_ref:
                            new_isoform = gene+'|'+str(sjc_num)+'.'+str(sjc_no)+'.'+str(pas_num)+'.'+'0'+'|'+id
                            import copy
                            info_copy = copy.deepcopy(info)
                            info_copy[3] = new_isoform
                            line = '\t'.join(info_copy+['#'])+'\n'
                            file.write(line)
    time = format_datetime()
    print(f'[{time}] Correcting bed file...')

    with open(args.output_path+'/added_pas.bed', 'r') as file,\
        open(args.output_path+'/pas_refined.bed', 'w') as out:
        for line in file:
            info = line.strip().split('\t')
            start,end,strand,blocks,starts,pas = int(info[1]),int(info[2]),info[5],info[10].split(',')[:-1],info[11].split(',')[:-1],info[12]
            if info[-1] != '#':
                if strand == '+':
                    new_end = int(pas)
                    info[2] = str(new_end)
                    blocks[-1] = str(new_end - end + int(blocks[-1]))
                if strand == '-':
                    new_start = int(pas) - 1
                    info[1] = str(new_start)
                    blocks[0] = str(start - new_start + int(blocks[0]))
                    starts = ['0'] + [str(start - new_start + int(i)) for i in starts[1:]]
            info[10] = ','.join(blocks)
            info[11] = ','.join(starts)
            outline = '\t'.join(info[:-1])
            out.write(outline+'\n')
    os.remove(args.output_path+'/added_pas.bed')
    print(f'[{time}] Done!')

if __name__ == '__main__':
	main()