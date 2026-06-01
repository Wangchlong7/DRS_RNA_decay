import sys
import argparse
import os
from datetime import datetime

def format_datetime():
    now = datetime.now()
    formatted_time = now.strftime('%Y-%m-%d %H:%M:%S')
    return formatted_time

def end_unmatch_dist(blocks1,blocks2,orientation):
    blocks = [blocks1,blocks2]
    block_num = [len(blocks1),len(blocks2)]
    if orientation == 'left':
        end= [[int(block.split('-')[0]) for block in blocks1],[int(block.split('-')[0]) for block in blocks2]]
        limit = [end[0][0],end[1][0]]
        if limit[0] > limit[1]:
            short,long = 0,1
        elif limit[0] < limit[1]:
            short,long = 1,0
        else:
            return 0
        seed = limit[short]
        forsearch = end[long]
        dist = [(seed-end)for end in forsearch]
        dist_correct = [float('inf') if i < 0 else i for i in dist]
        min_dist = min(map(abs,dist_correct))
        index = [i for i, x in enumerate(dist_correct) if abs(x) == min_dist]
        if index[0] >= 1:
            rest_blocks = [blocks[long][i] for i in range(index[0])]
            rest_len = sum([abs(int(block.split("-")[1]) - int(block.split("-")[0]) + 1) for block in rest_blocks])
        else:
            rest_len = 0
        if long == 0:
            unmatch_len = rest_len + dist[index[0]]
        else:
            unmatch_len = -(rest_len + dist[index[0]])

    if orientation == 'right':
        end= [[int(block.split('-')[1]) for block in blocks1],[int(block.split('-')[1]) for block in blocks2]]
        limit = [end[0][-1],end[1][-1]]
        if limit[0] > limit[1]:
            short,long = 1,0
        elif limit[0] < limit[1]:
            short,long = 0,1
        else:
            return 0
        seed = limit[short]
        forsearch = end[long]
        dist = [(seed-end )for end in forsearch]
        dist_correct = [float('inf') if i > 0 else i for i in dist]
        min_dist = min(map(abs,dist_correct))
        index = [i for i, x in enumerate(dist_correct) if abs(x) == min_dist]
        if index[0] <= block_num[long]-2:
            rest_blocks = [blocks[long][i] for i in range(index[0]+1,block_num[long])]
            rest_len = sum([abs(int(block.split("-")[1]) - int(block.split("-")[0]) + 1) for block in rest_blocks])
        else:
            rest_len = 0
        if long == 0:
            unmatch_len = rest_len + abs(dist[index[0]])
        else:
            unmatch_len = -(rest_len + abs(dist[index[0]]))   
    return unmatch_len

def create_dict(trmap):
     with open(trmap, 'r') as input:
        read_dict={}
        for line in input:
            if ">" in line.strip():
                info = line.strip().split(" ")
                read = info[0][1:]
                if not read in read_dict.keys():
                    read_dict[read] = {}
                    read_dict[read]['chrom'] = info[1].split(':')[0]
                    read_dict[read]['strand'] = info[2]
                    read_dict[read]['read_block'] = info[3]
                    read_dict[read]['ref_block'] = {}
            else:
                info = line.strip().split("\t")
                if info[0] in ['=','c','k']:
                    isoform,ref_block = info[5],info[6]
                    read_dict[read]['ref_block'][ref_block] = (info[0],isoform)
        return read_dict

def assignment(read_dict,read):
    code_list = []
    for ref in read_dict['ref_block']:
        code_list.append(read_dict['ref_block'][ref][0])
    read_block_list = read_dict['read_block'].split(',')
    ref_list = []
    end_unmatch_list = []
    for ref_block in read_dict['ref_block']:
        ref_block_list = ref_block.split(',')
        ref_list.append(ref_block)
        end_unmatch_list.append([end_unmatch_dist(read_block_list,ref_block_list,'left'),end_unmatch_dist(read_block_list,ref_block_list,'right')])
    
    end_unmatch_size = [abs(i[0]) + abs(i[1]) for i in end_unmatch_list]
    
    if '=' not in code_list:
        min_value = min(end_unmatch_size)
        min_index = [i for i, x in enumerate(end_unmatch_size) if x == min_value][0]
    else:
        min_index_list = [i for i, x in enumerate(code_list) if x == '=']
        strand = read_dict['strand']
        if strand == '+':
            end3_list = [abs(end_unmatch_list[i][1]) for i in min_index_list]
        else:
            end3_list = [abs(end_unmatch_list[i][0]) for i in min_index_list]

        min_end3_diff = min(end3_list)
        min_index = [min_index_list[i] for i, x in enumerate(end3_list) if x == min_end3_diff][0]
        min_value = end_unmatch_size[min_index]

    
    strand = read_dict['strand'] 
    if strand == '+':
        unmatch_info_end5 = (str(end_unmatch_list[min_index][0]))
        unmatch_info_end3 = (str(end_unmatch_list[min_index][1]))
    else:
        unmatch_info_end5 = (str(end_unmatch_list[min_index][1]))
        unmatch_info_end3 = (str(end_unmatch_list[min_index][0]))
    unmatch_info = '|'.join(['tss:' + unmatch_info_end5,'tes:' + unmatch_info_end3])

    transcript_id = read_dict['ref_block'][ref_list[min_index]][1]
    ref_block = ref_list[min_index].split(',')
    ref_length = sum([abs(int(block.split("-")[1]) - int(block.split("-")[0]) + 1) for block in ref_block])

    unmatch_info_dist = str(min_value/ref_length)
    
    gene_id = transcript_id.split("|")[0]
    selected_code = read_dict['ref_block'][ref_list[min_index]][0]
    chrom = read_dict['chrom']
    assignment_record = '\t'.join([read,chrom,strand,selected_code,transcript_id,gene_id,unmatch_info_dist,unmatch_info,read_dict['read_block'],ref_list[min_index],str(ref_length)])
    return assignment_record

def extract_fl_bed(read_assignment,corrected_reads,fl_read_assignment,fl_corrected_reads):
    with open(read_assignment,'r') as input, open(fl_read_assignment,'w') as output:
        fl_read_dict = {}
        for line in input:
            info = line.strip().split('\t')
            code = info[3]
            if not code == 'ambiguous':
                if code == '=':
                    if float(info[6]) < 0.4:
                        output.write(line)
                        if info[0] not in fl_read_dict:
                            fl_read_dict[info[0]] = info[4]
                elif len(info[-3].split(',')) == 1 and float(info[6]) < 0.2:
                    output.write(line)
                    if info[0] not in fl_read_dict:
                            fl_read_dict[info[0]] = info[4]
        with open(corrected_reads) as file,\
            open(fl_corrected_reads,'w') as out:
            for line in file:
                read= line.strip().split('\t')[3]
                if read in fl_read_dict:
                    info = line.strip().split('\t')
                    info[4] = fl_read_dict[read]
                    out.write('\t'.join(info)+'\n')

def extract_fl_bed_v2(read_assignment,corrected_reads,fl_read_assignment,fl_corrected_reads):
    with open(read_assignment,'r') as input, open(fl_read_assignment,'w') as output:
        fl_read_dict = {}
        for line in input:
            info = line.strip().split('\t')
            code = info[3]
            if not code == 'ambiguous':
                if code == '=':
                    strand = info[2]
                    blocks = info[-1]
                    if strand == '+':
                        exon1st = abs(int(blocks.split(',')[0].split('-')[0]) - int(blocks.split(',')[0].split('-')[1])) + 1
                    else:
                        exon1st = abs(int(blocks.split(',')[-1].split('-')[0]) - int(blocks.split(',')[-1].split('-')[1])) + 1
                    if exon1st >= 25:
                        output.write(line)
                        if info[0] not in fl_read_dict:
                            fl_read_dict[info[0]] = info[4]
                elif len(info[-3].split(',')) == 1 and float(info[6]) < 0.2:
                    output.write(line)
                    if info[0] not in fl_read_dict:
                            fl_read_dict[info[0]] = info[4]
        with open(corrected_reads) as file,\
            open(fl_corrected_reads,'w') as out:
            for line in file:
                read= line.strip().split('\t')[3]
                if read in fl_read_dict:
                    info = line.strip().split('\t')
                    info[4] = fl_read_dict[read]
                    out.write('\t'.join(info)+'\n')

def quantify(input):
    import pandas as pd
    full_length,human_isoform,dro_isoform,human_gene,dro_gene = input
    with open(full_length) as file,open(human_isoform,'w') as hi_out,open(dro_isoform,'w') as di_out,\
        open(human_gene,'w') as hg_out,open(dro_gene,'w') as dg_out:
        hi_dict = {}
        di_dict = {}
        hg_dict = {}
        dg_dict = {}
        for line in file:
            chrom,isoform = line.strip().split('\t')[1],line.strip().split('\t')[4]
            gene = isoform.split('|')[0]
            if chrom in ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']:
                if isoform not in hi_dict:
                    hi_dict[isoform] = 0
                hi_dict[isoform] += 1
                if gene not in hg_dict:
                    hg_dict[gene] = 0
                hg_dict[gene] += 1
            if chrom in ['Dro_2L','Dro_2R','Dro_3L','Dro_3R','Dro_4','Dro_mitochondrion_genome','Dro_rDNA','Dro_X']:
                if isoform not in di_dict:
                    di_dict[isoform] = 0
                di_dict[isoform] += 1
                if gene not in dg_dict:
                    dg_dict[gene] = 0
                dg_dict[gene] += 1
        hi_out.write('\t'.join(['isoform','count'])+'\n')
        di_out.write('\t'.join(['isoform','count'])+'\n')
        hg_out.write('\t'.join(['gene','count'])+'\n')
        dg_out.write('\t'.join(['gene','count'])+'\n')
        for isoform in hi_dict:
            hi_out.write('\t'.join([isoform,str(hi_dict[isoform])])+'\n')
        for isoform in di_dict:
            di_out.write('\t'.join([isoform,str(di_dict[isoform])])+'\n')
        for gene in hg_dict:
            hg_out.write('\t'.join([gene,str(hg_dict[gene])])+'\n')
        for gene in dg_dict:
            dg_out.write('\t'.join([gene,str(dg_dict[gene])])+'\n')

def main():
    time = format_datetime()
    print("---------------------------------------------------------------------------------------------------------------------------------")
    print(f'[{time}] python' ,' '.join(sys.argv), sep = '\t')

    parser = argparse.ArgumentParser(description='Exacting full length read-isoform pairs from trmap results')
    parser.add_argument('-1','--trmap',required = True, help="the path to load trmap txt results")
    parser.add_argument('-2','--corrected_reads_bed',required = True, help="the path to load corrected_reads.bed")
    parser.add_argument('-3','--prefix',required = True, help="the prefix to name output files")
    parser.add_argument('-o','--output_path',required = True, help="the path to output assignment txt results")
    
    args = parser.parse_args(sys.argv[1:])
    time = format_datetime()
    if os.path.exists(args.output_path+'/'+args.prefix+'.read_assignment.txt'):
        print(args.prefix+'.read_assignment.txt exist')
    else:
        print(f'[{time}] Extracting information from trmap results...')
        
        read_dict = create_dict(args.trmap)

        time = format_datetime()
        print(f'[{time}] Extract information from trmap results: Done')

        process_num = 0
        
        with open(args.output_path+'/'+args.prefix+'.read_assignment.txt','w') as output:      
            for read in read_dict:
                if len(read_dict[read]['ref_block']) > 0:
                    out_record = assignment(read_dict[read],read)
                    output.write(out_record + '\n')
                    process_num += 1
                if (process_num) % 5000 == 0:
                    time = format_datetime()
                    print(f'[{time}] Processed {process_num} lines')

        time = format_datetime()
        print(f'[{time}] {process_num} reads have been successfully processed')
    if os.path.exists(args.output_path+'/'+args.prefix+'.fl_read_assignment.txt'):
        print(args.prefix+'.fl_read_assignment.txt exist')
        print(args.prefix+'.fl_corrected_reads.txt exist')
    else:
        extract_fl_bed(args.output_path+'/'+args.prefix+'.read_assignment.txt',
                    args.corrected_reads_bed,
                    args.output_path+'/'+args.prefix+'.fl_read_assignment.txt',
                    args.output_path+'/'+args.prefix+'.fl_corrected_reads.bed')
    quantify([args.output_path+'/'+args.prefix+'.fl_read_assignment.txt',
              args.output_path+'/'+args.prefix+'.human_isoform_count.txt',
              args.output_path+'/'+args.prefix+'.dro_isoform_count.txt',
              args.output_path+'/'+args.prefix+'.human_gene_count.txt',
              args.output_path+'/'+args.prefix+'.dro_gene_count.txt'])
    time = format_datetime()
    print(f'[{time}] full length reads have been extracted')
    
if __name__ == '__main__':
	main()
