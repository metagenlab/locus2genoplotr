
def main():
    from locus2genoplotr import locus2genoplotr
    import argparse
    import sys
    import re
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-l",'--locus', type=str, help="locus_tag")
    parser.add_argument("-r",'--reference', type=str, help="reference genbank")
    parser.add_argument("-q",'--query', type=str, help="target genbank(s)", nargs='+')
    parser.add_argument("-ls",'--left_side_window', type=int, help="left siden window")
    parser.add_argument("-rs",'--right_side_window', type=int, help="right side window")
    parser.add_argument("-i",'--min_identity', type=int, help="minimum identity for blast", default=50)
    parser.add_argument("-sl",'--show_labels', action="store_false", help="do not show show labels")
    parser.add_argument("-s", '--samtools_depth', default=False, help="add depth plot from samtool depth (only for the last query). Should match the chromosome/contig names of the gbk.")
    parser.add_argument("-x",'--tblastx', action="store_true", help="execute tblastx and not blasn (6 frame translations)")
    parser.add_argument("-g",'--gc_plot', action="store_true", help="Show GC plot")


    args = parser.parse_args()

    L = Locus2genoplotR(args.locus,
                        args.reference,
                        args.query,
                        left_side=args.left_side_window,
                        right_side=args.right_side_window,
                        tblastx=args.tblastx)

    if args.query:
        start, end, flip_record = L.blast_target_genbank()

        all_records = [L.ref_sub_record] + L.sub_record_list
        names = [record.description for record in all_records]


        for i, name in enumerate(names):
            tmp_name = re.sub(', complete sequence.','', name)
            tmp_name = re.sub('strain ','', tmp_name)
            tmp_name = re.sub('Klebsiella pneumoniae ','K.p ', tmp_name)
            tmp_name = re.sub(', complete genome.','', tmp_name)
            tmp_name = re.sub('-contig_48','', tmp_name)
            tmp_name = re.sub('subsp. pneumoniae','', tmp_name)
            names[i] = tmp_name
        blast_result_files = L.record_list2blast(all_records, args.min_identity)
        gbk_list = L.write_genbank_subrecords(all_records)

        L.record2multi_plot(gbk_list,
                            blast_result_files,
                            names,
                            args.show_labels,
                            last_record_start=start,
                            last_record_end=end,
                            flipped_record=flip_record,
                            depth_file=args.samtools_depth)
    else:
        gbk_list = L.write_genbank_subrecords([L.ref_sub_record])
        L.record2single_plot(gbk_list[0], 'test', show_labels=args.show_labels, target_locus=args.locus)

