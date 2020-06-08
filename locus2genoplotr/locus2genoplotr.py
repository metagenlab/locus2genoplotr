#!/usr/bin/env python

import numpy as np
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri as numpy2ri
from rpy2.robjects.packages import importr
import rpy2.rlike.container as rlc
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()
#import pandas.rpy.common as com


class Locus2genoplotR():

    def __init__(self, query_locus,
                 reference_genbank,
                 target_genbank_list,
                 left_side=0,
                 right_side=0,
                 tblastx=False):

        from Bio import SeqRecord, SeqIO
        self.query_locus = query_locus
        self.target_genbank_list = target_genbank_list
        self.right_side = right_side
        self.left_side = left_side
        self.tblastx = tblastx

        if type(reference_genbank) == str:
            self.reference_record = [i for i in SeqIO.parse(open(reference_genbank), 'genbank')]
        elif isinstance(reference_genbank, SeqRecord.SeqRecord):
            self.reference_record = [reference_genbank]
        elif type(reference_genbank) == list and isinstance(reference_genbank[0], SeqRecord.SeqRecord):
            self.reference_record = reference_genbank
        else:
            raise IOError('wrong input reference')

        self.ref_locus_seqrecord, self.ref_sub_record, self.ref_feature, self.region_start, self.region_end = self.get_target_locus_region(self.query_locus,
                                                                                     self.reference_record,
                                                                                     right_side=right_side,
                                                                                     left_side=left_side,
                                                                                     flip_record=False)

        '''
        self.orf_leading_strand_count = 0
        self.orf_lagging_strand_count = 0
        for feature in self.ref_sub_record.features:
            if feature.type == 'CDS':
                if feature.location.strand == 1:
                    self.orf_leading_strand_count += 1
                else:
                    self.orf_lagging_strand_count += 1

        print 'orf count:', self.orf_leading_strand_count, self.orf_lagging_strand_count
        '''

    def get_target_locus_region(self, target_locus_tag,
                                records,
                                right_side=0,
                                left_side=0,
                                flip_record=False):

        from Bio.SeqRecord import SeqRecord
        match=False
        for record in records:
            #if record.name in ['NZ_JUBO01000032', 'NZ_JUBO01000028']:
            #    flip_record=True
            #    print 'match!-----------------------------------------------------'
            if flip_record:
                #print ('flipping record', record.id)
                name = record.description
                tmp_name = re.sub(', whole genome shotgun sequence.','', name)
                tmp_name = re.sub('strain ','', tmp_name)
                #tmp_name = re.sub('Klebsiella pneumoniae ','K.p ', tmp_name)
                tmp_name = re.sub(', complete genome.','', tmp_name)
                record.description = tmp_name
                record = record.reverse_complement(id=record.id,
                                                   name=record.id,
                                                   description=tmp_name)

            for feature in record.features:
                if 'locus_tag' in feature.qualifiers:
                    if feature.qualifiers['locus_tag'][0] == target_locus_tag:




                        query_start = int(feature.location.start)
                        query_end = int(feature.location.end)


                        region_start = query_start-left_side
                        region_end = query_end+right_side


                        if region_start<0:
                            region_start=0

                        if region_end > len(record.seq):
                            region_end=len(record.seq)
                        print ('extraction from %s to %s' % (region_start, region_end))

                        query_sub_record = record[region_start:region_end]

                        target_record = SeqRecord(feature.extract(record.seq).translate(),
                                                               id=self.query_locus,
                                                               name=self.query_locus,
                                                               description="target locus")
                        match=True
                        match_feature = feature
                        break
                if 'protein_id' in feature.qualifiers:
                    if feature.qualifiers['protein_id'][0] == target_locus_tag:



                        query_start = int(feature.location.start)
                        query_end = int(feature.location.end)
                        region_start = query_start-left_side
                        if region_start<0:
                            region_start=0
                        region_end = query_end+right_side
                        if region_end > len(record.seq):
                            region_end=len(record.seq)
                        print ('extraction from %s to %s' % (region_start, region_end))



                        query_sub_record = record[region_start:region_end]

                        target_record = SeqRecord(feature.extract(record.seq).translate(),
                                                               id=self.query_locus,
                                                               name=self.query_locus,
                                                               description="target locus")
                        match=True
                        match_feature = feature
                        break
        print('match,', match)
        if not match:
            return False
        else:
            return target_record, query_sub_record, match_feature, region_start, region_end

    def blast_target_genbank(self):
        # write target locus to fasta
        # for each target genbank
        #    gbk2faa
        #    blast
        #    read blast => keep best hit in memory
        from tempfile import NamedTemporaryFile
        from Bio import SeqRecord, SeqIO
        from io import StringIO
        import gbk2faa
        from TPutils import blast_utils

        self.record_id2best_blast_locus_tag = {}
        self.accession2record = {}
        self.record_order = []

        temp_ref = NamedTemporaryFile(delete=False, mode='w')
        fastastr = StringIO()

        SeqIO.write(self.ref_locus_seqrecord, fastastr, 'fasta')

        temp_ref.write(fastastr.getvalue())
        temp_ref.flush()

        self.best_hit_list = []

        self.sub_record_list = []

        for target_record in self.target_genbank_list:

            if type(target_record) == str:
                rec_list = [i for i in SeqIO.parse(open(target_record), 'genbank')]
                for rec in rec_list:
                    self.accession2record[rec.name] = rec
                    self.record_order.append(rec.name)

            elif isinstance(target_record, SeqRecord):
                self.accession2record[target_record.name] = target_record
                self.record_order.append(target_record.name)
                rec_list = target_record
            elif type(target_record) == list and isinstance(target_record[0], SeqRecord):
               for rec in target_record:
                   self.accession2record[rec.name] = rec
                   self.record_order.append(rec.name)
               rec_list = target_record

            else:
                raise IOError('wrong target gbk format')

            temp_target = NamedTemporaryFile(delete=False, mode='w')
            fastastr = StringIO()

            gbk2faa.gbk2faa(rec_list, lformat=True, output_handle=fastastr)
            #SeqIO.write(rec_list, fastastr, 'fasta')
            temp_target.write(fastastr.getvalue())
            temp_target.flush()

            B = blast_utils.Blast(temp_ref.name, temp_target.name, protein=True)
            print ('formatting...')
            B.format_database()
            print ('blasting...')
            B.run_blastp()
            best_hit = B.best_hit_list[0]
            print ('best hit locus:', best_hit[1])

            target_seq, target_sub_record, feature, start, end = self.get_target_locus_region(best_hit[1],
                                                                         rec_list,
                                                                         right_side=self.right_side,
                                                                         left_side=self.left_side,
                                                                         flip_record=False)
            flip_record = False

            '''
            orf_leading_strand_count = 0
            orf_lagging_strand_count = 0
            for feature in target_sub_record.features:
                if feature.type == 'CDS':

                    if feature.location.strand == 1:
                        orf_leading_strand_count += 1
                    else:
                        orf_lagging_strand_count += 1

            print 'orf count comp', self.orf_leading_strand_count > self.orf_lagging_strand_count, orf_leading_strand_count < orf_lagging_strand_count
            print self.orf_leading_strand_count, self.orf_lagging_strand_count, orf_leading_strand_count, orf_lagging_strand_count
            if (self.orf_leading_strand_count > orf_leading_strand_count) and (self.orf_lagging_strand_count < orf_lagging_strand_count):
            '''
            if (feature.location.strand != self.ref_feature.location.strand):
                print ('first trial')
                print (feature.location.strand != self.ref_feature.location.strand)
                target_seq, target_sub_record, feature, start, end = self.get_target_locus_region(best_hit[1],
                                                                             rec_list,
                                                                             right_side=self.right_side,
                                                                             left_side=self.left_side,
                                                                             flip_record=True)
                flip_record = True

            print (target_sub_record.id, '#############################################################################')
            '''
            if target_sub_record.id =="NC_009648.1":
                print 'fpli!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print feature.location.strand != self.ref_feature.location.strand
                target_seq, target_sub_record, feature = self.get_target_locus_region(best_hit[1],
                                                                             rec_list,
                                                                             right_side=self.right_side,
                                                                             left_side=self.left_side,
                                                                             flip_record=True)
            '''
            if target_sub_record:
                self.sub_record_list.append(target_sub_record)

        return start, end, flip_record
        #print 'sub record list', self.sub_record_list


    def write_genbank_subrecords(self, subrecord_list):
        from tempfile import NamedTemporaryFile
        from io import StringIO

        out_names = []

        for record in subrecord_list:

            temp_ref = NamedTemporaryFile(delete=False, mode='w')
            fastastr = StringIO()

            SeqIO.write(record, fastastr, 'genbank')

            temp_ref.write(fastastr.getvalue())
            temp_ref.flush()
            out_names.append(temp_ref.name)
        return out_names

    def record_list2blast(self, record_list, min_identity):
        '''

        :param record_list:
        :return:
        '''

        blast_result_files = []

        from TPutils import blast_utils

        for i in range(0, len(record_list)-1):
            # write fasta nucl to temp files
            # blast
            # keep output name into memory
            #print 'record 1\n--------------------', record_list[i]
            #print 'record 2\n--------------------', record_list[i+1]
            B = blast_utils.Blast(record_list[i], record_list[i+1], protein=False)
            print('formatting...')
            B.format_database()
            print('blasting...')
            if not self.tblastx:
                blast_file = B.run_blastn(min_identity=min_identity)
            else:
                blast_file = B.run_tblastx(evalue=0.0005)

            blast_result_files.append(blast_file)

        return blast_result_files

    def record_list2promer(self):
        pass

    def record2single_plot(self, 
                           genbank, 
                           name, 
                           show_labels=True, 
                           target_locus=False,
                           output_name="single_plot"):

        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri as numpy2ri
        numpy2ri.activate()
        from Bio.Graphics import GenomeDiagram
        from Bio import SeqIO

        record = SeqIO.read(genbank, "genbank")

        #robjects.r.assign('df3', com.convert_to_r_dataframe(table_comparison))

        gd_diagram = GenomeDiagram.Diagram("test",y=0.15)
        gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
        gd_feature_set = gd_track_for_features.new_set()
        for feature in record.features:
            if feature.type == "CDS":
                if len(gd_feature_set) % 2 == 0:
                    color = "blue"
                else:
                        color = "blue"
                    #print target_locus
                if target_locus:
                    if feature.qualifiers['locus_tag'][0] == target_locus:
                        color = "red"
            elif feature.type =="tRNA":
                color = "orange"
            elif feature.type =="rRNA":
                color = "pink"
            else:
                print(feature.type)
                continue
            gd_feature_set.add_feature(feature, sigil="ARROW", color=color, label=True, label_position="middle",label_strand=1, label_size=12, label_angle=25,arrowshaft_height=0.5)

        hauteur = 200
        largeur = len(record)/25
        gd_diagram.draw(format="linear", 
                            orientation="landscape", pagesize=(hauteur,largeur),
                    fragments=1, start=0, end=len(record))
        gd_diagram.write(f"{output_name}.svg", "SVG")


    def record2multi_plot(self,
                          genbank_list,
                          blast_file_list,
                          names,
                          show_labels=True,
                          depth_file=False,
                          last_record_start=False,
                          last_record_end=False,
                          flipped_record=False,
                          show_GC_last=False):

        import rpy2.robjects as robjects
        import rpy2.robjects.numpy2ri as numpy2ri
        from Bio import SeqIO
        numpy2ri.activate()

        #robjects.r.assign('df3', com.convert_to_r_dataframe(table_comparison))

        robjects.r('''
        genbank_list <- list()
        blast_list <- list()
        ''')
        robjects.r.assign('genbank_list', genbank_list)
        robjects.r.assign('blast_list', blast_file_list)
        robjects.r.assign('names', names)

        last_seq = SeqIO.read(genbank_list[-1], 'genbank').seq
        robjects.r.assign('last_seq', str(last_seq))
        robjects.r.assign('show_GC_last', show_GC_last)

        if depth_file:
            robjects.r.assign('last_record_start', last_record_start)
            robjects.r.assign('last_record_end', last_record_end)
            robjects.r.assign('last_record_flipped', flipped_record)
            robjects.r.assign('depth_file', depth_file)
        else:
            robjects.r.assign('depth_file', False)

        if show_labels:
            plot = '''
                p <- plot_gene_map(dna_segs=dna_seg_list,
                               comparisons=comp_list,
                               main="",
                               dna_seg_scale=TRUE,
                               scale=TRUE,
                               global_color_scheme=c("per_id", "increasing", "grey", 1),
                               annotations=annots, 
                               annotation_height=0.4, 
                               annotation_cex=0.5,
                               override_color_schemes=TRUE,
                               plot_new=FALSE)
            '''
        else:
            plot = '''
            plot_gene_map(dna_segs=dna_seg_list,
                          comparisons=comp_list,
                          main="",
                          dna_seg_scale=TRUE,
                          scale=TRUE,
                          override_color_schemes=TRUE,
                          global_color_scheme=c("per_id", "increasing", "grey", 1),
                          plot_new=FALSE)
            '''

        robjects.r('''
            library(genoPlotR, quietly = TRUE)
            library(Cairo, quietly = TRUE)
            library(squash, quietly = TRUE)
            library(gridBase, quietly = TRUE)
            library(seqinr, quietly = TRUE)
            library(genbankr, quietly = TRUE)


            slidingwindowplot <- function(windowsize, inputseq)
            {
               starts <- seq(1, length(inputseq)-windowsize, by = windowsize)
               n <- length(starts)    # Find the length of the vector "starts"
               chunkGCs <- numeric(n) # Make a vector of the same length as vector "starts", but just containing zeroes
               for (i in 1:n) {
                    chunk <- inputseq[starts[i]:(starts[i]+windowsize-1)]
                    chunkGC <- GC(chunk)
                    chunkGCs[i] <- chunkGC
               }
               plot(starts,chunkGCs,type="l", xaxt='n', ann=F, xaxs="i", yaxs="i")
               abline(h=GC(inputseq))
            }
            
            
            if (depth_file != FALSE) {
            

                
                depth_data <- read.table(depth_file, header=F)$V3
                
                if (last_record_flipped == TRUE){
                print('YES')
                    depth_data <- rev(depth_data)
                }else{
                print('NO')
                }
                depth_subset <- depth_data[as.numeric(last_record_start):as.numeric(last_record_end)]
            }else{
             print('NODEPTH')

            }
            
            if (depth_file != FALSE || show_GC_last != FALSE){
            print ('################# add viewport!!!!!!!!!!!!!!!!!')
                # for viewport and pdf dimension
                second_plot_height <- 0.6 
                pdf_height_scale <- 2.5#2.5            
            }else{
                 second_plot_height <- 0.1  
                 pdf_height_scale <- 2 
            }

            dna_seg_list <- list()
            for (i in 1:length(genbank_list)){
                #print(genbank_list[[i]])
                dna_seg_list[[i]] <- try(read_dna_seg_from_file(genbank_list[[i]], tagsToParse=c("CDS", "assembly_gap", "rRNA", "tRNA")))
                




                w <- which(dna_seg_list[[i]]$feature=='assembly_gap')
                dna_seg_list[[i]]$gene_type <- rep("arrows", length(dna_seg_list[[i]]$gene_type))
                
                dna_seg_list[[i]]$col <- rep("cadetblue1", length(dna_seg_list[[i]]$col))
		        dna_seg_list[[i]]$fill <- rep("blue", length(dna_seg_list[[i]]$fill))
                dna_seg_list[[i]]$lty <- rep(0, length(dna_seg_list[[i]]$lty))
                dna_seg_list[[i]]$lwd <- rep(2, length(dna_seg_list[[i]]$lwd))
                #dna_seg_list[[i]]$pch <- rep(3, length(dna_seg_list[[i]]$pch))
                #dna_seg_list[[i]]$cex <- rep(0.5, length(dna_seg_list[[i]]$cex))
                dna_seg_list[[i]]$gene_type[w] <- "bars"
                dna_seg_list[[i]]$col[w] <- "red"
                dna_seg_list[[i]]$lty[w] <- 1
                dna_seg_list[[i]]$proteinid[w] <- ""
                dna_seg_list[[i]]$synonym[w] <- ""
                dna_seg_list[[i]]$name[w] <- ""
                w <- which(abs(dna_seg_list[[i]]$length)<20000)
                dna_seg_list[[i]] <- dna_seg_list[[i]][w,]

                w <- which(dna_seg_list[[i]]$feature=='tRNA')
                dna_seg_list[[i]]$gene_type[w] <- "bars"
                dna_seg_list[[i]]$col[w] <- "orange"
                dna_seg_list[[i]]$lty[w] <- 1
		        #print(head(dna_seg_list[[i]]))	


                w <- which(dna_seg_list[[i]]$gene=='-')
                dna_seg_list[[i]]$name[w] <- ""

                #print(dna_seg_list[[i]])

            }

            names(dna_seg_list) <- names


            all_id <- c()

            comp_list <- list()
            for (i in 1:length(blast_list)){

                comp_list[[i]] <- try(read_comparison_from_blast(blast_list[[i]]))

                fake_row1 <- comp_list[[i]][1,]
                fake_row2 <- comp_list[[i]][1,]

                fake_row1$per_id <- 20
                fake_row2$per_id <- 100
                fake_row1$aln_len <- 1
                fake_row2$aln_len <- 1
                fake_row1$mism <- 0
                fake_row2$mism <- 0
                fake_row1$gaps <- 0
                fake_row2$gaps <- 0

                fake_row1$start1 <- 1
                fake_row2$start1 <- 1
                fake_row1$start2 <- 1
                fake_row2$start2 <- 1
                fake_row1$end1 <- 2
                fake_row2$end1 <- 2
                fake_row1$end2 <- 2
                fake_row2$end2 <- 2

                comp_list[[i]] <- rbind(comp_list[[i]], fake_row1, fake_row1)

                comp_list[[i]]$col <- apply_color_scheme(comp_list[[i]]$per_id, "grey")

                all_id <- c(all_id,comp_list[[i]]$per_id)
            }

            rng <- c(0.20,1)#range(all_id)

            #if (diff(rng) == 0){
            #  col <- rep(grey(0.5), length(all_id))
            #} else {
            #  level <- 0.75-((all_id-rng[1])/diff(rng))*0.5
            #  col <- grey(level)
            #}

            #col <- grey(sort(unique(level)))

            breaks <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
            col <- grey(breaks)
            transparency <- 0.9
            tpc <- floor(transparency*256)
            tpc <- sprintf("%%X", tpc)
            if (nchar(tpc) == 1) tpc <- paste("0", tpc, sep="")
            col <- paste(col, tpc, sep="")
            # color map
            map <- makecmap(sort(unique(rng)), n=10)
            map$colors<-col
            map$include.lowest<-TRUE
            map$breaks <- rep("",length(map$colors)+1)
            map$breaks[1] <-rng[2]
            map$breaks[length(map$colors)+1]<-rng[1]

            annots <- lapply(dna_seg_list, function(x){
               mid <- middle(x)
               annot <- annotation(x1=mid, text=x$name, rot=30)
               idx <- grep("^[^B]", annot$text, perl=TRUE)
               annot #[idx[idx %% 4 == 0],]
             })

            height <- length(blast_list)*pdf_height_scale
            len_seq <- (dna_seg_list[[1]]$end[length(dna_seg_list[[1]]$end)] - dna_seg_list[[1]]$start[1])
            len_kb <- len_seq/1000
            width <- len_kb/2
            #print("dimensions:")
            #print (height)
            #print (width)
            
            CairoPDF('test2.pdf',height=height,width=width)# 4,14 / 3.8 (yersinia)/2 (oxa)

                xlims <- list(c(1,50000), c(1,50000))
                plot.new()

                pushViewport(viewport(layout=grid.layout(2,1,
                                heights=unit(c(1, second_plot_height), rep("null", 2))), # 0.1
                                name="overall_vp"))
                    pushViewport(viewport(layout.pos.row=1, name="panelA"))
                    %s
                    upViewport()
                    pushViewport(viewport(layout.pos.row=2, name="panelB"))
                    hkey(map)
                    upViewport()
                    
                    if (depth_file != FALSE) {
                        pushViewport(viewport(layout.pos.row=2, name="panelC"))
                        par(new=TRUE, fig=gridFIG(), las=1,mar=c(5, (length(unlist(strsplit("K. penumoniae 1", "")))/2)+0.3, 0, 1.7))
                        plot(seq(1:length(depth_subset)),depth_subset,type='l',col='light grey', xaxt='n', ann=F, xaxs="i", yaxs="i", ylim=c(0, 3*median(depth_data)))
                        
                        mtext("Depth", side=2, line=3)
                        
                        xx <- c(seq(1:length(depth_subset)), rev(seq(1:length(depth_subset))))
                        yy <- yy <- c(rep(0, length(depth_subset)), rev(depth_subset)) 
                        polygon(xx, yy, col='light grey')
                        abline(h=median(depth_data), col="red")
                        axis(1, at = seq(0, length(depth_subset), 10000), labels = seq(0, length(depth_subset), 10000)/1000, las=2)
                        upViewport()
                    }
                
                    if (show_GC_last != FALSE) {
                        pushViewport(viewport(layout.pos.row=2, name="panelC"))
                        par(new=TRUE, fig=gridFIG(), las=1,mar=c(5, (length(unlist(strsplit("K. penumoniae 2", "")))/2)+0.3, 0, 2.0)) 
                        myseq <- unlist(strsplit(last_seq, ""))
                        slidingwindowplot(500, myseq)
                        # ,xlab="Nucleotide start position",ylab="GC content"
                        axis(1, at = seq(0, length(myseq), 10000), labels = seq(0, length(myseq), 10000)/1000, las=2)
                        upViewport()                
                    }

            dev.off()
        ''' % plot)


if __name__ == '__main__':
    import argparse
    import sys
    import re
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-l",'--locus', type=str, help="locus_tag")
    parser.add_argument("-r",'--reference', type=str, help="reference genbank")
    parser.add_argument("-q",'--query', type=str, help="target genbank(s)", nargs='+')
    parser.add_argument("-o",'--output_name', help="output name", default='out')
    parser.add_argument("-ls",'--left_side_window', type=int, help="left siden window", default=15000)
    parser.add_argument("-rs",'--right_side_window', type=int, help="right side window", default=15000)
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
            tmp_name = re.sub(',.*','', name)
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
                            depth_file=args.samtools_depth,
                            show_GC_last=args.gc_plot)
    else:
        gbk_list = L.write_genbank_subrecords([L.ref_sub_record])
        L.record2single_plot(gbk_list[0], 
                             'test', 
                             show_labels=args.show_labels, 
                             target_locus=args.locus,
                             output_name=args.output_name)

