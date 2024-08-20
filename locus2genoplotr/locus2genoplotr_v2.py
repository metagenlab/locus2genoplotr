#!/usr/bin/env python

import subprocess
import tempfile
import os
from Bio import SeqIO
from Bio import SeqFeature, SeqIO
from Bio.SeqRecord import SeqRecord
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.DEBUG, encoding='utf-8')
from TPutils import blast_utils
from Bio.Seq import Seq
import glob
from itertools import chain
from pygenomeviz import Genbank, GenomeViz, load_dataset

class Subrecord():
    def __init__(self, gbk_record, out_prefix, basename=False):
        '''
        From gbk record
        - create a faa record
        - write faa to <faa_name>
        '''
                
        # write CDS to fasta
        if basename:
            self.basename = basename
        else:
            self.basename = None
        self.faa_file_path = out_prefix + '.faa'
        self.fna_file_path = out_prefix + '.fna'
        self.gbff_file_path = out_prefix + '.gbff'
        self.faa_record = self.write_coding_features_to_fasta(gbk_record, self.faa_file_path)
        SeqIO.write(gbk_record, self.fna_file_path, "fasta")
        SeqIO.write(gbk_record, self.gbff_file_path, "genbank")
        self.gbk_record = gbk_record
        
        
    def write_coding_features_to_fasta(self, seqrecord, faa_name):
        """
        Takes a list of SeqRecords as input, writes coding features as a temporary FASTA file with locus_tag
        (or protein_id if the locus_tag is missing) as header, and returns the path to that file.

        Args:
            genbank_file (str): Path to the input GenBank file.

        Returns:
            str: Path to the output FASTA file.
        """

        if isinstance(seqrecord, list):
            seqrecord_list = seqrecord
        elif isinstance(seqrecord, SeqRecord):
            seqrecord_list = [seqrecord]


        prot_record_list = []
        for record in seqrecord_list:
            # Write the coding features to the output FASTA
            for feature in record.features:
                if feature.type == "CDS" and "translation" in feature.qualifiers:
                    # Use locus_tag as header if available, otherwise use protein_id
                    if "locus_tag" in feature.qualifiers:
                        seqid = feature.qualifiers["locus_tag"][0]
                    elif "protein_id" in feature.qualifiers:
                        seqid = feature.qualifiers["protein_id"][0]
                        logging.warning(f"Missing locus_tag for protein_id {feature.qualifiers['protein_id'][0]}")
                    else:
                        logging.warning(f"Missing locus_tag and protein_id for feature at location {feature.location}")
                    if "product" in feature.qualifiers:
                        description = feature.qualifiers["product"][0]
                    else:
                        description = ""
                    # Write the sequence to the output file
                    seq = feature.qualifiers["translation"][0]
                    
                    prot_record = SeqRecord(
                        Seq(seq),
                        id=seqid,
                        name=seqid,
                        description=description,
                    )
                    prot_record_list.append(prot_record)
                    
                elif feature.type == "CDS" and "translation" not in feature.qualifiers:
                    logging.warning(f"Missing translation for feature at location {feature.location}")

        # Close the output file and return its path
        SeqIO.write(prot_record_list, faa_name, "fasta")

        return prot_record_list

        
class GenomeComp():
    def __init__(self, 
                 query_locus,
                 reference_gbk,
                 target_records_list,
                 upstream_bp=0,
                 downstream_bp=0,
                 tblastx=False,
                 output_name='out',
                 output_folder='l2p',
                 svg=False,
                 min_identity=50,
                 force_data_dir=False,
                 filp_records=True):

        from Bio import SeqRecord, SeqIO
        import shutil
         
        
        self.query_locus = query_locus
        self.upstream_bp = upstream_bp
        self.upstream_bp = upstream_bp
        self.tblastx = tblastx
        self.min_identity = min_identity
        self.cluster2highlights = {}
        self.basename2label = {}
        self.contig2basename = {}
        
        self.output_folder = output_folder
        
        if os.path.exists(self.output_folder):
            if not force_data_dir:
                raise IOError("Output folder already present in the working directory, remove it or use -f to force rewrite")
            else:
                shutil.rmtree(self.output_folder)
        
        os.mkdir(self.output_folder)
        
        self.reference_basename = os.path.basename(os.path.splitext(reference_gbk)[0])
        
        self.reference_records = { "basename": self.reference_basename, 
                                   "records": [i for i in SeqIO.parse(reference_gbk, "genbank")]}
        self.target_records_list = []
        for record in target_records_list:
            basename = os.path.basename(os.path.splitext(record)[0])
            records = [i for i in SeqIO.parse(record, "genbank")]
            self.target_records_list.append({"basename": basename, "records": records})
        
        
        if svg:
            self.output_format = 'svg'
            self.output_name = f'{output_name}.svg'
        else:
            self.output_format = 'pdf'
            self.output_name = f'{output_name}.pdf'

        logging.info("Extracting reference region")
        self.ref_subset = self.get_sequence_subset(self.reference_records["records"], 
                                                   query_locus, 
                                                   upstream_bp, 
                                                   downstream_bp, 
                                                   flip_record=False,
                                                   basename=self.reference_basename)
        
        self.reference_orientation = self.ref_subset["target_locus_strand"]
    
        self.ref_locus_record = SeqRecord.SeqRecord(self.ref_subset["seq"],
                                    id=self.query_locus,
                                    name=self.query_locus,
                                    description="reference locus")
        
        logging.info("Writing reference locus faa")
        self.ref_locus_faa = f'{self.output_folder}/reference.faa'
        SeqIO.write(self.ref_locus_record, self.ref_locus_faa, "fasta")
        
        logging.info("Writing reference subrecord faa & fna")
        self.ref_subrecord = Subrecord(self.ref_subset["subset_record"], 
                                       f'{self.output_folder}/{self.reference_records["basename"]}.sub',
                                       )
        
        #logging.info("Clustering reference proteins (highlight repeats)")
        #self.ref_record_clusters, self.locus2cluster = self.cluster_protein_sequences(self.ref_subrecord.faa_file_path)
               
        logging.info("Writing target faa(s)")
        self.target_full_records = [Subrecord(target["records"], 
                                    f'{self.output_folder}/{target["basename"]}.full',
                                    target["basename"]) for target in self.target_records_list]
        
        logging.info("Extracting bbh from targets")
        self.target_BBH_locus = [self.blast_best_hit(self.ref_locus_faa, i.faa_file_path) for i in self.target_full_records]
        logging.info("Extracting subrecords targets")
        self.target_subrecords = [self.get_sequence_subset(seqrecords["records"], 
                                                           target_id["subject_id"], 
                                                           upstream_bp, 
                                                           downstream_bp, 
                                                           basename=seqrecords["basename"],
                                                           flip_record=filp_records) for target_id, seqrecords in zip(self.target_BBH_locus, 
                                                                                                                             self.target_records_list)]

        self.target_faa_subrecords = [Subrecord(target_sub["subset_record"], 
                                               f'{self.output_folder}/{target_sub["basename"]}.sub',
                                               basename=target_sub["basename"]) for target_sub in self.target_subrecords]

        
        logging.info("Perform blast comparisons of the reference and targets")
        self.blast_comparaisons = self.record_list2blast([self.ref_subset["subset_record"]] + [i["subset_record"] for i in self.target_subrecords], min_identity)


        # create a combined fasta with all sequences
        records = chain.from_iterable([target_sub.faa_record for target_sub in self.target_faa_subrecords] + [self.ref_subrecord.faa_record])
        SeqIO.write(records, f'{self.output_folder}/combined_sub.faa', "fasta")

        logging.info("Clustering protein sequences of all records combined")
        
        self.clusters, self.locus2cluster = self.cluster_protein_sequences(f'{self.output_folder}/combined_sub.faa')
        self.cluster2annot = {}


    def record_list2blast(self, record_list, min_identity):
        '''
        :param record_list:
        :return:
        '''

        blast_result_files = []

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

        
    def parse_basename2label(self, path):
        import pandas 
        df = pandas.read_csv(path, sep="\t",names=["basename", "label"])
        self.basename2label = {row.basename: row.label for n, row in df.iterrows()}
        



    def blast_best_hit(self, query_file, db_file):
        """
        Formats a database file using makeblastdb, performs a blastp search of the input query file
        against the database, extracts the best hit, and returns it.

        Args:
            query_file (str): Path to the input query FASTA file.
            db_file (str): Path to the input database FASTA file.

        Returns:
            dict: A dictionary containing information about the best blast hit, including the
            query ID, subject ID, percent identity, e-value, and alignment length.
        """
        # Create a temporary file for the output BLAST results
        db_basename = os.path.basename(os.path.splitext(db_file)[0])
        blast_output_file = f'{self.output_folder}/{db_basename}.m6'

        # Format the database file using makeblastdb
        try:
            subprocess.run(["makeblastdb", "-in", db_file, "-dbtype", "prot"])
        except subprocess.CalledProcessError as e:
            raise Exception(f"Error executing makeblastdb: {e.stderr.decode().strip()}")
        # Execute the blastp search
        try:
            subprocess.run(["blastp", "-query", query_file, "-db", db_file, "-out", blast_output_file, "-outfmt", "6"])
        except subprocess.CalledProcessError as e:
            raise Exception(f"Error executing blastp: {e.stderr.decode().strip()}")
        # Parse the BLAST results and extract the best hit
        with open(blast_output_file, "r") as blast_output:
            data = [i.strip().split("\t") for i in blast_output.readlines()]          
        
            bbh =  {"query_id": data[0][0], 
                    "subject_id": data[0][1],
                    "percent_identity": float(data[0][2]),
                    "evalue": float(data[0][10]), 
                    "alignment_length": int(data[0][3])}
            
            tmp_files = glob.glob(f'{self.output_folder}/{db_basename}.faa.p*')
            for tmp in tmp_files:
                os.remove(tmp)
            
            return bbh


    def get_sequence_subset(self, 
                            seqrecords, 
                            target_id, 
                            upstream_bp, 
                            downstream_bp, 
                            basename=False,
                            flip_record=True):
        """
        Reads a GenBank file, searches for a specific locus_tag or protein_id in the record,
        and returns a subset of the record with the matching feature and a specified number of
        bp upstream and downstream of the target feature.

        Args:
            genbank_file (str): Path to the input GenBank file.
            target_id (str): The target locus_tag or protein_id to search for.
            upstream_bp (int): Number of base pairs upstream of the target feature to include in the subset.
            downstream_bp (int): Number of base pairs downstream of the target feature to include in the subset.

        Returns:
            Bio.SeqRecord.SeqRecord: A new SeqRecord object containing the sequence subset.
        """

        # Search for the feature with the specified locus_tag or protein_id
        feature = None
        flip = False
        for record in seqrecords:
            for f in record.features:
                if "protein_id" in f.qualifiers and target_id == f.qualifiers["protein_id"][0]:
                    feature = f
                    match_record = record
                    seq = feature.extract(record.seq).translate()
                    break
                if "locus_tag" in f.qualifiers and target_id == f.qualifiers["locus_tag"][0]:
                    feature = f
                    match_record = record
                    seq = feature.extract(record.seq).translate()
                    break
        
        if feature is None:
            raise ValueError(f"Target ID '{target_id}' not found in record.")

        # Determine the start and end positions of the feature
        feature_start = feature.location.start.position
        feature_end = feature.location.end.position
        strand = feature.location.strand
        # Calculate the start and end positions of the desired sequence subset
        subset_start = max(0, feature_start - upstream_bp)
        subset_end = min(len(match_record), feature_end + downstream_bp)

        # Create a new SeqRecord object for the sequence subset
        subset_record = match_record[subset_start:subset_end]

        if flip_record:
            if strand != self.reference_orientation:
                subset_record = subset_record.reverse_complement(id=subset_record.id,
                                                                 name=subset_record.id,
                                                                 description=subset_record.description,
                                                                 annotations=True)

        # track mapping of contig names with file name
        # => allow renaming of tracks
        self.contig2basename[subset_record.id] = basename

        return {"subset_record": subset_record, 
                "seq":seq, 
                "basename": basename,
                "feature": feature, 
                "subset_start": subset_start,
                "subset_end": subset_end,
                "target_locus_strand": strand}

    def get_colors_from_clusters(self,):
        #self.clusters, self.locus2cluster
        
        keep = []
        for cluster in self.clusters:
            if len(self.clusters[cluster]) > 1:
                keep.append(cluster)
        n_clusters = len(keep)
        logging.info(f"Non singletons clusters: {n_clusters}")
        colors = self.get_spaced_colors(n_clusters)
        
        self.cluster2color = {cluster:color for cluster,color in zip(keep,colors)}
        

    def get_spaced_colors(self, n):
        
        '''
        Return n spaces colors color
        '''

        max_value = 16581375 #255**3
        interval = int(max_value / n)
        colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

        return ['#%02x%02x%02x' % (int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]


    def get_genome_label(self, contig_name):
        '''
        get genome label from contig label
        WARNING: require unique contig names 
        TODO: make it work with non unique contig name
        '''
        if self.contig2basename[contig_name] in self.basename2label:
            label = self.basename2label[self.contig2basename[contig_name]]
        else:
            label = self.contig2basename[contig_name]
        return label

    def add_gbk_track(self, 
                      gv, 
                      gbk, 
                      reg2color = {},
                      color_clusters=False,
                      count=0,
                      highligh_line=True):
        import re

        label = self.get_genome_label(gbk.id)

        print(f"ADDING {label}")
        track = gv.add_feature_track(label, 
                                     len(gbk.seq))

        for feature in gbk.features:
            if feature.type == 'CDS':
                # default
                color='grey'
                labelcolor = 'black'
                if 'gene' in feature.qualifiers:
                    label = feature.qualifiers["gene"][0]
                    
                else:
                    label = ""
                # either color by cluster or by dict data
                if color_clusters:
                    if "locus_tag" in feature.qualifiers:
                        seqid = feature.qualifiers["locus_tag"][0]
                    elif "protein_id" in feature.qualifiers:
                        seqid = feature.qualifiers["protein_id"][0]
                    else: continue 
                    
                    if seqid in  self.locus2cluster:
                        cluster = self.locus2cluster[seqid]
                    else: continue
                    
                    if cluster in self.cluster2color:
                        color = self.cluster2color[cluster]
                        
                        
                    # color and label specific clusters based on color file 
                    if cluster in self.cluster2highlights:
                        color = self.cluster2highlights[cluster]["color"]
                        label = self.cluster2highlights[cluster]["label"]
                        labelcolor = 'red'
                # color particular locus
                linecolor = color
                for reg in reg2color:
                    if 'product' in feature.qualifiers:
                        m = re.match(reg, feature.qualifiers["product"][0])
                        if m:
                            if reg2color[reg]["label"]:
                                label = reg2color[reg]["label"]
                            else:
                                label = m.group(0)
                            linecolor = reg2color[reg]["color"]
                            break 
                    if 'gene' in feature.qualifiers:
                        m = re.match(reg, feature.qualifiers["gene"][0])
                        if m:
                            if reg2color[reg]["label"] is not None:
                                label = reg2color[reg]["label"]
                            else:
                                label = m.group(0)
                            linecolor = reg2color[reg]["color"]
                            break 
      
                # if we have a match, increase linewidth and set line color
                if highligh_line and linecolor != color:
                    linecolor = linecolor
                    linewidth = 3
                else:
                    color = linecolor
                    linecolor = None 
                    linewidth = 0

                track.add_feature(feature.location.start, 
                                  feature.location.end, 
                                  feature.location.strand, 
                                  label=label, 
                                  labelrotation=45, 
                                  facecolor=color,
                                  labelcolor = labelcolor,
                                  linewidth = linewidth,
                                  edgecolor=linecolor)
            
            if feature.type in ["tRNA", 'rRNA']:                
                track.add_feature(feature.location.start, 
                                feature.location.end, 
                                feature.location.strand, 
                                label=feature.type, 
                                linewidth=0, 
                                labelrotation=45, 
                                facecolor='orange')
    
    
    def add_link(self, gv, blast_table):
        import pandas
        names = ["qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
        blast_data = pandas.read_csv(blast_table, sep="\t", names=names)
        for n,link in blast_data.iterrows():
            qseqid = self.contig2basename[link.qseqid]
            sseqid =  self.contig2basename[link.sseqid]
            
            if qseqid in self.basename2label:
                query_label = self.basename2label[qseqid]
            else:
                query_label = qseqid
            if sseqid in self.basename2label:
                subject_label = self.basename2label[sseqid]
            else:
                subject_label = sseqid
            print()
            link_data1 = (query_label, link.qstart, link.qend)
            link_data2 = (subject_label, link.sstart, link.send)
            gv.add_link(link_data1, link_data2, v=link.pident, curve=True,normal_color="skyblue", inverted_color="lime")

    def parse_highlights(self, path):
        import pandas
        df = pandas.read_csv(path, sep="\t")
        
        
        for n,row in df.iterrows():
            print(row)
            cluster = self.locus2cluster[row["locus"]]
            self.cluster2highlights[cluster] = {"color": row["color"], "label": row["label"]}
        


    def draw_alignment(self, outname, color_clusters=False, reg_highlights={}, highlight_line=True):
        gv = GenomeViz(
            fig_track_height=0.7,
            feature_track_ratio=0.2,
            tick_track_ratio=0.4,
            tick_style="bar",
            align_type="center",
        )

        if color_clusters:
            self.get_colors_from_clusters()

        print('Adding ref track...........')
        self.add_gbk_track(gv, 
                           self.ref_subrecord.gbk_record, 
                           color_clusters=color_clusters, 
                           reg2color=reg_highlights)
        
        print('Adding target tracks...........')
        for n,sub in enumerate(self.target_faa_subrecords):
            print(f"######### Track {n}")
            self.add_gbk_track(gv, sub.gbk_record, color_clusters=color_clusters, reg2color=reg_highlights, count=n)
        
        for blast in self.blast_comparaisons:
            self.add_link(gv, blast)
        
        fig = gv.plotfig()
        gv.set_colorbar(fig, vmin=self.min_identity, bar_height=0.4,bar_colors=["skyblue", "lime"] )
        fig.savefig(outname)

    def get_cluster_consensus_annot(self, records):
        from collections import Counter

        cluster2n_genomes = {}
        for record in records:
            clusters = set()
            for feature in record.features:
                gene = '-'
                product = '-'
                if "locus_tag" in feature.qualifiers:
                    seqid = feature.qualifiers["locus_tag"][0]
                elif "protein_id" in feature.qualifiers:
                    seqid = feature.qualifiers["protein_id"][0]
                else: continue 
                if 'gene' in feature.qualifiers:
                    gene = feature.qualifiers["gene"][0]
                if 'product' in feature.qualifiers:
                    product = feature.qualifiers["product"][0]
                if seqid in  self.locus2cluster:
                    cluster = self.locus2cluster[seqid]
                    clusters.add(cluster)
                    if cluster not in self.cluster2annot:
                        self.cluster2annot[cluster] = {"gene": [], "product": []}
                    self.cluster2annot[cluster]["gene"].append(gene)
                    self.cluster2annot[cluster]["product"].append(product)
            for cluster in clusters:
                if cluster not in cluster2n_genomes:
                    cluster2n_genomes[cluster] = 1
                else:
                    cluster2n_genomes[cluster] += 1
        
        self.cluster2annot_consensus = {}
        for cluster in self.clusters:
            gene_c = dict(Counter(self.cluster2annot[cluster]["gene"]))
            product_c = dict(Counter(self.cluster2annot[cluster]["product"]))
            self.cluster2annot_consensus[cluster] = {"freq": cluster2n_genomes[cluster], "gene": gene_c, "product": product_c }
        


    def extract_clusters(self, record):

        clusters = set()
        
        for feature in record.features:
            gene = '-'
            product = '-'
            if "locus_tag" in feature.qualifiers:
                seqid = feature.qualifiers["locus_tag"][0]
            elif "protein_id" in feature.qualifiers:
                seqid = feature.qualifiers["protein_id"][0]
            else: continue 
            if seqid in  self.locus2cluster:
                cluster = self.locus2cluster[seqid]
                clusters.add(cluster)
            else: continue
        return clusters

    def compare_neigborhood(self,):
        '''
        Determine the number of shared OG between pairs of genome
        Calculate jaccard distance
        '''
        from itertools import combinations
        # {"subset_record": subset_record, 
        #        "seq":seq, 
        #        "basename": basename,
        #        "feature": feature, 
        #        "subset_start": subset_start,
        #        "subset_end": subset_end,
        #        "target_locus_strand": strand}
        all_records = [self.ref_subset] + self.target_subrecords
        # get consensus cluster annotation
        self.get_cluster_consensus_annot([i["subset_record"] for i in all_records])
        with open(f"{self.output_folder}/cluster_summary.tsv", "w") as f:
            f.write("cluster_id\tfreq\tgene\tproduct\n")
            for cluster in self.cluster2annot_consensus:
                f.write(f'{cluster}\t{self.cluster2annot_consensus[cluster]["freq"]}\t{self.cluster2annot_consensus[cluster]["gene"]}\t{self.cluster2annot_consensus[cluster]["product"]}\n')

        # {"basename": basename, "records": records})
        with open(f"{self.output_folder}/cluster_comparaison.tsv", "w") as f:
            f.write("genome_1\tgenome_2\tn_clusters_1\tn_clusters_2\tn_shared\tn_unique_1\tn_unique_2\tjaccard\n")
            for pair in combinations(all_records, 2):
                genome_label_1 = self.get_genome_label(pair[0]["subset_record"].id)
                genome_label_2 = self.get_genome_label(pair[1]["subset_record"].id)
                record_1 = pair[0]["subset_record"]
                record_2 = pair[1]["subset_record"]

                clusters_1 = self.extract_clusters(record_1)
                clusters_2 = self.extract_clusters(record_2)

                shared = clusters_1.intersection(clusters_2)
                unique_1 = clusters_1.difference(clusters_2)
                unique_2 = clusters_2.difference(clusters_1)
                jaccard = len(shared) / (len(shared) + len(unique_1) + len(unique_2))
                f.write(f"{genome_label_1}\t{genome_label_2}\t{len(clusters_1)}\t{len(clusters_2)}\t{len(shared)}\t{len(unique_1)}\t{len(unique_2)}\t{1-jaccard}\n")



    def cluster_protein_sequences(self, 
                                  faa_file, 
                                  min_identity=0.5, 
                                  cluster_mode="0"):
        """
        Clusters protein sequences from a GenBank file using mmseqs.

        Args:
            faa_file (str): Path to the input faa file.
            min_identity (float): Minimum sequence identity threshold for clustering (default: 0.9).
            cluster_mode (str): Clustering mode for mmseqs (default: 'cluster').

        Returns:
            dict: A dictionary where the keys are cluster IDs and the values are lists of protein sequences.
        """

        try:
            cmd = f"mmseqs easy-cluster {faa_file} {self.output_folder}/mmseqs {faa_file}.tmp --min-seq-id {min_identity} --cluster-mode {cluster_mode}"
            subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            raise Exception(f"Error executing mmseqs: {e.stderr.decode().strip()}")

        clusters = {}
        locus2cluster = {}
        with open(f"{self.output_folder}/mmseqs_cluster.tsv") as f:
            for line in f:
                cluster_id, seq_id = line.strip().split("\t")
                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append(seq_id)
                locus2cluster[seq_id] = cluster_id

        #os.remove(tmp_file_handle.name)
        #os.remove(f"{tmp_file_handle.name}.tmp")

        return clusters, locus2cluster




if __name__ == '__main__':
    import argparse
    import sys
    import re
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-a",'--input table', type=str, help="input table")
    parser.add_argument("-l",'--locus', type=str, help="locus_tag")
    parser.add_argument("-r",'--reference', type=str, help="reference genbank")
    parser.add_argument("-t",'--targets', type=str, help="target genbank(s)", nargs='+')
    parser.add_argument("-o",'--output_name', help="output name", default='out')
    parser.add_argument("-v",'--svg', help="output svg rather than pdf", action='store_true')
    parser.add_argument("-ls",'--left_side_window', type=int, help="left siden window", default=15000)
    parser.add_argument("-rs",'--right_side_window', type=int, help="right side window", default=15000)
    parser.add_argument("-i",'--min_identity', type=int, help="minimum identity for blast", default=50)
    parser.add_argument("-sl",'--show_labels', action="store_false", help="do not show show labels")
    parser.add_argument("-s", '--samtools_depth', default=False, help="add depth plot from samtool depth (only for the last query). Should match the chromosome/contig names of the gbk.")
    parser.add_argument("-x",'--tblastx', action="store_true", help="execute tblastx and not blasn (6 frame translations)")
    parser.add_argument("-g",'--gc_plot', action="store_true", help="Show GC plot")
    parser.add_argument("-f", '--force', action="store_true", help="Don't prompt before output folder removal")
    parser.add_argument("-c",'--color', type=str, help="color table")
    parser.add_argument("-b",'--basename', type=str, help="file basename 2 label")

    args = parser.parse_args()

    if args.reference in args.targets:
        print("Removing reference from targets")
        args.targets.pop(args.targets.index(args.reference))
    
    G = GenomeComp(args.locus,
                   args.reference,
                   args.targets,
                   upstream_bp=args.right_side_window,
                   downstream_bp=args.left_side_window,
                   tblastx=args.tblastx,
                   output_name='out',
                   svg=False,
                   min_identity=args.min_identity,
                   force_data_dir=args.force)

    if args.color:
        G.parse_highlights(args.color)
    if args.basename:
        G.parse_basename2label(args.basename)

    reg_highlights = {r'IS[0-9]+' : { "color": 'blue', "label": None},
                    r'OXA-[0-9]+': { "color": 'red', "label": None},
                    r'CTX-M-[0-9]+': { "color": 'red', "label": None},
                    r'TEM-[0-9]+': { "color": 'red', "label": None},
                    r'NDM-[0-9]+': { "color": 'red', "label": None},
                    r'blaNDM-[0-9]+': { "color": 'red', "label": None},
                    r'mecA': { "color": 'red', "label": None},
                    r'bla.*': { "color": 'red', "label": None},
                    r'transposase': { "color": 'blue', "label": 'transp.'},
                    r'lactamase': { "color": 'red', "label": 'BL'},
                    #r'hypothetical': { "color": 'black', "label": ''},
                    }

    G.draw_alignment("test.pdf", 
                     color_clusters=True,
                     reg_highlights=reg_highlights,
                     highlight_line=True)

    G.compare_neigborhood()