import sys
import os
import subprocess
import argparse
import traceback

try:
    from Bio import SeqIO
except ImportError as e:
    print(e)
    sys.exit("ERROR...unable to load required python modules\ntry reinstalling and/or activating TIDAL environment:\n\tconda env create -f TIDAL.yml --name TIDAL\n\tconda activate TIDAL\n")


class Insertion:
    def __init__(self):
        self.info = {}
        self.is_ref = False

ANNOTATION_FILES = {
    "consensus" : "Tidalbase_transposon_gene_sequence.fa",
    "annotation" : "refflat_dm6.txt",
    "gem" : "gem_mappability_dm6_100mer.mappability.gz",
    "virus" : "fly_virus_structure_repbase.fa",
    "repeatmasker" : "repmasker_dm6_track.txt",
    "table" : "Tidalbase_Dmel_TE_classifications_2015.txt"
}

def main():
    code_dir = os.path.dirname(os.path.abspath(__file__))
    annotation_dir = code_dir+"/../annotation/"
    print(code_dir)
    args = parse_args(annotation_dir)
    args = setup_input_files(args)
    chrom_len_file = make_chrom_files(args)
    read_length = estimate_read_length(args.fastq)
    depletion_tbl, insertion_tbl  = run_tidal(args, chrom_len_file, read_length, code_dir)
    consensus_repmask = repeatmask(args.reference, args.consensus, args.processors, args.out)
    write_output(depletion_tbl, insertion_tbl, consensus_repmask, args.table, args.reference, args.sample_name, args.out)


def parse_args(annot_dir):
    parser = argparse.ArgumentParser(prog='TIDAL', description="Pipeline built to identify Transposon Insertion and Depletion in flies")

    ## required ##
    parser.add_argument("-f", "--fastq", type=str, help="A FASTQ file containing the NGS reads", required=True)
    parser.add_argument("-r", "--reference", type=str, help="A reference genome sequence in fasta format", required=True)
    parser.add_argument("-m", "--masked", type=str, help="A reference genome sequence in fasta format masked by RepeatMasker", required=True)

    parser.add_argument("-c", "--consensus", type=str, help="The consensus sequences of the TEs in fasta format", required=False)
    parser.add_argument("-a", "--annotation", type=str, help="The RefSeq annotation from UCSC genome browser", required=False)
    parser.add_argument("-g", "--gem", type=str, help="The gem mappability file for use with FREEC", required=False)
    parser.add_argument("-v", "--virus", type=str, help="Manunally curated sequence from fly viruses, structural and repbase sequences", required=False)
    parser.add_argument("-n", "--repeatmasker", type=str, help=" Repeat masker track from UCSC genome browser (table browser, track: Repeatmasker, table: rmsk, output format: all fields from table)", required=False)
    parser.add_argument("-t", "--table", type=str, help="Custom table for repbase to flybase lookup", required=False)
    parser.add_argument("-p", "--processors", type=int, help="Number of CPU threads to use for compatible tools (default=1)", required=False)
    parser.add_argument("-s", "--sample_name", type=str, help="Sample name to use for output files (default=FASTQ name)", required=False)
    parser.add_argument("-o", "--out", type=str, help="Directory to create the output files (must be empty/creatable)(default='.')", required=False)

    args = parser.parse_args()

    args.fastq = get_abs_path(args.fastq)
    args.reference = get_abs_path(args.reference)
    args.masked = get_abs_path(args.masked)
    
    if args.consensus is None:
        args.consensus = get_abs_path(annot_dir+ANNOTATION_FILES['consensus'])
    else:
        args.consensus = get_abs_path(args.consensus)


    if args.annotation is None:
        args.annotation = get_abs_path(annot_dir+ANNOTATION_FILES['annotation'])
    else:
        args.annotation = get_abs_path(args.annotation)
    
    if args.gem is None:
        args.gem = get_abs_path(annot_dir+ANNOTATION_FILES['gem'])
    else:
        args.gem = get_abs_path(args.gem)

    if args.virus is None:
        args.virus = get_abs_path(annot_dir+ANNOTATION_FILES['virus'])
    else:
        args.virus = get_abs_path(args.virus)

    if args.table is None:
        args.table = get_abs_path(annot_dir+ANNOTATION_FILES['table'])
    else:
        args.table = get_abs_path(args.table)

    if args.sample_name is None:
        args.sample_name = os.path.basename(args.fastq)
        args.sample_name = ".".join(args.sample_name.split(".")[:-1])
    else:
        if "/" in args.sample_name:
            sys.exit(args.sample_name+" is not a valid sample name...\n")
    
    if args.processors is None:
        args.processors = 1

    if args.out is None:
        args.out = os.path.abspath(".")
    else:
        args.out = os.path.abspath(args.out)
        try:
            mkdir(args.out)
        except Exception as e:
            track = traceback.format_exc()
            print(track, file=sys.stderr)
            print("cannot create output directory: ",args.out,"exiting...", file=sys.stderr)
            sys.exit(1)

    if args.repeatmasker is None:
        args.repeatmasker = get_abs_path(annot_dir+ANNOTATION_FILES['repeatmasker'])
        # args.repeatmasker = get_abs_path(args.out+"/repeatmasker.out")
    else:
        args.repeatmasker = get_abs_path(args.repeatmasker)


    return args


def mkdir(indir, log=None):
    if os.path.isdir(indir) == False:
        os.mkdir(indir)

def remove_file_ext(infile):
    no_ext = ".".join(infile.split(".")[:-1])
    return no_ext

def get_abs_path(in_file, log=None):
    if os.path.isfile(in_file):
        return os.path.abspath(in_file)
    else:
        msg = " ".join(["Cannot find file:",in_file,"exiting....\n"])
        sys.stderr.write(msg)
        sys.exit(1)

def get_base_name(path):
    no_path = os.path.basename(path)
    if no_path[-3:] == ".gz":
        no_path = no_path[:-3]
    no_ext = ".".join(no_path.split(".")[:-1])

    return no_ext

def copy(infile, outdir, outfilename=None):
    if outfilename is None:
        outfile = outdir+"/"+os.path.basename(infile)
    else:
        outfile = outdir+"/"+outfilename
    subprocess.call(["cp", infile, outfile])
    return outfile

def run_command(cmd_list, log=None, fatal=True):
    msg = ""
    if log is None:
        try:
            subprocess.check_call(cmd_list)
        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            sys.stderr.write(msg)
            if fatal:
                sys.exit(1)
            else:
                return False
    
    else:
        try:
            out = open(log,"a")
            out.write(" ".join(cmd_list)+"\n")
            subprocess.check_call(cmd_list, stdout=out, stderr=out)
            out.close()

        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            writelog(log, msg)
            sys.stderr.write(msg)
            if fatal:
                sys.exit(1)
            else:
                return False
    
    return True

def run_command_stdout(cmd_list, out_file, log=None, fatal=True):
    msg = ""
    if log is None:
        try:
            # print(" ".join(cmd_list)+" > "+out_file)
            out = open(out_file,"w")
            subprocess.check_call(cmd_list, stdout=out)
            out.close()
        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            sys.stderr.write(msg)
            if fatal:
                sys.exit(1)
            else:
                return False
    
    else:
        try:
            out_log = open(log,"a")
            out_log.write(" ".join(cmd_list)+" > "+out_file+"\n")
            out = open(out_file,"w")
            subprocess.check_call(cmd_list, stdout=out, stderr=out_log)
            out.close()
            out_log.close()

        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            writelog(log, msg)
            sys.stderr.write(msg)
            if fatal:
                sys.exit(1)
            else:
                return False
    
    return True

def writelog(log, msg):
    if log is not None:
        with open(log, "a") as out:
            out.write(msg)

def setup_input_files(args):
    ref_dir = args.out+"/reference"
    mkdir(ref_dir)
    os.chdir(ref_dir)
    args.reference = copy(args.reference, ref_dir)
    args.consensus = copy(args.consensus, ref_dir)
    args.annotation = copy(args.annotation, ref_dir)
    args.gem = copy(args.gem, ref_dir)
    if args.gem[-3:] == ".gz":
        run_command(["gunzip", args.gem])
        args.gem = args.gem[:-3]
    args.virus = copy(args.virus, ref_dir)
    args.masked = copy(args.masked, ref_dir, outfilename="masked."+ get_base_name(args.masked))
    args.repeatmasker = copy(args.repeatmasker, ref_dir)
    args.table = copy(args.table, ref_dir)

    run_command(["bowtie-build", args.virus, get_base_name(args.virus)])

    run_command(["bowtie-build", args.consensus, get_base_name(args.consensus)])
    run_command(["bowtie2-build", args.consensus, get_base_name(args.consensus)])

    run_command(["bowtie-build", args.reference, get_base_name(args.reference)])
    run_command(["bowtie2-build", args.reference, get_base_name(args.reference)])
    run_command(["bowtie-build", args.masked, get_base_name(args.masked)])

    mkdir(args.out+"/TIDAL_out")
    if args.fastq[-3:] == ".gz":
        args.fastq = copy(args.fastq, args.out+"/TIDAL_out", outfilename=args.sample_name+".fastq.gz")
        run_command(["gunzip", args.out+"/TIDAL_out/"+args.sample_name+".fastq.gz"])
        args.fastq = args.fastq[:-3]
    else:
        args.fastq = copy(args.fastq, args.out+"/TIDAL_out", outfilename=args.sample_name+".fastq")
    
    return args


def make_chrom_files(args):
    mkdir(args.out+"/reference/chroms")
    chrom_len_file = args.out+"/reference/chrom_lengths.tsv"
    with open(args.reference, "r") as infa, open(chrom_len_file,"w") as out_lens:
        for record in SeqIO.parse(infa, "fasta"):
            seq_name = str(record.id)
            seq_len = len(str(record.seq))
            out_lens.write(seq_name+"\t"+str(seq_len)+"\n")

            chrom_file = args.out+"/reference/chroms/"+seq_name+".fa"
            with open(chrom_file,"w") as outfa:
                outfa.write(">"+seq_name+"\n")
                outfa.write(str(record.seq)+"\n")
    
    return chrom_len_file


def estimate_read_length(fq, reads=10000):
    lengths = []
    with open(fq,"r") as f:
        for x, line in enumerate(f):
            if x%4 == 1:
                lengths.append(len(line.replace('\n',"")))
            
            if x >= reads:
                break
    
    length = sum(lengths)//len(lengths)

    return length

def repeatmask(reference, consensus, threads, out):
    tmp_dir = out+"/rm/"
    os.mkdir(tmp_dir)
    run_command(["cp", reference, tmp_dir+"/reference.fasta"])
    
    te_consensus = tmp_dir+"consensus_tes.fasta"
    with open(te_consensus, "w") as outfa, open(consensus, "r") as infa:
        for record in SeqIO.parse(infa, "fasta"):
            if "|TE@" in str(record.id):
                outfa.write(">"+str(record.id)+"\n")
                outfa.write(str(record.seq)+"\n")


    command = ["RepeatMasker","-pa", str(threads), "-lib", te_consensus, "-dir", tmp_dir, "-s", "-nolow", "-no_is", reference]
    run_command(command)

    rm_out = ""
    for f in os.listdir(tmp_dir):
        if ("fasta.out" in f and f[-9:] == "fasta.out") or ("fa.out" in f and f[-6:] == "fa.out"):
            rm_out = tmp_dir+"/"+f
    
    if rm_out == "":
        sys.exit("can't find Repeatmasker output in:"+tmp_dir+"\n")
    

    rm_formatted = out+"/rm/repmasker_track.txt"
    with open(rm_formatted,"w") as out, open(rm_out, "r") as inf:
        header = ["#bin","swScore","milliDiv","milliDel","milliIns","genoName","genoStart","genoEnd","genoLeft","strand","repName","repClass","repFamily","repStart","repEnd","repLeft","id"]
        out.write("\t".join(header)+"\n")

        for ln,line in enumerate(inf):
            if ln > 2:
                line = line.replace("\n","")
                vals = []
                split_line = line.split(" ")
                for val in split_line:
                    if val != "":
                        vals.append(val)

                genoLeft = vals[7]
                if "(" in genoLeft:
                    genoLeft = genoLeft.replace("(","")
                    genoLeft = genoLeft.replace(")","")
                    genoLeft = "-"+genoLeft

                strand = vals[8]
                if strand == "C":
                    strand = "-"

                te_name = vals[9]
                te_name = te_name.split("@")[-1]
                te_name = te_name.split("=")[0]

                repStart = vals[11]
                if "(" in repStart:
                    repStart = repStart.replace("(","")
                    repStart = repStart.replace(")","")
                    repStart = "-"+repStart
                
                repEnd = vals[12]
                if "(" in repEnd:
                    repEnd = repEnd.replace("(","")
                    repEnd = repEnd.replace(")","")
                    repEnd = "-"+repEnd

                repLeft = vals[13]
                if "(" in repLeft:
                    repLeft = repLeft.replace("(","")
                    repLeft = repLeft.replace(")","")
                    repLeft = "-"+repLeft

                vals[7] = genoLeft
                vals[8] = strand
                vals[9] = te_name
                vals[11] = repStart
                vals[12] = repEnd
                vals[13] = repLeft

                out_line = ["0", vals[0], vals[1], vals[2], vals[3], vals[4], vals[5], vals[6], vals[7], vals[8], vals[9], vals[10], vals[10], vals[11], vals[12], vals[13], vals[14]]

                out.write("\t".join(out_line)+"\n")

    return rm_formatted

def run_tidal(args, chrom_len, read_length, codedir):
    os.chdir(args.out+"/TIDAL_out")

    run_command(["bash", codedir+"/data_prep.sh", os.path.basename(args.fastq), codedir, str(args.processors)])

    if not os.path.exists(os.path.basename(args.fastq)+".uq.polyn"):
        sys.exit("data_prep.sh failed...exiting...\n")

    run_command([
        "bash", codedir+"/insert_pipeline.sh",
            os.path.basename(args.fastq)+".uq.polyn",
            str(read_length),
            codedir,
            remove_file_ext(args.reference),
            remove_file_ext(args.masked),
            remove_file_ext(args.consensus),
            args.reference,
            args.annotation,
            chrom_len,
            args.out+"/reference/chroms/",
            args.gem,
            remove_file_ext(args.virus),
            args.consensus,
            str(args.processors)
    ])

    if not os.path.exists(args.out+"/TIDAL_out/insertion/"+args.sample_name+"_Inserts_Annotated.txt"):
        sys.exit("insert_pipeline.sh failed...exiting...\n")

    run_command(["bash", codedir+"/setup.sh", os.path.basename(args.fastq)])

    run_command([
        "bash", codedir+"/depletion_pipeline.sh",
            os.path.basename(args.fastq)+".uq.polyn",
            str(read_length),
            codedir,
            remove_file_ext(args.reference),
            remove_file_ext(args.masked),
            remove_file_ext(args.consensus),
            args.reference,
            args.masked,
            args.repeatmasker,
            args.annotation,
            args.table,
            chrom_len,
            str(args.processors)
    ])

    if not os.path.exists(args.out+"/TIDAL_out/depletion/"+args.sample_name+"_Depletion_Annotated.txt"):
        sys.exit("insert_pipeline.sh failed...exiting...\n")

    run_command(["bash", codedir+"/last_part.sh", os.path.basename(args.fastq), codedir])

    if not os.path.exists(args.out+"/TIDAL_out/"+args.sample_name+"_result/"+args.sample_name+"_Depletion_Annotated.txt"):
        sys.exit("last_part.sh failed...exiting...\n")

    depletion_table = args.out+"/TIDAL_out/"+args.sample_name+"_result/"+args.sample_name+"_Depletion_Annotated_TEonly.txt"
    insertion_table = args.out+"/TIDAL_out/"+args.sample_name+"_result/"+args.sample_name+"_Inserts_Annotated.txt"

    return depletion_table, insertion_table 


def read_table(intbl):
    insertions = []
    with open(intbl,"r") as tbl:
        for x,line in enumerate(tbl):
            if x == 0:
                headers = line.replace("\n","").split("\t")
            else:
                insert = Insertion()
                split_line = line.replace("\n","").split("\t")
                for x,val in enumerate(split_line):
                    if headers[x] == "TE":
                        val = val.split("@")[-1]
                        val = val.split("=")[0]
                    insert.info[headers[x]] = val
                
                insertions.append(insert)
    
    return insertions

def get_repbase_families(table):
    families = {}
    
    with open(table,"r") as tbl:
        for x,line in enumerate(tbl):
            split_line = line.split("\t")
            families[split_line[7]] = split_line[9]
    
    return families
    

def write_output(depletion_tbl, insertion_tbl, rm_out, rep_fly_table, reference, sample_name, out_dir):

    # make insertions bed
    insertions = read_table(insertion_tbl)
    tmp_bed = out_dir+"/TIDAL_out/"+sample_name+"_result/"+sample_name+"_TIDAL_tmp.bed"
    with open(tmp_bed, "w") as outbed:
        for x,insert in enumerate(insertions):
            if not insert.is_ref:
                name = insert.info['TE']+"|non-reference|"+insert.info['Coverage_Ratio']+"|"+sample_name+"|sr|"
                out_line = [insert.info["Chr"], insert.info["Chr_coord_5p"], insert.info["Chr_coord_3p"], name, "0", "."]

                outbed.write("\t".join(out_line) + "\n")
    

    all_repeats = read_table(rm_out)
    depletions = read_table(depletion_tbl)
    # families = get_repbase_families(rep_fly_table)

    chroms = []
    with open(reference, "r") as infa:
        for record in SeqIO.parse(infa, "fasta"):
            chroms.append(str(record.id))


    # make all ref TE bed
    ref_insert_bed = out_dir+"/TIDAL_out/"+sample_name+"_result/"+sample_name+"_TIDAL_ref.bed"
    with open(ref_insert_bed,"w") as ref_bed:
        for repeat in all_repeats:
            if repeat.info['genoName'] in chroms:
                chrom = repeat.info['genoName']
                start = repeat.info['genoStart']
                end = repeat.info['genoEnd']
                name = repeat.info['repName']+"|reference|NA|"+sample_name+"|sr|"
                out_line = [chrom, start, end, name, "0", repeat.info['strand']]
                ref_bed.write("\t".join(out_line) + "\n")
    

    # make TE depletion bed
    depletion_bed = out_dir+"/TIDAL_out/"+sample_name+"_result/"+sample_name+"_TIDAL_depletion.bed"
    with open(depletion_bed, "w") as dep_bed:
        for depletion in depletions:
            chrom = depletion.info['Chr_5p']
            start = depletion.info['Chr_coord_5p_start']
            end = depletion.info['Chr_coord_3p_end']
            name = depletion.info['repName']+"|depletion|NA|"+sample_name+"|sr|"
            out_line = [chrom, start, end, name, "0", "."]
            dep_bed.write("\t".join(out_line) + "\n")
    

    # remove ref TEs that are depleted
    nonabs_bed = out_dir+"/TIDAL_out/"+sample_name+"_result/"+sample_name+"_TIDAL_nonabs.bed"
    run_command_stdout(["bedtools", "intersect", "-v", "-a", ref_insert_bed, "-b", depletion_bed], nonabs_bed)

    # combine nonabs and nonref TEs
    with open(nonabs_bed, "r") as inbed, open(tmp_bed, "a") as outbed:
        for line in inbed:
            outbed.write(line)
    
    # get uniq predictions
    uniq_inserts = {}
    with open(tmp_bed, "r") as inbed:
        for line in inbed:
            split_line = line.split("\t")
            te_type = split_line[3].split("|")[1]
            key = "_".join([split_line[0], split_line[1], split_line[2], te_type])
            if key not in uniq_inserts.keys() or te_type == "reference":
                uniq_inserts[key] = split_line
            else:
                ratio = float(split_line[3].split("|")[2])
                existing_ratio = float(uniq_inserts[key][3].split("|")[2])
                if ratio > existing_ratio:
                    uniq_inserts[key] = split_line

    unsorted_bed = out_dir+"/TIDAL_out/"+sample_name+"_result/"+sample_name+"_TIDAL_nonredundant_unsorted.bed"
    with open(unsorted_bed, "w") as outbed:
        for x,key in enumerate(uniq_inserts.keys()):
            split_line = uniq_inserts[key]
            split_line[3] = split_line[3]+str(x)
            outbed.write("\t".join(split_line)+"\n")


    nonredundant_bed = out_dir+"/TIDAL_out/"+sample_name+"_result/"+sample_name+"_TIDAL_nonredundant_unlabeled.bed"
    run_command_stdout(["bedtools", "sort", "-i", unsorted_bed], nonredundant_bed)

    labeled_bed = out_dir+"/TIDAL_out/"+sample_name+"_result/"+sample_name+"_TIDAL_nonredundant.bed"
    with open(nonredundant_bed, "r") as inbed, open(labeled_bed, "a") as outbed:
        outbed.write('track name="'+sample_name+'_tidal" description="'+sample_name+'_tidal"\n')
        for x,line in enumerate(inbed):
            line = line.replace("|reference|","|nonabs|")
            split_line = line.split("\t")
            split_line[3] += str(x)
            outbed.write("\t".join(split_line))
    

    run_command(["rm", tmp_bed, ref_insert_bed, depletion_bed, nonabs_bed, unsorted_bed, nonredundant_bed])





if __name__ == "__main__":                
    main()