#!/usr/bin/env python
'''
This is a modified version of the code present in:
https://github.com/CGATOxford/UMI-tools/blob/77e186c6b51917136fe5b4faf3f01ead7eb75aff/umi_tools/group.py
'''
import sys
from collections import Counter
import os

# required to make iteritems python2 and python3 compatible
from builtins import dict
from future.utils import iteritems

import pysam

import umi_tools.Utilities as U
import umi_tools.Documentation as Documentation
import umi_tools.network as network
import umi_tools.sam_methods as sam_methods

# add the generic docstring text
__doc__ = __doc__ + Documentation.GENERIC_DOCSTRING_GDC
__doc__ = __doc__ + Documentation.GROUP_DEDUP_GENERIC_OPTIONS

usage = '''
group - Group reads based on their UMI

Usage: umi_tools group --output-bam [OPTIONS] [--stdin=INFILE.bam] [--stdout=OUTFILE.bam]

       note: If --stdout is ommited, standard out is output. To
             generate a valid BAM file on standard out, please
             redirect log with --log=LOGFILE or --log2stderr '''

def parse_bed(bed_regions):
    with open(bed_regions) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                print("Ignoring BED entry: {}".format(line))
                continue

            region = {
                "chr": cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
                "name": cols[3],
            }
            yield region

def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = U.OptionParser(version="%prog version: $Id$",
                            usage=usage,
                            description=globals()["__doc__"])

    if len(argv) == 1:
        parser.print_usage()
        print ("Required options missing, see --help for more details")
        return 1

    group = U.OptionGroup(parser, "group-specific options")

    group.add_option("--keep", dest="keep", type="string",
                    help="Keep sequences that contain Ds and Ns",
                    default=None)

    group.add_option("--bed", dest="bed", type="string",
                    help="Bed file",
                    default=None)
    
    group.add_option("--sizethreshold", dest="thres", type="int",
                    help="UMI family size threshold",
                    default=None)
    
    group.add_option("--size-out", dest="tsv", type='string',
                     help="Outfile name for file containing distribution of umi family sizes",
                     default=None)
    

    group.add_option("--output-bam", dest="output_bam", action="store_true",
                     default=False,
                     help=("output a bam file with read groups tagged using the UG tag"
                           "[default=%default]"))

    parser.add_option("--umi-group-tag", dest="umi_group_tag",
                      type="string", help="tag for the outputted umi group",
                      default='BX')

    parser.add_option_group(group)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = U.Start(parser, argv=argv, add_group_sam_options=True)

    U.validateSamOptions(options, group=True)

    if options.stdin != sys.stdin:
        in_name = options.stdin.name
        options.stdin.close()
    else:
        raise ValueError("Input on standard in not currently supported")

    if not options.no_sort_output:  # need to determine the output format for sort
        if options.out_sam:
            sort_format = "sam"
        else:
            sort_format = "bam"

    if options.in_sam:
        in_mode = "r"
    else:
        in_mode = "rb"

    if options.out_sam:
        out_mode = "wh"
    else:
        out_mode = "wb"

    infile = pysam.Samfile(in_name, in_mode)

    gene_tag = options.gene_tag
    metacontig2contig = None

    if options.unmapped_reads in ["use", "output"]:
        output_unmapped = True
    else:
        output_unmapped = False

    if options.chrom:
        inreads = infile.fetch(reference=options.chrom)
    else:
        if options.per_gene and options.gene_transcript_map:
            metacontig2contig = sam_methods.getMetaContig2contig(
                infile, options.gene_transcript_map)
            metatag = "MC"
            inreads = sam_methods.metafetcher(infile, metacontig2contig, metatag)
            gene_tag = metatag

        else:
            inreads = infile.fetch(until_eof=output_unmapped)

    bundle_iterator = sam_methods.get_bundles(
        options,
        all_reads=True,
        return_read2=True,
        return_unmapped=output_unmapped,
        metacontig_contig=metacontig2contig)

    # set up UMIClusterer functor with methods specific to
    # specified options.method
    processor = network.UMIClusterer(options.method)

    if options.tsv:
        mapping_outfile = U.openFile(options.tsv, "w")
        
    nInput, nOutput, unique_id, input_reads, output_reads = 0, 0, 0, 0, 0

    for region in parse_bed(options.bed):


        if options.no_sort_output:
            out_name = region['name']
        else:
            out_name = U.getTempFilename(dir=options.tmpdir)
            sorted_out_name = region['name'] + '.bam'
        options.stdout.close()
        assert options.output_bam, (
            "To output a bam you must include --output-bam option")

        if options.output_bam:
            outfile = pysam.Samfile(out_name, out_mode, template=infile)
        else:
            outfile = None

        inreads = infile.fetch(contig=region["chr"], start=region["start"], stop=region["end"])

        groupcounter = Counter()

        for bundle, key, status in bundle_iterator(inreads):

            # write out read2s and unmapped/chimeric (if these options are set)
            if status == 'single_read':
                # bundle is just a single read here
                nInput += 1

                if outfile:
                    outfile.write(bundle)

                nOutput += 1
                continue

            umis = bundle.keys()
            counts = {umi: bundle[umi]["count"] for umi in umis}

            nInput += sum(counts.values())

            while nOutput >= output_reads + 10000:
                output_reads += 10000
                U.info("Written out %i reads" % output_reads)

            while nInput >= input_reads + 1000000:
                input_reads += 1000000
                U.info("Parsed %i input reads" % input_reads)

            # group the umis
            groups = processor(
                counts,
                threshold=options.threshold)
            
            for umi_group in groups:
                top_umi = umi_group[0]

                group_count = sum(counts[umi] for umi in umi_group)

                groupcounter[group_count] +=1
                
                if group_count < options.thres:
                    continue

                for umi in umi_group:
                    reads = bundle[umi]['read']
                    for read in reads:
                        if outfile:
                            # Add the 'UG' tag to the read
                            read.set_tag('UG', unique_id)
                            read.set_tag(options.umi_group_tag, top_umi)
                            outfile.write(read)

                        nOutput += 1

                unique_id += 1

        if outfile:
            outfile.close()
            if not options.no_sort_output:
                # sort the output
                pysam.sort("-o", sorted_out_name, "-O", sort_format, "--no-PG", out_name)
                os.unlink(out_name)  # delete the tempfile
                pysam.index(sorted_out_name)

        if options.tsv:
            for i in sorted(groupcounter.items()):
                mapping_outfile.write("%s\n" % "\t".join(map(str, (
                    i[0],
                    i[1],                
                    region['name']))))                
            mapping_outfile.close()

        # write footer and output benchmark information.
        U.info(
            "Reads: %s" % ", ".join(["%s: %s" % (x[0], x[1]) for x in
                                    bundle_iterator.read_events.most_common()]))
        U.info("Number of reads out: %i, Number of groups: %i" %
            (nOutput, unique_id))

        U.info("Total number of positions deduplicated: %i" %
            processor.positions)
        if processor.positions > 0:
            U.info("Mean number of unique UMIs per position: %.2f" %
                (float(processor.total_umis_per_position) /
                    processor.positions))
            U.info("Max. number of unique UMIs per position: %i" %
                processor.max_umis_per_position)
        else:
            U.warn("The BAM did not contain any valid "
                "reads/read pairs for deduplication")

        U.Stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
