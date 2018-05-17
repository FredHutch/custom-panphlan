#!/usr/bin/env python3
"""Download reference data from PATRIC."""

import os
import sys
import time
import gzip
import json
import logging
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from urllib.request import urlopen
from urllib.error import URLError
from collections import defaultdict


def get_patric_data(output_folder=None, test=False, MAX_RETRIES=10):
    # Make sure the output folder exists
    if not os.path.exists(output_folder):
        logging.info(
            "Making directory for PATRIC data ({})".format(output_folder))
        os.mkdir(output_folder)

    # Get the genome metadata
    genome_metadata_fp = os.path.join(output_folder, "genome_metadata.tsv")
    if os.path.exists(genome_metadata_fp):
        logging.info("Genome metadata file exists, skipping ({})".format(genome_metadata_fp))
    else:
        logging.info("Writing PATRIC genome metadata to " + genome_metadata_fp)

        with open(genome_metadata_fp, "wt") as fo:
            with urlopen("ftp://ftp.patricbrc.org/RELEASE_NOTES/genome_metadata") as fi:
                fo.write(fi.read().decode("utf-8"))
        logging.info("Done writing PATRIC genome metadata to " + genome_metadata_fp)

    # Set up the folders for genomes, transcripts, and pathways
    output_subfolders = {
        data_type: os.path.join(output_folder, data_type)
        for data_type in ["transcripts", "pathways"]
    }
    for f in output_subfolders.values():
        if not os.path.exists(f):
            logging.info("Making subfolder " + f)
            os.mkdir(f)

    # Read in all of the genome metadata
    logging.info("Downloading data for individual genomes")

    for row_ix, genome_entry in enumerate(read_tsv(genome_metadata_fp)):
        for retry_ix in range(MAX_RETRIES):
            try:
                get_genome_entry(genome_entry["genome_id"], output_subfolders)
                break
            except:
                logging.info("Problem getting data for {}, waiting 10 seconds to retry".format(
                    genome_entry["genome_id"]
                ))
                time.sleep(10)

        # If we're running in testing mode, stop after 10 genomes
        if test and row_ix == 10:
            logging.info("TESTING MODE: stopping after 10 genomes")
            break

    # Now make the input files needed for PICRUSt
    make_picrust_inputs(
        output_subfolders["transcripts"],
        output_subfolders["pathways"],
        output_folder
    )

def read_tsv(fp):
    with open(fp, "rt") as f:
        # Get the header line
        header = f.readline()
        # Parse as tab-delimited
        header = header.rstrip("\n").split("\t")
        # Get all of the lines for each individual genome
        for line in f:
            # Add the headers to each field
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def make_picrust_inputs(transcript_folder, pathway_folder, output_folder):
    """
    
    Make all of the files needed by PICRUSt to make a custom database
    
    Output files:
        * Marker gene sequences (FASTA .fasta)
        * Marker gene copy number (tab-delimited .tab)
        * Functional trait copy number (tab-delimited .tab)
        * Mapping file for tree tip ids to genome ids (tab-delimited .tab)
    """

    marker_genes = {}
    marker_gene_copy_number = {}
    functional_trait_copy_number = {}
    genome_mapping = {}

    transcript_files = {
        genome_id: [
            os.path.join(transcript_folder, genome_id, f)
            for f in os.listdir(os.path.join(transcript_folder, genome_id))
            if f.endswith(".frn")
        ][0]
        for genome_id in os.listdir(transcript_folder)
    }

    pathway_files = {
        genome_id: [
            os.path.join(pathway_folder, genome_id, f)
            for f in os.listdir(os.path.join(pathway_folder, genome_id))
            if f.endswith(".pathway.tab")
        ][0]
        for genome_id in os.listdir(pathway_folder)
    }

    for genome_id in transcript_files:
        if genome_id not in pathway_files:
            logging.info("Transcripts but no pathway information for {}, skipping".format(genome_id))
            continue

        logging.info("Processing " + genome_id)

        logging.info("Reading in transcripts for " + genome_id)

        ssu_rrna = {
            header.split(" ")[0].replace("|", "_").rstrip("_"): seq
            for header, seq in SimpleFastaParser(open(transcript_files[genome_id], "rt"))
            if " 16S " in header or " SSU " in header
        }
        logging.info("Found {} 16S / SSU rRNA for {}".format(
            len(ssu_rrna), genome_id
        ))
        if len(ssu_rrna) == 0:
            logging.info("Skipping " + genome_id)
            continue

        logging.info("Reading in pathway information for " + genome_id)
        pathways = defaultdict(float)
        for entry in read_tsv(pathway_files[genome_id]):
            pathway_name = "{}: {}".format(entry["pathway_id"], entry["pathway_name"])
            pathways[pathway_name] += 1

        logging.info("Read in {:,} pathway annotations for {}".format(
            len(pathways), genome_id
        ))
        if len(pathways) == 0:
            logging.info("Skipping " + genome_id)
            continue

        # At this point we're happy that this record (1) has 16S and (2) has pathways

        # Record the marker genes
        for header, seq in ssu_rrna.items():
            assert header not in marker_genes
            marker_genes[header] = seq

            # Map the 16S headers to the genome
            assert header not in genome_mapping
            genome_mapping[header] = genome_id

        # Set the marker gene copy number
        marker_gene_copy_number[genome_id] = len(ssu_rrna)

        # Set the functional trait copy number
        functional_trait_copy_number[genome_id] = pathways

    logging.info("Writing out information for all genomes")

    logging.info("Writing out 16S copy number information")
    with open(os.path.join(output_folder, "16S_copy_number.tab"), "wt") as fo:
        fo.write("OTU_IDs\t16S_rRNA_Count\n")
        for header, copy_number in marker_gene_copy_number.items():
            fo.write("{}\t{}\n".format(header, float(copy_number)))

    logging.info("Writing out 16S sequences")
    with open(os.path.join(output_folder, "16S_sequences.fasta"), "wt") as fo:
        for header, seq in marker_genes.items():
            fo.write(">{}\n{}\n".format(header, seq))

    logging.info("Writing out pathway abundances")
    all_pathway_ids = list(set([
        pathway_id
        for genome_id in functional_trait_copy_number
        for pathway_id in functional_trait_copy_number[genome_id]
    ]))
    with open(os.path.join(output_folder, "pathway_abundance.tab"), "wt") as fo:
        fo.write("OTU_IDs\t{}\n".format("\t".join(all_pathway_ids)))
        for genome_id in functional_trait_copy_number:
            fo.write("{}\t{}\n".format(
                genome_id,
                "\t".join(map(str, [
                    functional_trait_copy_number[genome_id][pathway_id]
                    for pathway_id in all_pathway_ids
                ]))
            ))

    logging.info("Writing out the 16S -> genome mapping")
    with open(os.path.join(output_folder, "genome_mapping.tab"), "wt") as fo:
        fo.write("#OTU_IDs\tgenome_ID\n")
        for otu_id, genome_id in genome_mapping.items():
            fo.write("{}\t{}\n".format(otu_id, genome_id))

    logging.info("Done")


def get_genome_entry(
    genome_id,
    output_subfolders,
    suffix={
        "transcripts": [".PATRIC.frn", ".RefSeq.frn"],
        "pathways": [".PATRIC.pathway.tab", ".RefSeq.pathway.tab"],
}):
    """Get information for a single genome from PATRIC."""

    # Make a subfolder for each of the types of data to retrieve
    for data_type, data_suffixes in suffix.items():
        assert data_type in output_subfolders
        genome_subfolder = os.path.join(output_subfolders[data_type], genome_id)
        if not os.path.exists(genome_subfolder):
            os.mkdir(genome_subfolder)

        # Keep track of whether this piece of data was able to be downloaded
        any_written = False

        # Try all possible suffixes
        for data_suffix in data_suffixes:
            # Download the data for this type of data for this genome
            data_url = "ftp://ftp.patricbrc.org/genomes/{}/{}{}".format(
                genome_id, genome_id, data_suffix
            )
            local_path = os.path.join(
                genome_subfolder, genome_id + data_suffix
            )
            if os.path.exists(local_path):
                any_written = True
                logging.info("Found local copy of " + local_path)
                break

            logging.info("Downloading " + data_url)
            try:
                with open(local_path, "wt") as fo:
                    with urlopen(data_url) as fi:
                        fo.write(fi.read().decode("utf-8"))
                any_written = True
                break
            except URLError:
                os.remove(local_path)
                logging.info("URL not found, trying next option")

        assert any_written, "Was not able to download {} for {}".format(data_type, genome_id)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Download reference data from PATRIC"""
    )

    parser.add_argument("--output-folder",
                        type=str,
                        required=True,
                        help="""Folder for downloaded genome data.""")
    parser.add_argument("--test",
                        action="store_true",
                        help="""Use this flag to download a subset of the data for testing.""")
    
    args = parser.parse_args(sys.argv[1:])

    # Set up logging
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [get-patric-data] %(message)s'
    )
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write logs to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    get_patric_data(
        **args.__dict__
    )
