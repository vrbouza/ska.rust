//! Split k-mer analysis (version 2) uses exact matching of split k-mer sequences to align closely related
//! sequences, typically small haploid genomes such as bacteria and viruses.
//!
//! SKA can only align SNPs further than the k-mer length apart,
//! and does not use a gap penalty approach or give alignment scores.
//! But the advantages are speed and flexibility, particularly the ability to
//! run on a reference-free manner (i.e. including accessory genome variation)
//! on both assemblies and reads.
//!
//! ## Details
//!
//! A split k-mer is a sequence of bases with the middle position removed. For
//! example, the split 9-mer pattern is `XXXX-XXXX`, and the split 9-mers of the
//! sequence `ACGAGAGTCTT` are:
//!
//! | Split k-mer | Middle base |
//! |-------------|-------------|
//! | `ACGAAGTC`  | `G`         |
//! | `CGAGGTCT`  | `A`         |
//! | `GAGATCTT`  | `G`         |
//!
//! Which is broadly how SKA represents sequences in `.skf` files, as a dictionary with split
//! k-mers as keys and the middle base as values. These dictionaries are merged
//! by exactly matching the split k-mers, such that the middle positions are aligned (but
//! unordered). Split k-mers can also be matched to those from a reference sequence
//! to give an ordered pseudoalignment.
//!
//! Various optimisations are used to make this as fast as possible. For a more thorough comparison with version 1.0 of SKA, see the
//! [github description](https://github.com/bacpop/ska.rust/blob/master/README.md).
//!
//! *NB As split k-mers are even lengths, it is possible that they are their
//! own reverse complement. The original version of ska would randomly pick a strand,
//! possibly introducing a SNP across samples. This version uses an ambiguous middle
//! base (W for A/T; S for C/G) to represent this case.*
//!
//! Tutorials:
//! - [From genomes to trees](https://www.bacpop.org/guides/building_trees_with_ska/).
//! - [Filtering options](https://www.bacpop.org/guides/snp_alignment_with_ska/).
//!
//! Papers:
//! - [Seamless, rapid, and accurate analyses of outbreak genomic data using split k-mer analysis ](https://genome.cshlp.org/content/34/10/1661.abstract).
//! - [skalo: using SKA split k-mers with coloured de Brujin graphs to genotype indels](https://academic.oup.com/mbe/article/42/4/msaf077/8103706).
//!
//! Command line usage follows. For API documentation and usage, see the [end of this section](#api-usage).
//!
//! # Usage
//!
//! `.skf` files represent merged split k-mers from multiple sequence files. They
//! are created with `ska build`. You can subsequently use `ska align` to write
//! out an unordered alignment from these files, or `ska map` to write an alignment
//! ordered versus a reference sequence.
//!
//! Alternatively, both `ska align` and `ska map` can take input sequence directly to obtain output alignments
//! in a single command and without saving an `.skf` file. NB: This uses the default options
//! of `ska build`, so to change these you will need to run the alignment in two steps.
//!
//! Output from `ska align` and `ska map` is to STDOUT, so you can use a redirect `>` to save to a file or pipe `|`
//! to stream into another compatible program on STDIN. You can also add an output
//! file prefix directly with `-o` (for `ska build` this is required).
//!
//! ## Important notes
//!
//! - In this version of ska the k-mer length is the _total_ split k-mer size,
//!   whereas in ska1 the k-mer length is half the split length. So the default of
//!   k=15 in ska1 corresponds to choosing `-k 31` in this version.
//! - If you are using FASTQ input files (reads), these must be provided as two
//!   deinterleaved files using the `-f` option to `ska build`. Providing them as
//!   a single file will treat them as FASTA, and not correct sequencing errors.
//! - If you are running `ska map`, it is more memory efficient to run these one at
//!   a time, rather than merging a single `.skf` file. A command for doing this
//!   in parallel is listed below.
//!
//! ## Common options
//!
//! Version can be viewed by running `ska -V`.
//!
//! Details and progress messages are written on STDERR. You can see more logging
//! information by adding the verbose flag `-v`.
//!
//! ## ska build
//!
//! This command creates an `.skf` file from sequence (.fasta/.fasta.gz/.fastq/.fastq.gz) input.
//! K-mer size must be odd, greater than 5, and a maximum of 63. Smaller k-mers
//! are more sensitive and can find closer positions, but are less specific so
//! may lead to more repeated split k-mers and ambiguous bases. Using k <= 31 uses
//! 64-bit integers and may be faster than 31 < k <= 63, which uses 128-bits.
//!
//! This is typically the most computationally intensive step of `ska`, and
//! multiple `--threads` can be used to split the work over multiple CPU cores.
//!
//! Build from two input FASTA files with a k-mer size of 31:
//! ```bash
//! ska build -o seqs -k 31 assemblies/seq1.fa assemblies/seq2.fa
//! ```
//! This will assume sample names of `seq1` and `seq2` –- the base path and with known fastx
//! file extensions removed. If you know the strand,
//! for example when exclusively using reference sequences or single stranded viruses, add `--single-strand`
//! to ignore reverse complements.
//!
//! ### Read files
//!
//! To use FASTQ files, or specify sample names, or more easily input a larger number of input files,
//! you can provide a tab separated file list via `-f` instead of listing files. For example
//! a file called `input_sequence.txt` to describe `sample_1` (an assembly) and `sample_2` (paired reads):
//! ```text
//! sample_1    assemblies/seq1.fa
//! sample_2    reads/seq2_1.fastq.gz   reads/seq2_2.fastq.gz
//! ```
//! You'd run:
//! ```bash
//! ska build -o seqs -f input_sequence.txt
//! ```
//! Options for filtering/error correction are:
//! - `--min-count`. Specify a minimum number of appearances a k-mer must have
//!   to be included. This is an effective way of filtering sequencing errors if set
//!   to at least three, but higher may be appropriate for higher coverage data.
//!   A two-step blocked bloom and hash table filter is used for memory efficiency.
//! - `--proportion-reads`. Specify a proportion of reads to use to build the .skf file.
//!   The value of this parameter must be between 0 and 1. If you have high coverage samples
//!   this can be used to reduce the build time.
//! - `--qual-filter`. `none` do not filter based on quality scores.
//!   `middle` (default) filter k-mers where the middle base is below the minimum quality.
//!   `strict` filter k-mers where any base is below the minimum quality.
//! - `--min-qual`. Specify a minimum PHRED score to use in the filter.
//!
//! FASTQ files must be paired end. If you'd like to request more flexibility in
//! this regard please contact us.
//!
//! ### Threads
//!
//! The maximum threads actually used will be a power of two, so if you provided
//! `--threads 6` only four would be used. Additionally, at least ten samples per
//! thread are required so maximums are:
//!
//! | Samples | Maximum threads |
//! |---------|-----------------|
//! | 1-19    | 1               |
//! | 20-39   | 2               |
//! | 40-79   | 4               |
//!
//! and so on. Use `-v` to see a message with the number being used.
//!
//! Using more threads will increase memory use.
//!
//! You can also run blocks of samples independently (e.g. with snakemake or a
//! job array) then use `ska merge` to combine results.
//!
//! ## ska align
//!
//! Creates an alignment from a `.skf` file or sequence files. Sites (columns) are
//! in an arbitrary order. Two basic filters are available: `--min-freq` which
//! sets the maximum number of missing sites; `--filter` which can be set to
//! one of three options:
//! * `no-filter` -- no extra filtering
//! * `no-const` -- no constant sites
//! * `no-ambig-or-const` -- no constant sites, or sites where the only variable base is ambiguous
//!
//! `no-filter` may be useful for ascertainment bias correction in phylogenetic algorithms,
//! but note the flanking split k-mers will never be included. `no-const` is the default.
//! `no-ambig-or-const` can be used when you want to treat any ambiguity as an `N`.
//!
//! With an `.skf` file from `ska build`, constant sites, and no missing variants:
//! ```bash
//! ska align --min-freq 1 --filter no-filter -o seqs seqs.skf
//! ```
//!
//! Another example: directly from FASTA files, with default `build` and `align` settings,
//! and gzipping the output alignment
//! ```bash
//! ska align seq*.fa --threads 8 | gzip -c - > seqs.aln.gz
//! ```
//!
//! ## ska map
//!
//! Creates an alignment from a `.skf` file or sequence files by mapping its
//! split k-mers to split k-mers of a linear reference sequence. This produces pseudoalignment
//! with the same sites/columns as the reference sequence. Sites not mapped will
//! be output as missing '-'.
//!
//! Pass the FASTA file of the reference as the first argument, and the skf file as
//! the second argument:
//! ```bash
//! ska map ref.fa seqs.skf -o ref_mapped.aln
//! ```
//! Add `--repeat-mask` to mask any repeated split k-mers in the reference with 'N'.
//!
//! You can also get a VCF as output, which has rows as variants, and only has the
//! variable sites (but will include unmapped bases as missing).
//! An example command, also demonstrating that everything can be done from input sequences in a single command:
//! ```bash
//! ska map ref.fa seq1.fa seq2.fa -f vcf --threads 2 | bgzip -c - > seqs.vcf.gz
//! ```
//!
//! ## ska lo
//!
//! Converts split k-mers from a `.skf` file into a colored De Bruijn graph and infers indels from graph bubbles and SNPs from variant groups in
//! reference-free mode (as with `ska align`). SNPs are only composed of ATGC variants (no ambigous nucleotides). The same
//! filtering applies to indels.
//! Multithreading ('-t' argument) is not fully optimised - it usually takes 4 threads to halve runtimes.
//!
//! To generate a SNP alignment and an indel VCF file (here named 'test_snps.fas' and 'test_indels.vcf'):
//! ```bash
//! ska lo seqs.skf test
//! ```
//!
//! It can also position SNPs on a reference genome if provided using the '-r' argument. The reference genome should
//! be in FASTA format and composed of a unique sequence. skalo lo will then generate, in addition to the SNP alignment and the indel
//! VCF file, a SNP VCF file and a pseudo-genome alignment (as with `ska map`) that can be used for recombination analyses.
//! In such use case, we recommmend to increase the maximum proportion of missing data allowed per variant ('-m' argument), but not above 0.5.
//! ```bash
//! ska lo seqs.skf test -r reference.fas -m 0.4
//! ```
//! Please note that at the moment indels cannot be positioned on a reference genome.
//!
//!
//! ### Efficiency
//!
//! As `ska map` is independent for each input file, the most efficient way to
//! run this on a number of samples is with gnu parallel, for example:
//! ```bash
//! cat ${FILE_LIST} | parallel -j 8 "ID=\$(basename {} | sed 's/_a.*fasta/_temp_skf/g');
//! ska build -o \$ID -k 31 {} ;
//! MAPOUT=\$(basename {} | sed 's/_a.*fasta/_aln/g');
//!  ska map $REFERENCE_LOC \${ID}.skf -o \${MAPOUT}.aln"
//! cat *.aln > ska_map.aln
//! ```
//!
//! ## ska distance
//!
//! Use to calculate distances between all samples within an `.skf` file. The output
//! will contain the number of SNP differences between all pairs of samples, as
//! well as the proportion of matching k-mers.
//!
//! ```bash
//! ska distance -o distances.txt seqs.skf
//! ```
//!
//! Mismatches is an estimate of the alignable genome fraction based on the unfiltered
//! k-mers.
//!
//! Also consider ambiguous bases by adding `--allow-ambiguous` flag, which will give partial SNP distances.
//! Note that ambiguous bases may overestimate distances due to repeat k-mers. Use `--min-freq` to
//! ignore k-mers only found in some samples at the population level (default = 0.0). F
//! or finer control over filtering, first run `ska weed`
//! on the input .skf.
//!
//! Multiple threads can be used, but this will only be effective with large numbers of samples.
//!
//! The companion script in `scripts/cluster_dists.py` (requires `networkx`) can
//! be used to make single linkage clusters from these distances at given thresholds,
//! and create a visualisation in [Microreact](https://microreact.org/):
//! ```bash
//! ska distance seqs.skf > distances.txt
//! python scripts/cluster_dists.py distances.txt --snps 20 --mismatches 0.05
//! ```
//! If you install `rapidnj` and `biopython` you can also draw an NJ tree from
//! these distances, which will be displayed in Microreact. Use your `--api-key`
//! to directly upload and get the URL printed on the terminal.
//!
//! ## ska merge
//!
//! Use to combine multiple `.skf` files into one, typically for subsequent use in `align` or `map`.
//! This may be particularly useful if `ska build` was run on multiple input files
//! one at a time (for example in a job array).
//!
//! ```bash
//! ska merge -o all_samples sample1.skf sample2.skf sample3.skf
//! ```
//!
//! ## ska delete
//!
//! Remove samples by name from an `.skf` file. All sample names must exist in the
//! file, or an error will be thrown. After removing samples, the input `.skf` file will be overwritten,
//! unless you provide a new one with `-o`.
//! Any split k-mers which have no observed missing bases after removal will also be deleted.
//!
//! ```bash
//! ska delete --skf-file all_samples.skf sample_1 sample_3
//! ```
//! If you wish to delete many samples, you can use `-f` to provide their names
//! by line via an input file.
//!
//! ## ska weed
//!
//! Remove (weed) split k-mers from an `.skf` file. The split k-mers to be removed
//! are generated from a FASTA file, which may for example contain known transposons or
//! other repeat sequences. As with `delete`, the input `.skf` file is overwritten by default.
//! Use `-o` to write the results to a new file instead of overwriting the input `.skf`.
//!
//! ```bash
//! ska weed all_samples.skf MGEs.fa
//! ```
//!
//! In addition, you can also use `ska weed` to filter split k-mers by proportion of samples
//! they appear in, and constant or amibguous middle cases, which does
//! not require a list of split k-mers to be removed (but both can be used together, if you wish).
//!
//! For example, you may want to reduce the size of an `.skf` file for online use.
//! You can do this by removing any constant sites or ambiguous-only s0tes (which are typically unused), and by hard-filtering
//! by frequency (i.e. before writing output):
//! ```bash
//! ska weed --filter no-ambig-or-const --min-freq 0.9 all_samples.skf
//! ```
//!
//! ## ska nk
//! Return information on the *n*umber of *k*-mers in an `.skf` file. This will
//! print on STDOUT:
//! - The version of ska used to build the file.
//! - The k-mer size.
//! - Number of bits used for the split k-mer (64 or 128).
//! - Whether reverse complements were used.
//! - The total number of split k-mers.
//! - The total number of samples.
//! - The sample names.
//! - The number of split k-mers found in each sample.
//!
//! ```bash
//! ska nk all_samples.skf
//! ska_version=0.3.1
//! k=21
//! k_bits=64
//! rc=true
//! k-mers=3228084
//! samples=28
//! sample_names=["19183_4#73", "12673_8#29", "12673_8#31", "12754_5#61", "12754_5#89", "12754_5#47", "12754_5#32", "12754_5#78", "12754_4#85", "12754_5#74", "19183_4#57", "12754_5#36", "19183_4#49", "19183_4#79", "19183_4#60", "12754_5#24", "12754_5#22", "12754_5#71", "12673_8#26", "12754_5#95", "12754_5#86", "12673_8#24", "19183_4#61", "12673_8#41", "12754_4#66", "12754_5#80", "12754_5#84", "19183_4#78"]
//! sample_kmers=[2872587, 2997448, 2949719, 2997496, 2997178, 2912749, 2996491, 2997221, 2949102, 2997454, 2914109, 2912237, 2872518, 2869957, 2872470, 2997992, 2997647, 2958512, 2998099, 2997290, 2950253, 3027707, 2997881, 2907920, 2911447, 2997644, 2944830, 2915080]
//! ```
//!
//! If you add `--full-info`, the split k-mer dictionary will be decoded and printed
//! to STDOUT after the baseline information, for example:
//! ```bash
//! TAAATATC        TAAACGCC        A,-
//! AGACTCTC        TACAGCTA        G,G
//! AAACCATC        AAACACTC        A,-
//! TTAAAAGA        TCTCGTAC        C,C
//! ```
//! This is tab-separated, showing the first and second half of the split k-mer
//! and the middle bases of each sample (comma seperated).
//!
//! NB: Only one split k-mer is shown even if the reverse complement was used.
//! These are not precisely canonical k-mers, as the encoding order `{A, C, T, G}` is used internally.
//! But if you can't find a sequence in your input file, you will find its reverse complement.
//!
//! ## ska cov
//!
//! Estimate a coverage cutoff for use with read data. This will count the split
//! k-mers in a pair of FASTQ samples, and create a histogram of these counts.
//! A mixture model is then fitted to this histogram using maximum likelihood,
//! which can give a suggested cutoff with noisy data.
//!
//! The cutoff will be printed to STDERR. Use `-v` to get additional information on the
//! optimisation progress and result. A table of counts and the fit will be printed
//! to STDOUT, which can then be plotted by the companion script in
//! `scripts/plot_cov.py` (requires `matplotlib`):
//! ```bash
//! ska cov reads_1.fastq.gz reads_2.fastq.gz -k 31 -v > cov_plot.txt
//! python scripts/plot_cov.py cov_plot.txt
//! ```
//!
//! The cutoff can be used with the `--min-count` parameter of `ska build`. For
//! a set of sequence experiments with similar characteristics it may be sufficient
//! to use the same cutoff. Alternatively `ska cov` can be run on every sample
//! independently (`gnu parallel` would be an efficient way to do this).
//!
//! ## ska lo
//!
//! Runs "skalo", a graph-based algorithm for inferring variants from WGS outbreak data.
//! Skalo takes as input a `.skf` file generated by the `ska build` command.
//!
//! To run skalo on a file named `my_file.skf` use the following command `ska lo my_file.skf skalo_results`.
//!
//! The skalo method
//! can be modified using the following flags:
//! * `-r` the path to a reference genome that the detected SNPs can be mapped onto
//! * `-m` the maximum fraction of missing data permissible for a site where a possible SNP has been detected (default = 0.2)
//! * `-d` the depth of the recursive path allowed during graph traversal (default = 4)
//! * `-n` maximum number of internal indel k-mers (default = 2)
//! * `-t` number of threads to use for parallel computation (default = 1)
//!
//!
//! # API usage
//!
//! See the submodule documentation linked below.
//!
//! To use `ska build` in other rust code:
//! ```rust
//! use ska::merge_ska_dict::{InputFastx, build_and_merge};
//! use ska::merge_ska_array::MergeSkaArray;
//! use ska::{QualOpts, QualFilter};
//!
//! // Build, merge
//! let input_files: [InputFastx; 2] = [("test1".to_string(),
//!                                      "tests/test_files_in/test_1.fa".to_string(),
//!                                      None),
//!                                     ("test2".to_string(),
//!                                      "tests/test_files_in/test_2.fa".to_string(),
//!                                      None)];
//! let rc = true;
//! let k = 17;
//! let quality = QualOpts {min_count: 1, min_qual: 0, qual_filter: QualFilter::NoFilter};
//! let threads = 2;
//! let proportion_reads = None;
//! // NB u64 for k<=31, u128 for k<=63
//! let merged_dict =
//!     build_and_merge::<u64>(&input_files, k, rc, &quality, threads, proportion_reads);
//!
//! // Save
//! let ska_array = MergeSkaArray::new(&merged_dict);
//! ska_array
//!     .save(format!("merged.skf").as_str())
//!     .expect("Failed to save output file");
//! ```
//!
//! To use `ska align` in other rust code:
//! ```rust
//! use ska::io_utils::{load_array, set_ostream};
//! use ska::cli::FilterType;
//!
//! // Load an .skf file
//! let threads = 4;
//! let input = vec!["tests/test_files_in/merge.skf".to_string()];
//! let mut ska_array = load_array::<u64>(&input, threads).expect("Could not open input as u64");
//!
//! // Apply filters
//! let min_count = 2;
//! let filter_ambig_as_missing = false;
//! let update_kmers = false;
//! let filter = FilterType::NoConst;
//! let ignore_const_gaps = false;
//! let ambig_mask = false;
//! ska_array.filter(min_count, filter_ambig_as_missing, &filter, ambig_mask, ignore_const_gaps, update_kmers);
//!
//! // Write out to stdout
//! let mut out_stream = set_ostream(&None);
//! ska_array
//!     .write_fasta(&mut out_stream)
//!     .expect("Couldn't write output fasta");
//! ```
//!
//! To use `ska map` in other rust code, see the example for [`ska_ref`].
//!
//! ## Other APIs
//!
//! It would be possible to make a python API using [`maturin`](https://github.com/PyO3/maturin).
//! Please submit a feature request if you would find this useful.
//!

#![warn(missing_docs)]
use std::fmt;
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

use clap::ValueEnum;
extern crate num_cpus;

#[cfg(not(target_arch = "wasm32"))]
pub mod merge_ska_dict;
pub mod ska_dict;
#[cfg(not(target_arch = "wasm32"))]
use crate::merge_ska_dict::build_and_merge;

pub mod ska_ref;
use crate::ska_ref::RefSka;
#[cfg(not(target_arch = "wasm32"))]
pub mod merge_ska_array;
#[cfg(not(target_arch = "wasm32"))]
use crate::merge_ska_array::MergeSkaArray;

#[cfg(not(target_arch = "wasm32"))]
pub mod generic_modes;
#[cfg(not(target_arch = "wasm32"))]
use crate::generic_modes::*;

#[cfg(not(target_arch = "wasm32"))]
pub mod cli;
#[cfg(not(target_arch = "wasm32"))]
use crate::cli::*;

#[cfg(not(target_arch = "wasm32"))]
pub mod io_utils;
#[cfg(not(target_arch = "wasm32"))]
use crate::io_utils::*;

#[cfg(not(target_arch = "wasm32"))]
pub mod coverage;
#[cfg(not(target_arch = "wasm32"))]
use crate::coverage::CoverageHistogram;

#[cfg(not(target_arch = "wasm32"))]
pub mod skalo;
#[cfg(not(target_arch = "wasm32"))]
use crate::io_utils::load_array;
#[cfg(not(target_arch = "wasm32"))]
use crate::skalo::utils::Config;

#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;
#[cfg(target_arch = "wasm32")]
use wasm_bindgen_file_reader::WebSysFile;
#[cfg(target_arch = "wasm32")]
extern crate console_error_panic_hook;
#[cfg(target_arch = "wasm32")]
pub mod fastx_wasm;
#[cfg(target_arch = "wasm32")]
use json;
#[cfg(target_arch = "wasm32")]
pub mod ska_map;
#[cfg(target_arch = "wasm32")]
use crate::ska_map::SkaMap;
#[cfg(target_arch = "wasm32")]
pub mod ska_align;
#[cfg(target_arch = "wasm32")]
use crate::ska_align::SkaAlign;
#[cfg(target_arch = "wasm32")]
use crate::ska_dict::bit_encoding::UInt;

/// Possible quality score filters when building with reads
#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum)]
pub enum QualFilter {
    /// Ignore quality scores in reads
    NoFilter,
    /// Filter middle bases below quality threshold
    Middle,
    /// Filter entire k-mer when any base below quality threshold
    Strict,
}

impl fmt::Display for QualFilter {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Self::NoFilter => write!(f, "No quality filtering"),
            Self::Middle => write!(f, "Middle base quality filtering"),
            Self::Strict => write!(f, "Whole k-mer quality filtering"),
        }
    }
}

/// Quality filtering options for FASTQ files
pub struct QualOpts {
    /// Minimum k-mer count across reads to be added
    pub min_count: u16,
    /// Minimum base quality to be added
    pub min_qual: u8,
    /// [`QualFilter`]: apply quality across whole k-mer, just middle base, or not at all
    pub qual_filter: QualFilter,
}

impl fmt::Display for QualOpts {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "min count: {}; minimum quality {} ({}); filter applied: {}",
            self.min_count,
            self.min_qual,
            (self.min_qual + 33) as char,
            self.qual_filter
        )
    }
}

#[cfg(not(target_arch = "wasm32"))]
#[doc(hidden)]
pub fn main() {
    let args = cli_args();
    if args.verbose {
        simple_logger::init_with_level(log::Level::Info).unwrap();
    } else {
        simple_logger::init_with_level(log::Level::Warn).unwrap();
    }

    eprintln!("SKA: Split K-mer Analysis (the alignment-free aligner)");
    let start = Instant::now();
    match &args.command {
        Commands::Build {
            seq_files,
            file_list,
            output,
            k,
            single_strand,
            min_count,
            min_qual,
            qual_filter,
            threads,
            proportion_reads,
        } => {
            check_threads(*threads);

            // Read input
            let input_files = get_input_list(file_list, seq_files);
            let rc = !*single_strand;
            // Build, merge
            // check for 64 or 128 bit representation
            let mut quality = QualOpts {
                min_count: 0,
                min_qual: *min_qual,
                qual_filter: *qual_filter,
            };
            if k.le(&31) {
                log::info!("k={}: using 64-bit representation", *k);
                quality.min_count =
                    kmer_min_cutoff::<u64>(min_count, &input_files, k, rc, args.verbose);
                let merged_dict = build_and_merge::<u64>(
                    &input_files,
                    *k,
                    rc,
                    &quality,
                    *threads,
                    *proportion_reads,
                );

                // Save
                save_skf(&merged_dict, output);
            } else {
                log::info!("k={}: using 128-bit representation", *k);
                quality.min_count =
                    kmer_min_cutoff::<u128>(min_count, &input_files, k, rc, args.verbose);
                let merged_dict = build_and_merge::<u128>(
                    &input_files,
                    *k,
                    rc,
                    &quality,
                    *threads,
                    *proportion_reads,
                );

                // Save
                save_skf(&merged_dict, output);
            }
        }
        Commands::Align {
            input,
            output,
            min_freq,
            filter_ambig_as_missing,
            filter,
            ambig_mask,
            no_gap_only_sites,
            threads,
        } => {
            check_threads(*threads);
            if let Ok(mut ska_array) = load_array::<u64>(input, *threads) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                align(
                    &mut ska_array,
                    output,
                    filter,
                    *ambig_mask,
                    *no_gap_only_sites,
                    *min_freq,
                    *filter_ambig_as_missing,
                );
            } else if let Ok(mut ska_array) = load_array::<u128>(input, *threads) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                align(
                    &mut ska_array,
                    output,
                    filter,
                    *ambig_mask,
                    *no_gap_only_sites,
                    *min_freq,
                    *filter_ambig_as_missing,
                );
            } else {
                panic!("Could not read input file(s): {input:?}");
            }
        }
        Commands::Map {
            reference,
            input,
            output,
            format,
            ambig_mask,
            repeat_mask,
            threads,
        } => {
            check_threads(*threads);

            log::info!("Loading skf as dictionary");
            if let Ok(ska_array) = load_array::<u64>(input, *threads) {
                log::info!(
                    "Making skf of reference k={} rc={}",
                    ska_array.kmer_len(),
                    ska_array.rc()
                );
                let mut ska_ref = RefSka::<u64>::new(
                    ska_array.kmer_len(),
                    reference,
                    ska_array.rc(),
                    *ambig_mask,
                    *repeat_mask,
                );
                map(&ska_array, &mut ska_ref, output, format, *threads);
            } else if let Ok(ska_array) = load_array::<u128>(input, *threads) {
                log::info!(
                    "Making skf of reference k={} rc={}",
                    ska_array.kmer_len(),
                    ska_array.rc()
                );
                let mut ska_ref = RefSka::<u128>::new(
                    ska_array.kmer_len(),
                    reference,
                    ska_array.rc(),
                    *ambig_mask,
                    *repeat_mask,
                );
                map(&ska_array, &mut ska_ref, output, format, *threads);
            } else {
                panic!("Could not read input file(s): {input:?}");
            }
        }
        Commands::Distance {
            skf_file,
            output,
            min_freq,
            allow_ambiguous,
            threads,
        } => {
            check_threads(*threads);
            let filter_ambiguous = !*allow_ambiguous;
            if let Ok(mut ska_array) = MergeSkaArray::<u64>::load(skf_file) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                distance(&mut ska_array, output, *min_freq, filter_ambiguous);
            } else if let Ok(mut ska_array) = MergeSkaArray::<u128>::load(skf_file) {
                // In debug mode (cannot be set from CLI, give details)
                log::debug!("{ska_array}");
                distance(&mut ska_array, output, *min_freq, filter_ambiguous);
            } else {
                panic!("Could not read input file(s): {skf_file}");
            }
        }
        Commands::Merge { skf_files, output } => {
            if skf_files.len() < 2 {
                panic!("Need at least two files to merge");
            }

            log::info!("Loading first alignment");
            if let Ok(first_array) = MergeSkaArray::<u64>::load(&skf_files[0]) {
                merge(&first_array, &skf_files[1..], output);
            } else if let Ok(first_array) = MergeSkaArray::<u128>::load(&skf_files[0]) {
                merge(&first_array, &skf_files[1..], output);
            } else {
                panic!("Could not read input file: {}", skf_files[0]);
            }
        }
        Commands::Delete {
            skf_file,
            output,
            file_list,
            names,
        } => {
            let input_files = get_input_list(file_list, names);
            let input_names: Vec<&str> = input_files.iter().map(|t| &*t.0).collect();
            let output_file = output.clone().unwrap_or(skf_file.to_string());
            log::info!("Loading skf file");
            if let Ok(mut ska_array) = MergeSkaArray::<u64>::load(skf_file) {
                delete(&mut ska_array, &input_names, &output_file);
            } else if let Ok(mut ska_array) = MergeSkaArray::<u128>::load(skf_file) {
                delete(&mut ska_array, &input_names, &output_file);
            } else {
                panic!("Could not read input file: {skf_file}");
            }
        }
        Commands::Weed {
            skf_file,
            weed_file,
            output,
            reverse,
            min_freq,
            filter_ambig_as_missing,
            ambig_mask,
            no_gap_only_sites,
            filter,
        } => {
            log::info!("Loading skf file");
            if let Ok(mut ska_array) = MergeSkaArray::<u64>::load(skf_file) {
                weed(
                    &mut ska_array,
                    weed_file,
                    *reverse,
                    *min_freq,
                    *filter_ambig_as_missing,
                    filter,
                    *ambig_mask,
                    *no_gap_only_sites,
                    if output.is_none() {
                        skf_file
                    } else {
                        output.as_ref().unwrap().as_str()
                    },
                );
            } else if let Ok(mut ska_array) = MergeSkaArray::<u128>::load(skf_file) {
                weed(
                    &mut ska_array,
                    weed_file,
                    *reverse,
                    *min_freq,
                    *filter_ambig_as_missing,
                    filter,
                    *ambig_mask,
                    *no_gap_only_sites,
                    if output.is_none() {
                        skf_file
                    } else {
                        output.as_ref().unwrap().as_str()
                    },
                );
            } else {
                panic!("Could not read input file: {skf_file}");
            }
        }
        Commands::Nk {
            skf_file,
            full_info,
        } => {
            if let Ok(ska_array) = MergeSkaArray::<u64>::load(skf_file) {
                println!("{ska_array}");
                if *full_info {
                    log::info!("Printing full info");
                    println!("{ska_array:?}");
                }
            } else if let Ok(ska_array) = MergeSkaArray::<u128>::load(skf_file) {
                println!("{ska_array}");
                if *full_info {
                    log::info!("Printing full info");
                    println!("{ska_array:?}");
                }
            } else {
                panic!("Could not read input file: {skf_file}");
            }
        }
        Commands::Cov {
            fastq_fwd,
            fastq_rev,
            k,
            single_strand,
        } => {
            // Build, merge
            let rc = !*single_strand;
            let cutoff;
            if *k <= 31 {
                log::info!("k={}: using 64-bit representation", *k);
                let mut cov =
                    CoverageHistogram::<u64>::new(fastq_fwd, fastq_rev, *k, rc, args.verbose);
                cutoff = cov.fit_histogram().expect("Couldn't fit coverage model");
                cov.plot_hist();
            } else {
                log::info!("k={}: using 128-bit representation", *k);
                let mut cov =
                    CoverageHistogram::<u128>::new(fastq_fwd, fastq_rev, *k, rc, args.verbose);
                cutoff = cov.fit_histogram().expect("Couldn't fit coverage model");
                cov.plot_hist();
            }
            eprintln!("Estimated cutoff\t{cutoff}");
        }
        Commands::Lo {
            input_skf,
            reference,
            output,
            missing,
            depth,
            indel_kmers,
            threads,
        } => {
            let config = Config {
                input_file: input_skf.clone(),
                output_name: output.clone(),
                max_missing: *missing,
                max_depth: *depth,
                max_indel_kmers: *indel_kmers,
                nb_threads: *threads,
                reference_genome: reference.clone(),
            };

            if let Ok(ska_array) = load_array::<u64>(std::slice::from_ref(input_skf), *threads) {
                log::info!("Reading file {input_skf}");
                log::info!("Using 64-bit representation");
                skalo(ska_array, config);
            } else if let Ok(ska_array) =
                load_array::<u128>(std::slice::from_ref(input_skf), *threads)
            {
                log::info!("Reading file {input_skf}");
                log::info!("Using 128-bit representation");
                skalo(ska_array, config);
            } else {
                panic!("Could not read input file(s): {input_skf:?}");
            }
        }
    }
    let end = Instant::now();

    eprintln!("SKA done in {}s", end.duration_since(start).as_secs());
    eprintln!("⬛⬜⬛⬜⬛⬜⬛");
    eprintln!("⬜⬛⬜⬛⬜⬛⬜");
    log::info!("Complete");
}


// WASM implementation
#[cfg(target_arch = "wasm32")]
#[doc(hidden)]
pub fn main() {
    panic!("You've compiled Ska2 for WebAssembly support, you cannot use it as a normal binary anymore!");
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
extern {
    #[wasm_bindgen(js_namespace = console)]
    fn log(s: &str);
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
/// Function that allows to propagate panic error messages when compiling to wasm, see https://github.com/rustwasm/console_error_panic_hook
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

/// Logging wrapper function for the WebAssembly version
#[cfg(target_arch = "wasm32")]
pub fn logw(text : &str, typ : Option<&str>) {
    if typ.is_some() {
        log((String::from("ska.rust::") + typ.unwrap() + "::" + text).as_str());
    } else {
        log(text);
    }
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
/// Struct to interact with JS when working with WebAssembly
pub struct SkaData {
    k: usize,
    rc: bool,
    reference64: Option<RefSka<u64>>,
    reference128: Option<RefSka<u128>>,
    reference_string: Vec<String>,
    mapped64: Option<Vec<SkaMap<u64>>>,
    mapped128: Option<Vec<SkaMap<u128>>>,
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl SkaData {
    /// Constructor of the SkaData struct
    pub fn new(ref_file: web_sys::File, k: usize) -> Self {
        if cfg!(debug_assertions) {
            init_panic_hook();
        }

        let rc = true;
        let ambig_mask = false;
        let repeat_mask = false;

        let mut wf = WebSysFile::new(ref_file);
        if k < 32 {
            let reference = RefSka::<u64>::new(k, &mut wf, rc, ambig_mask, repeat_mask);
            logw(&format!(
                "Indexed reference: {} split k-mers",
                reference.ksize()
            ), None);
            let reference_string = reconstruct_sequence(&reference);

            Self {
                k: k,
                rc: rc,
                reference64: Some(reference),
                reference128: None,
                reference_string: reference_string,
                mapped64: Some(Vec::new()),
                mapped128: None,
            }
        } else if k < 64 {
            let reference = RefSka::<u128>::new(k, &mut wf, rc, ambig_mask, repeat_mask);
            logw(&format!(
                "Indexed reference: {} split k-mers",
                reference.ksize()
            ), None);
            let reference_string = reconstruct_sequence(&reference);

            Self {
                k: k,
                rc: rc,
                reference64: None,
                reference128: Some(reference),
                reference_string: reference_string,
                mapped64: None,
                mapped128: Some(Vec::new()),
            }
        } else {
            panic!("k values larger than 64 not supported");
        }
    }

    /// Mapping function.
    pub fn map(&mut self, input_file: web_sys::File, rev_reads: Option<web_sys::File>, proportion_reads: Option<f64>) -> String {
        if rev_reads.is_some() {
            logw(&format!("Detected paired fastq input files"), None);
        }
        logw(&format!("Mapping reads to reference"), None);
        let mut input_file_ref = input_file;
        if self.k < 32 {
            if let Some(mut second_file) = rev_reads {
                self.mapped64.as_mut().unwrap().push(SkaMap::<u64>::new(
                    self.reference64.as_ref().unwrap(),
                    &mut input_file_ref,
                    Some(&mut second_file),
                    proportion_reads,
                    self.rc,
                ));
            } else {
                self.mapped64.as_mut().unwrap()
                    .push(SkaMap::<u64>::new(&self.reference64.as_ref().unwrap(), &mut input_file_ref, None, proportion_reads, self.rc));
            };
        } else {
            if let Some(mut second_file) = rev_reads {
                self.mapped128.as_mut().unwrap().push(SkaMap::<u128>::new(
                    &self.reference128.as_ref().unwrap(),
                    &mut input_file_ref,
                    Some(&mut second_file),
                    proportion_reads,
                    self.rc,
                ));
            } else {
                self.mapped128.as_mut().unwrap()
                    .push(SkaMap::<u128>::new(&self.reference128.as_ref().unwrap(), &mut input_file_ref, None, proportion_reads, self.rc));
            };
        }

        logw("Reads mapped successfully!", None);

        let mut results = json::JsonValue::new_array();

        results["Mapped sequences"] = json::JsonValue::new_array();
        let whole_mapping;
        if self.k < 32 {
            whole_mapping = self.reference64.as_ref().unwrap().pseudoalignment(self.mapped64.as_ref().unwrap()[self.mapped64.as_ref().unwrap().len() - 1].get_mapped_bases())[0].clone();
        } else {
            whole_mapping = self.reference128.as_ref().unwrap().pseudoalignment(self.mapped128.as_ref().unwrap()[self.mapped128.as_ref().unwrap().len() - 1].get_mapped_bases())[0].clone();
        }

        let mut current_length = 0;
        for chr in 0..self.reference_string.len() {
            let chr_mapping: String = whole_mapping[current_length..current_length + self.reference_string[chr].len()].to_string();
            let _ = results["Mapped sequences"].push(chr_mapping);
            current_length += self.reference_string[chr].len();
        }

        if self.k < 32 {
            results["Number of variants"] = self.mapped64.as_ref().unwrap()[self.mapped64.as_ref().unwrap().len() - 1].get_mapped_bases().len().into();
        } else {
            results["Number of variants"] = self.mapped128.as_ref().unwrap()[self.mapped128.as_ref().unwrap().len() - 1].get_mapped_bases().len().into();
        }
        // The result is a concatenated string of all the mapped bases

        let mut count_mapped_bases = 0;
        let mut count_total_bases = 0;
        for base in whole_mapping.chars() {
            if base != '-' {
                count_mapped_bases += 1;
                count_total_bases += 1;
            } else {
                count_total_bases += 1;
            }
        }

        results["Coverage"] = (count_mapped_bases as f64 / count_total_bases as f64).into();

        logw("Results computed successfully!", None);

        return results.dump();
    }

    /// Retrieves the reference sequence
    pub fn get_reference(&self) -> String {
        return self.reference_string.join("\n");
    }
}

#[cfg(target_arch = "wasm32")]
/// Reconstructs a sequence beginning from a reference
pub fn reconstruct_sequence<IntT: for<'a> UInt<'a>>(reference: &RefSka<IntT>) -> Vec<String> {
    let sequence_u8 = reference.get_seq();
    let mut sequence_string = Vec::new();

    for sequence in sequence_u8 {
        let mut current_seq= "".to_string();
        for base in sequence {
            if *base != 10 {
                current_seq.push(*base as char);
            }

        }
        sequence_string.push(current_seq);
    }

    return sequence_string;
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
/// AlignData struct for doing alignment while in WebAssembly
pub struct AlignData {
    k: usize,
    alignment64: Option<SkaAlign<u64>>,
    alignment128: Option<SkaAlign<u128>>,
    file_names: Vec<String>,
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
impl AlignData {
    /// Constructor of the AlignData struct
    pub fn new(k: usize) -> Self {
        if k < 32 {
            Self {
                k : k,
                alignment64: Some(SkaAlign::<u64>::new(k)),
                alignment128: None,
                file_names: Vec::new(),
            }
        } else if k < 64 {
            Self {
                k: k,
                alignment64: None,
                alignment128: Some(SkaAlign::<u128>::new(k)),
                file_names: Vec::new(),
            }
        } else {
            panic!("k values larger than 64 are not supported.");
        }
    }

    /// Alignment function.
    pub fn align(&mut self, input_files: Vec<web_sys::File>, proportion_reads: Option<f64>) -> String {
        logw(&format!("Aligning reads"), None);

        // let mut wf1: WebSysFile;

        for input_file in input_files {
            let file_name = input_file.name();
            self.file_names.push(file_name.clone());
            // let mut file_type = file_name.split('.').nth(file_name.split('.').count() - 1).unwrap();
            // if file_type == "gz" {
            //     file_type = file_name.split('.').nth(file_name.split('.').count() - 2).unwrap();
            // }
            // wf1 = WebSysFile::new(input_file);
            if self.k < 32 {
                self.alignment64.as_mut().unwrap().add_file(
                    // &mut wf1,
                    & input_file,
                    proportion_reads,
                );
            } else {
                self.alignment128.as_mut().unwrap().add_file(
                    // &mut wf1,
                    & input_file,
                    proportion_reads,
                );
            }
        }

        if (self.k < 32 && self.alignment64.as_ref().unwrap().get_size() <= 2) || (self.k > 32 && self.alignment128.as_ref().unwrap().get_size() <= 2) {
            let mut results = json::JsonValue::new_array();
            results["newick"] = "Not enough sequences to align".into();
            results["names"] = json::JsonValue::new_array();
            for name in &self.file_names {
                let _ = results["names"].push(name.to_string());
            }

            return results.dump();
        }

        let newick;

        if self.k < 32 {
            newick = self.alignment64.as_mut().unwrap().align(&self.file_names);
        } else {
            newick = self.alignment128.as_mut().unwrap().align(&self.file_names);
        }

        let mut results = json::JsonValue::new_array();

        results["newick"] = newick.into();
        results["names"] = json::JsonValue::new_array();
            for name in &self.file_names {
                let _ = results["names"].push(name.to_string());
            }
        results.dump()
    }
}
