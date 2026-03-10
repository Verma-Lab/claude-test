#!/usr/bin/env Rscript
# plink_allele_freq.R
# Calculate allele frequencies for a selected list of variants from PLINK files.
#
# Usage:
#   Rscript plink_allele_freq.R \
#     --bfile <plink_prefix>       # prefix for .bed/.bim/.fam (mutually exclusive with --vcf)
#     --vcf <vcf_file>             # input VCF file (mutually exclusive with --bfile)
#     --variants <variant_list>    # file with one variant ID per line
#     --out <output_prefix>        # output file prefix
#     [--plink2 <path>]            # path to plink2 binary (default: plink2)
#     [--maf-filter <value>]       # optional MAF filter threshold (e.g. 0.01)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--bfile",      type = "character", default = NULL,
              help = "PLINK binary file prefix (.bed/.bim/.fam)"),
  make_option("--vcf",        type = "character", default = NULL,
              help = "Input VCF/BCF file"),
  make_option("--variants",   type = "character", default = NULL,
              help = "File containing variant IDs to extract (one per line) [required]"),
  make_option("--out",        type = "character", default = NULL,
              help = "Output file prefix [required]"),
  make_option("--plink2",     type = "character", default = "plink2",
              help = "Path to plink2 binary [default: plink2]"),
  make_option("--maf-filter", type = "double",    default = NULL,
              help = "Minor allele frequency filter threshold")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------------------------------------------------------
# Validate required arguments
# ---------------------------------------------------------------------------
errors <- character(0)

if (is.null(opt$variants)) errors <- c(errors, "--variants is required")
if (is.null(opt$out))      errors <- c(errors, "--out is required")
if (is.null(opt$bfile) && is.null(opt$vcf))
  errors <- c(errors, "Either --bfile or --vcf must be provided")
if (!is.null(opt$bfile) && !is.null(opt$vcf))
  errors <- c(errors, "--bfile and --vcf are mutually exclusive")

if (!is.null(opt$bfile)) {
  for (ext in c(".bed", ".bim", ".fam")) {
    f <- paste0(opt$bfile, ext)
    if (!file.exists(f)) errors <- c(errors, sprintf("Missing file: %s", f))
  }
}
if (!is.null(opt$vcf) && !file.exists(opt$vcf))
  errors <- c(errors, sprintf("VCF file not found: %s", opt$vcf))
if (!is.null(opt$variants) && !file.exists(opt$variants))
  errors <- c(errors, sprintf("Variant list file not found: %s", opt$variants))
if (length(errors) > 0) {
  cat("ERROR(s):\n")
  cat(paste0("  ", errors, "\n"), sep = "")
  quit(status = 1)
}

# ---------------------------------------------------------------------------
# Build plink2 command
# ---------------------------------------------------------------------------
tmpdir  <- tempfile("plink_freq_")
dir.create(tmpdir, showWarnings = FALSE)
tmp_out <- file.path(tmpdir, "extract")

args <- character(0)

if (!is.null(opt$bfile)) {
  args <- c(args, "--bfile", opt$bfile)
} else {
  args <- c(args, "--vcf", opt$vcf)
}

args <- c(args,
  "--extract",   opt$variants,
  "--freq",
  "--out",       tmp_out,
  "--no-psam-pheno"
)

if (!is.null(opt[["maf-filter"]])) args <- c(args, "--maf", opt[["maf-filter"]])

cmd <- paste(c(shQuote(opt$plink2), args), collapse = " ")
cat("Running:\n  ", cmd, "\n\n")

# ---------------------------------------------------------------------------
# Run plink2
# ---------------------------------------------------------------------------
ret <- system(cmd)
if (ret != 0) {
  cat(sprintf("ERROR: plink2 exited with status %d\n", ret))
  quit(status = ret)
}

# plink2 --freq writes <out>.afreq
freq_file <- paste0(tmp_out, ".afreq")
if (!file.exists(freq_file)) {
  cat("ERROR: Expected output file not found:", freq_file, "\n")
  quit(status = 1)
}

# ---------------------------------------------------------------------------
# Load and filter results
# ---------------------------------------------------------------------------
freq <- fread(freq_file)

# Read requested variant IDs
requested_ids <- fread(opt$variants, header = FALSE, col.names = "ID")

# Standardise column name for variant ID (plink2 uses #CHROM/ID or ID)
id_col <- intersect(c("ID", "SNP"), colnames(freq))
if (length(id_col) == 0) {
  cat("ERROR: Cannot identify variant ID column in plink2 output. Columns:",
      paste(colnames(freq), collapse = ", "), "\n")
  quit(status = 1)
}
id_col <- id_col[1]

freq_filtered <- freq[get(id_col) %in% requested_ids$ID]

missing_vars <- setdiff(requested_ids$ID, freq_filtered[[id_col]])
if (length(missing_vars) > 0) {
  cat(sprintf("WARNING: %d variant(s) from the list were not found in the data:\n",
              length(missing_vars)))
  cat(paste0("  ", missing_vars, "\n"), sep = "")
}

# ---------------------------------------------------------------------------
# Write output
# ---------------------------------------------------------------------------
out_dir <- dirname(opt$out)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- paste0(opt$out, ".afreq.tsv")
fwrite(freq_filtered, out_file, sep = "\t")
cat(sprintf("Output written to: %s\n", out_file))
cat(sprintf("Variants in output: %d / %d requested\n",
            nrow(freq_filtered), nrow(requested_ids)))

# ---------------------------------------------------------------------------
# Cleanup temp files
# ---------------------------------------------------------------------------
unlink(tmpdir, recursive = TRUE)
