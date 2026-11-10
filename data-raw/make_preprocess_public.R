# =============================================================================
# make_preprocess_public.R: PUBLIC replacements for Data_preprocessI and Data_preprocessS.
# Ray, 2026-11-09: the earlier versions held data from a controlled-access cohort (a GWAS, an
# in-sample LD matrix and a cohort-trained model) and are removed. Same structure and column names, so the
# preprocessing vignette (Appendix_Preprocess.Rmd) runs unchanged.
#
#   Rscript make_preprocess_public.R <out_dir>        (on psoriasis; needs the public resources)
#
# Every input is public:
#   SNPs       the public HapMap3 AFR set (resources/public_resource/zenodo/panAFR_example/data/
#              analysis_set.tsv: PRS-CS UKBB HapMap3, MAF >= 1% in 1000G and UKBB AFR, positions in
#              GRCh38): the 10,000 with the smallest Pan-UKB AFR p-values for the target; every 8th SNP
#              plus those for the larger external (about 157,000)
#   target.ss  Pan-UKB AFR height (continuous-50-both_sexes-irnt): BETA for the A1 (= ALT) allele,
#              STAT = BETA / SE, P, NMISS = the Pan-UKB AFR sample size
#   target.ld  the UKBB AFR LD panel (PRS-CS blocks, signs for ALT) on the target SNPs
#   external.coef1, external.coef2  two PGS Catalog height scores (GRCh38 harmonized files) in THEIR
#              OWN allele orientation and per-allele scale, so alignment has real flips to resolve;
#              scores developed on UKB AFR are excluded; coef1 on the target SNPs, coef2 on the
#              larger set (merging and zero padding)
#   target.info (Data_preprocessI)  the target SNPs as PLINK-style annotation (A1 = ALT)
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(Matrix) })
ROOT = "/net/psoriasis/psorgenom/zrayw/BRIER_software"
PUB = file.path(ROOT, "resources", "public_resource")
Z = file.path(PUB, "zenodo", "panAFR_example", "data")
PGS_SRC = file.path(ROOT, "manuscript_experiments", "MGI_application", "traits", "height", "externals_pgs")
OUT = commandArgs(trailingOnly = TRUE)[1]; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

a = fread(file.path(Z, "analysis_set.tsv")); tg = fread(file.path(Z, "target_AFR.tsv"))
stopifnot(identical(a$SNP, tg$SNP))
meta = jsonlite::fromJSON(file.path(Z, "target_meta.json"))
# target: the 10,000 SNPs with the smallest Pan-UKB AFR p-values (as the original example's "top
# 10,000 SNPs by marginal association"; an evenly thinned set carried no signal: GIC chose the empty
# model at every eta), in genome order; larger external set: every 8th SNP plus those 10,000
pz = 2 * stats::pnorm(-abs(tg$beta / tg$se))
small = sort(order(pz)[1:10000]); large = sort(unique(c(seq(1L, nrow(a), by = 8L), small)))
cat(sprintf("target SNPs %d, larger external set %d\n", length(small), length(large)))

# target summary statistics and SNP annotation (A1 = ALT, the effect allele)
ts = data.frame(SNP = a$SNP[small], CHR = a$CHR[small], BP = a$BP_hg38[small], A2 = a$REF[small], A1 = a$ALT[small],
                BETA = tg$beta[small], STAT = tg$beta[small] / tg$se[small],
                P = 2 * stats::pnorm(-abs(tg$beta[small] / tg$se[small])), NMISS = as.integer(meta$N))
ld_info = ts[, c("SNP", "CHR", "BP", "A2", "A1")]

# UKBB AFR LD on the target SNPs
blocks = lapply(sort(unique(ts$CHR)), function(cc) {
  R = readRDS(file.path(PUB, "public_ld", "panels", "AFR", "UKBB", sprintf("LD_chr%d.rds", cc)))
  j = which(ts$CHR == cc); k = match(ts$SNP[j], rownames(R)); stopifnot(!anyNA(k))
  R[k, k, drop = FALSE]
})
LD = methods::as(methods::as(Matrix::bdiag(blocks), "generalMatrix"), "CsparseMatrix")
dimnames(LD) = NULL
cat(sprintf("LD %d x %d, nonzero %d\n", nrow(LD), ncol(LD), length(LD@x)))

# two PGS Catalog scores, in their own orientation, restricted by position to a SNP subset
excl = fread(file.path(PUB, "panAFR_example", "data", "pgs_exclusion.tsv"))
keep = excl[exclude == 0, pgs_id]
read_score = function(id, rows) {
  f = fread(file.path(PGS_SRC, paste0(id, "_GRCh38.txt.gz")))
  k = match(paste(a$CHR[rows], a$BP_hg38[rows]), paste(f$CHR, f$BP))
  out = data.frame(CHR = as.character(a$CHR[rows]), BP = as.character(a$BP_hg38[rows]),
                   REF = ifelse(is.na(k), a$REF[rows], f$REF[k]), ALT = ifelse(is.na(k), a$ALT[rows], f$ALT[k]),
                   coef = ifelse(is.na(k), 0, f$coef[k]))
  out$coef[is.na(out$coef)] = 0
  out
}
ids = strsplit(Sys.getenv("PP_PGS", ""), ",")[[1]]
if (length(ids) == 2L) {
  id1 = ids[1]; id2 = ids[2]            # 2026-11-09 build: PGS004995, PGS004998 (the two best-covered)
} else {
  cov_small = sapply(keep, function(id) { sum(read_score(id, small)$coef != 0) })
  id1 = names(which.max(cov_small)); id2 = names(sort(cov_small, decreasing = TRUE))[2]
}
stopifnot(id1 %in% keep, id2 %in% keep)
c1 = read_score(id1, small); c2 = read_score(id2, large)
# The harmonized GRCh38 files already share the panel's REF / ALT, so they hold no allele flips. To
# keep the alignment demonstration, external.coef1 is written in the common "effect allele = the
# trait-increasing allele" convention: where the weight is negative, REF and ALT are swapped and the
# sign is flipped. The same model, written another way; preprocessI / preprocessS must undo it.
neg = c1$coef < 0
c1[neg, c("REF", "ALT")] = c1[neg, c("ALT", "REF")]
c1$coef[neg] = -c1$coef[neg]
cat(sprintf("external.coef1 = %s (%d nonzero of %d, %d allele swaps vs the target), external.coef2 = %s (%d nonzero of %d)\n",
            id1, sum(c1$coef != 0), nrow(c1), sum(c1$ALT != ts$A1), id2, sum(c2$coef != 0), nrow(c2)))
attr(c1, "source") = paste0("PGS Catalog ", id1, " (GRCh38 harmonized)")
attr(c2, "source") = paste0("PGS Catalog ", id2, " (GRCh38 harmonized)")

Data_preprocessS = list(target.ss = ts, target.ld = ld_info, target.ld.mat = LD, external.coef1 = c1, external.coef2 = c2)
Data_preprocessI = list(target.info = ld_info, external.coef1 = c1, external.coef2 = c2)
save(Data_preprocessS, file = file.path(OUT, "Data_preprocessS.rda"), compress = "xz")
save(Data_preprocessI, file = file.path(OUT, "Data_preprocessI.rda"), compress = "xz")
writeLines(c(id1, id2), file.path(OUT, "pgs_ids.txt"))
cat(sprintf("sizes: S %.0f KB, I %.0f KB\n", file.size(file.path(OUT, "Data_preprocessS.rda")) / 1024,
            file.size(file.path(OUT, "Data_preprocessI.rda")) / 1024))
