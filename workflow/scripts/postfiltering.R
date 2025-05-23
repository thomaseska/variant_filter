# Genomenexus postfiltering
# clinvar, oncoKB + new output format

library(data.table)
library(dplyr)
library(biomaRt)
library(xlsx)

tbl <- fread(snakemake@input[["csv"]])

filterterms <- c("Conflicting_classifications_of_pathogenicity",
                 "Likely Loss-of-function,Likely Oncogenic",
                 "Likely_pathogenic",
                 "Pathogenic",
                 "Gain-of-function,Oncogenic",
                 "Pathogenic/Likely_pathogenic",
                 "Uncertain_significance",
                 "Loss-of-function,Likely Oncogenic",
                 "Likely Loss-of-function,Oncogenic"
)

tbl <- tbl[(clinvar %in% filterterms | oncokb %in% filterterms),]

tbl <- tbl[, head(.SD, 1), by = "hgvsg"]

form_tbl <- tbl[, c("hugo_gene_symbol", "hgvsg", "hgvsc", "hgvsp", "consequence_terms", "allele_frequency")]

form_tbl[, variant_hg19 := paste(hgvsg, hgvsc, hgvsp, consequence_terms, sep = "; ")]


# get ensembl

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
hugo_ens <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters = "hgnc_symbol",
  values = unique(form_tbl$hugo_gene_symbol),
  mart = ensembl)

form_tbl <- merge(form_tbl, hugo_ens,
                  by.x="hugo_gene_symbol",
                  by.y="hgnc_symbol",
                  all.x=TRUE)

form_tbl[, gene := paste0(hugo_gene_symbol, " (", ensembl_gene_id, ")")]


form_tbl <- form_tbl[, c("gene", "variant_hg19", "allele_frequency")]
form_tbl$interpretation <- NA

fwrite(form_tbl, file=snakemake@output[["csv"]], sep = ",")

write.xlsx(form_tbl, snakemake@output[["xlsx"]], row.names = FALSE, showNA = FALSE)



