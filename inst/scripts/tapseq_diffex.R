## code to prepare `tapseq_diffex` dataset goes here

## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135497

## Creates a SingleCellExperiment of TAP_DIFFEX samples 1 and 2
## from the TAPseq paper (Schraivogel et al 2020)
## (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135497)

## The data is downloaded via the Python crispat package which
## bundles it, which in turn originally downloaded it from the
## following URLs:
## http://steinmetzlab.embl.de/TAPdata/TAP.nods.RDS
## http://steinmetzlab.embl.de/TAPdata/Whole.nods.RDS.

## Before running this script, clone the crispat repo as follows:
## git clone https://github.com/velten-group/crispat ~/src/crispat

library(magrittr)

library(Matrix)

library(SingleCellExperiment)

# Load gene expression counts

cnt_gex <- read.csv(
    '~/src/crispat/example_data/Schraivogel/gene_expression_counts.csv'
)

rownames(cnt_gex) <- cnt_gex[, 1]
cnt_gex <- cnt_gex[, -1]

cnt_gex <- as.matrix(cnt_gex)
cnt_gex <- as(cnt_gex, 'CsparseMatrix')

cnt_gex <- t(cnt_gex)

# Load gRNA counts

cnt_grna <- read.csv(
    '~/src/crispat/example_data/Schraivogel/gRNA_counts.csv',
    check.names = F
)

rownames(cnt_grna) <- cnt_grna[, 1]
cnt_grna <- cnt_grna[, -1]

cnt_grna <- as.matrix(cnt_grna)
cnt_grna <- as(cnt_grna, 'CsparseMatrix')

# Generate gRNA metadata from rownames

grna_rowdat <- stringr::str_match(
    rownames(cnt_grna),
    "CROPseq_dCas9_DS_(.*)"
)

colnames(grna_rowdat) <- c("grna_id", "grna_name")

grna_rowdat %<>%
    as.data.frame() %>%
    dplyr::mutate(
        # simplify NTC names
        grna_name = dplyr::case_when(
            startsWith(grna_name, 'non-targeting') ~ sprintf(
                "NTC-%02d",
                as.integer(stringr::str_match(
                    grna_name,
                    'non-targeting_(\\d+)'
                )[, 2])
            ),
            TRUE ~ grna_name
        )
    ) %>%
    dplyr::mutate(
        target_gene = stringr::str_match(
            grna_name,
            "([^-_]+)[-_].*"
        )[, 2]
    ) %>%
    dplyr::mutate(
        # https://www.nature.com/articles/s41592-020-0837-5
        # Supplemental Table 2, Chr 8 control library
        target_element_type = dplyr::case_when(
            target_gene == 'NTC' ~ 'NTC',
            target_gene %in% c('GATA1', 'HS2', 'MYC', 'ZFPM2') ~ 'enhancer',
            TRUE ~ 'promoter'
        )
    ) %>%
    dplyr::mutate(
        # simplify promoter names
        grna_name = dplyr::if_else(
            target_element_type == 'promoter',
            sprintf(
                "%s-%d",
                target_gene,
                as.integer(stringr::str_match(
                    grna_name,
                    "^[^_]+_[+-]_(\\d+)(\\.23)?(-P1P2)?"
                )[, 2])
            ),
            grna_name
        )
    ) %>%
    `rownames<-`(.$grna_name)

stopifnot(rownames(cnt_grna) == grna_rowdat$grna_id)
rownames(cnt_grna) <- rownames(grna_rowdat)

stopifnot(colnames(cnt_gex) == colnames(cnt_grna))

# Extract cell metadata from column names

coldat <- stringr::str_match(
    colnames(cnt_gex),
    "^(TAP\\d+)-[ACTG]+$"
)

colnames(coldat) <- c('cell_id', 'sample')
coldat %<>% as.data.frame()
rownames(coldat) <- coldat$cell_id

# Create and save the object

tapseq_diffex <- SingleCellExperiment(
    assays = list(counts = cnt_gex),
    colData = coldat,
    altExps = list(
        grna = SingleCellExperiment(
            list(counts = cnt_grna),
            rowData = grna_rowdat
        )
    )
)

usethis::use_data(tapseq_diffex, overwrite = TRUE)
