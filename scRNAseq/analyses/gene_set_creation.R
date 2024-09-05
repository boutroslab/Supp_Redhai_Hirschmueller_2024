library(here)
library(tidyverse)
library(data.table)
# for both ORA and GSEA we need predefined gene sets.
# we have downloaded them from here: https://maayanlab.cloud/FlyEnrichr/#stats
# in each raw_data directory there is a link to the original source file.

####################
# ENRICHR DATASETS #
####################
files <- list.files(
    path = here("raw_data"),
    pattern = "enrichR_*", recursive = T, full.names = F
)

# they have generated them a bit oddly... for some reason sometimes there are multiple \t
gene_sets <- lapply(files, function(x) {
    name <- str_remove_all(x, ".txt")
    print(x)

    f <- data.table::fread(here("raw_data", x), sep = "\n", header = F)$V1
    df <- data.frame(
        term = trimws(str_match(f, "(.+?)\t")[, 2]),
        genes = trimws(str_match(f, ".+?\t(.+)")[, 2])
    )

    df$genes <- strsplit(df$genes, "\t+")

    df <- unnest(df, genes)
    df$term <- paste(name, df$term, sep = "_")

    return(df)
})
names(gene_sets) <- str_remove_all(files, ".txt")

####################
# FLYBASE PATHWAYS #
####################
flybase_pathways <- data.table::fread(here("raw_data", "pathway_group_data_fb_2023_01.tsv"),
    skip = "Group_member_FB_gene_symbol"
)

flybase_pathways <- flybase_pathways[, c(2, 3, 6, 7)]
colnames(flybase_pathways) <- c(
    "pathway_abbreviation",
    "pathway_name",
    "FBid",
    "symbol"
)
is_empty <- function(x) {
    return(x == "")
}
flybase_pathways <- flybase_pathways %>%
    filter(!if_any(everything(), ~ is_empty(.))) %>% # remove lines with empty entry. they are not NA but ""
    group_by(pathway_name) %>%
    mutate(n = n()) %>%
    filter(n > 5) %>% # remove pathways with less than 5 genes associated with it
    dplyr::select(pathway_name, symbol)

flybase_pathways$symbol <- case_when(
    flybase_pathways$symbol == "Parp" ~ "Parp1",
    flybase_pathways$symbol == "GlyP" ~ "Glyp",
    T ~ flybase_pathways$symbol
)

gene_sets$flybase_pathways <- flybase_pathways

#####################
# HALLMARK GENESETS #
#####################
msigdb_gene_sets <- msigdbr::msigdbr(species = "Drosophila melanogaster")
hallmark_geneset <- msigdb_gene_sets %>%
    filter(gs_cat == "H") %>% # hallmark gene set
    select(pathway_name = gs_name, symbol = gene_symbol)

gene_sets$hallmark <- hallmark_geneset
saveRDS(gene_sets, here("output", "gene_sets.rds"))
