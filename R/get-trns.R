# remotes::install_github("maia-sh/tidypubmed")
# remotes::install_github("maia-sh/ctregistries")

library(dplyr)
library(here)
library(readr)
library(fs)
library(tidypubmed)
library(ctregistries)

# Set up data directories
dir_raw <- dir_create(here("data", "raw"))
dir_processed <- dir_create(here("data", "processed"))

dir_raw_pubmed <- dir_create(path(dir_raw, "pubmed"))
dir_processed_pubmed <- dir_create(path(dir_processed, "pubmed"))

dir_processed_trn <- dir_create(path(dir_processed, "trn"))

dir_raw_ft_xml <- dir_create(path(dir_raw, "fulltext", "xml"))


# Get trials/publications csv ---------------------------------------------

# Manually downloaded trials/publications list: pilot-result_oddpub.csv
# "https://charitede-my.sharepoint.com/:x:/r/personal/tamarinde_haven_bih-charite_de/Documents/Responsible%20Supervision/pilot-result_oddpub.csv?d=wa2f09e107f8346588d1e0b792ecf152c&csf=1&web=1&e=eDVt5H"

df <- read_csv(path(dir_raw, "pilot-result_oddpub.csv"))

pmids <-
  df %>%
  distinct(pmid) %>%
  pull()

pmids_dois <-
  df %>%
  distinct(pmid, doi)

# Download PubMed ---------------------------------------------------------
# Adapted from https://github.com/maia-sh/intovalue-data/blob/main/scripts/02_get-pubmed.R

# If pmids already downloaded, remove those from list to download
if (dir_exists(dir_raw_pubmed)){

  pmids_downloaded <-
    dir_ls(dir_raw_pubmed) %>%
    path_file() %>%
    path_ext_remove() %>%
    as.numeric()

  # Check whether pmids downloaded which aren't needed and manually review and remove
  pmids_downloaded_unused <- setdiff(pmids_downloaded, pmids)
  if (length(pmids_downloaded_unused) > 0) {
    rlang::warn(glue::glue("Unused pmid downloaded: {pmids_downloaded_unused}"))
  }

  pmids <- setdiff(pmids, pmids_downloaded)
}

# Download remaining pmids, if any
if (length(pmids) > 0) {

  # Use pubmed api key locally stored as "ncbi-pubmed", if available
  # Else ask user and store
  pubmed_api_key <-
    ifelse(
      nrow(keyring::key_list("ncbi-pubmed")) == 1,
      keyring::key_get("ncbi-pubmed"),
      keyring::key_set("ncbi-pubmed")
    )

  pmids %>%
    purrr::walk(download_pubmed,
                dir = dir_raw_pubmed,
                api_key = pubmed_api_key
    )

  # Log query date
  loggit::set_logfile(here::here("queries.log"))
  loggit::loggit("INFO", "PubMed")
}


# Extract TRNs from PubMed metadata and abstracts -------------------------


# Adapted from https://github.com/maia-sh/intovalue-data/blob/main/scripts/05_extract-trns.R

# PubMed xmls should be available here
pubmed_xmls <- dir_ls(dir_raw_pubmed)

# Extract secondary identifiers from pubmed
si <-
  pubmed_xmls %>%
  purrr::map_dfr(extract_pubmed, datatype = "databanks", quiet = FALSE) %>%
  ctregistries::mutate_trn_registry(accession_number)

write_csv(si, path(dir_processed_pubmed, "pubmed-si.csv"))

# Visually inspect mismatching trns and registries
si_trn_mismatches <-
  si %>%
  filter(!accession_number %in% trn |
           !databank %in% registry)

si_trns <-
  si %>%
  tidyr::drop_na(trn) %>%
  select(pmid, registry, trn_detected = trn) %>%
  distinct() %>%
  group_by(pmid) %>%
  mutate(n_detected = row_number()) %>%
  ungroup() %>%
  mutate(
    source = "secondary_id",
    # Note: `tidypubmed` does not currently clean CTRI so this doesn't run
    # trn_cleaned = purrr::map_chr(trn_detected, ctregistries::clean_trn)
  ) %>%
  left_join(pmids_dois, by  = "pmid")

write_csv(si_trns, path(dir_processed_trn, "trn-si.csv"))

# Extract abstracts from pubmed
abs <-
  pubmed_xmls %>%
  purrr::map_dfr(extract_pubmed, datatype = "abstract", quiet = FALSE) %>%
  ctregistries::mutate_trn_registry(abstract)

write_csv(abs, path(dir_processed_pubmed, "pubmed-abstract.csv"))

abs_trns <-
  abs %>%
  tidyr::drop_na(trn) %>%
  distinct(pmid, registry, trn_detected = trn) %>%
  group_by(pmid) %>%
  mutate(n_detected = row_number()) %>%
  ungroup() %>%
  mutate(
    source = "abstract",
    # Note: `tidypubmed` does not currently clean CTRI so this doesn't run
    # trn_cleaned = purrr::map_chr(trn_detected, ctregistries::clean_trn)
  ) %>%
  left_join(pmids_dois, by  = "pmid")

write_csv(abs_trns, path(dir_processed_trn, "trn-abstract.csv"))


# Extract TRNs from full-text (xmls) --------------------------------------
# Manually downloaded from:
# https://charitede-my.sharepoint.com/:f:/r/personal/tamarinde_haven_bih-charite_de/Documents/Responsible%20Supervision/2022-02-01-resp-sup/2022-02-01-resp-sup-xml?csf=1&web=1&e=pzkuGl

# Copied from https://github.com/maia-sh/intovalue-data/blob/main/scripts/functions/get_grobid_ft_trn.R
source(here("R", "get_grobid_ft_trn.R"))

ft_xmls <- dir_ls(dir_raw_ft_xml)

ft_trns <-
  ft_xmls %>%
  purrr::map_dfr(get_grobid_ft_trn) %>%
  select(-pmid) %>%
  rename(trn_detected = trn, n_detected = n) %>%
  mutate(
    source = "ft",
    # trn_cleaned = purrr::map_chr(trn_detected, ctregistries::clean_trn)
  ) %>%
  left_join(pmids_dois, by  = "doi")

write_csv(ft_trns, path(dir_processed_trn, "trn-ft.csv"))

# Combine reported TRNs (secondary id, abstract, full-text) ---------------

trn_combined <-
  bind_rows(si_trns, abs_trns, ft_trns) %>%

  distinct(pmid, doi, trn = trn_detected, registry, source) %>%

  # All records should have a trn
  assertr::assert(assertr::not_na, trn)

write_csv(trn_combined, path(dir_processed_trn, "trn-reported-long.csv"))

# Pivot wider to for one row per TRN with sources as columns --------------

trn_combined_wide <-
  trn_combined %>%

  # Note: `value_fill` is FALSE for all but some will be replaced with NA, e.g., because no full-text
  mutate(value = TRUE) %>%
  tidyr::pivot_wider(
    names_from = source, names_prefix = "has_trn_",
    values_from = value, values_fill = FALSE
  )

write_csv(trn_combined_wide, path(dir_processed_trn, "trn-reported-wide.csv"))
