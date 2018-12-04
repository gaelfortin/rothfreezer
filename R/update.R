#' Update a CSV file from a google spreadsheet
#' 
#' @param file Where to save the CSV.
#' @param from Where to download the CSV. (Must be an accessible google spreadsheet URL)
#' 
#' @importFrom googlesheets gs_url gs_download
#' @export

update_from_googlesheet <- function(file, from) {
  googlesheets::gs_download(googlesheets::gs_url(from), to = file, overwrite = TRUE)
}


#' Update gene annotation table
#' 
#' Update *S.cerevisiae* gene annotations from 
#' [SGD](http://www.yeastgenome.org)
#' 
#' @param file Where to write the CSV.
#' @param from Where to download gene annotations from.
#' 
#' @importFrom readr read_tsv write_csv cols_only col_character
#' @importFrom rlang .data
#' @export

update_genes <- function(file = 'genes/sgd.csv', 
                         from = 'http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab') {
  from %>%
    readr::read_tsv(
      col_names = c(
        'gene_id', 'type', 'qualifier', 'orf_id', 'gene_name', 'gene_alias',
        'parent_feature_name', 'secondary_sgd_id', 'chromosome', 'start', 'stop',
        'strand', 'genetic_position', 'coordinate_version', 'sequence_version', 
        'gene_description'
      ), 
      col_types = readr::cols_only(
        gene_id          = readr::col_character(), 
        orf_id           = readr::col_character(), 
        gene_name        = readr::col_character(), 
        gene_alias       = readr::col_character(), 
        gene_description = readr::col_character()
      )
    ) %>%
    mutate( 
      gene_name = ifelse(is.na(.data$gene_name), .data$orf_id, .data$gene_name)
    ) %>%
    select_at(c('gene_id', 'orf_id', 'gene_name', 'gene_alias', 'gene_description')) %>%
    write_csv(file)
}

#' Update genetic distances
#' 
#' Update *S.cerevisiae* genetic distances from 
#' [SGD](http://www.yeastgenome.org)
#' 
#' @param dir Where to write the CSVs. Each chromosome is saved in a separate file.
#' @param from Where to download genetic distances from.
#' 
#' @importFrom readr read_tsv write_csv
#' @importFrom rlang .data
#' @export
#' @md

update_genetic_distances <- function(dir = 'genetic_distances',
                                     from = 'http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab') {
  
  # Get location information for all ORFs on SGD
  sgd <- 
    from %>%
    read_tsv(
      col_types = 'cc_cc___ciic____',
      col_names = c(
        'gene_id', 'type', 'orf_id', 'gene_name', 'chromosome', 'start', 'stop',
        'strand')
    ) %>%
    filter(.data$type %in% c('ORF', 'centromere')) %>%
    mutate(
      gene_name = ifelse(is.na(.data$gene_name), .data$orf_id, .data$gene_name),
      max = ~pmax(.data$start, .data$stop),
      min = ~pmin(.data$start, .data$stop)
    )
  
  # Assign each ORF to a chromosome arm
  with_chromosome_arm <-
    left_join(
      sgd %>% filter(.data$type != 'centromere'),
      sgd %>% filter(.data$type == 'centromere') %>% select('chromosome', cen_min = 'min', cen_max = 'max')
    ) %>%
    mutate(arm = ifelse(.data$max < .data$cen_min, 'left', 'right'))
  
  # Pair each ORF to all ORFs on the same chromosome arm and compute the distance
  gene_pairs <-
    full_join(
      select('with_chromosome_arm', 'chromosome', 'arm', gene_id_a = 'gene_id', orf_id_a = 'orf_id', gene_name_a = 'gene_name', a_max = 'max', a_min = 'min'),
      select('with_chromosome_arm', 'chromosome', 'arm', gene_id_b = 'gene_id', orf_id_b = 'orf_id', gene_name_b = 'gene_name', b_max = 'max', b_min = 'min')
    ) %>%
    mutate(
      end_to_start = .data$a_max - .data$b_min,
      start_to_end = .data$a_min - .data$b_max,
      distance = ~pmin(abs(.data$end_to_start), abs(.data$start_to_end)),
      distance = ~ifelse(sign(.data$end_to_start) != sign(.data$start_to_end), 0, .data$distance)
    ) %>%
    select('chromosome', 'arm', 'gene_id_a', 'gene_id_b', 'orf_id_a', 'orf_id_b', 'gene_name_a', 'gene_name_b', 'distance')
  
  # Write to file
  by_chromosome <- split(gene_pairs, paste(gene_pairs$chromosome, gene_pairs$arm))
  lapply(by_chromosome, function(x) { 
    write_csv(x, paste0(dir, '/', x$chromosome[1], '-', x$arm[1], '.csv'))
  })
}


# ---- Shared resources ----
#' @rdname update_from_googlesheet
#' @export
update_antibodies <- function(file = 'antibodies/antibodies.csv', 
                              from = 'https://docs.google.com/spreadsheets/d/1gyiZth6awGqdB6YiQqDLL9WSha7KxI525GYY6aoxcCI/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_primers <- function(file = 'primers/primers.csv', 
                           from = 'https://docs.google.com/spreadsheets/d/1tYxhHsEN4mBGJRZXqJS7vK7TJJdP7Wayv4P14R9wTMI/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_media <- function(file = 'media/media.csv', 
                         from = 'https://docs.google.com/spreadsheets/d/1mG2Vr1sMnNAIX-NCFC33O0WoUVDrI3F6no78Czs6Y8U/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}


# ---- Experiment annotation ----
#' @rdname update_from_googlesheet
#' @export
update_queries <- function(file = 'queries/queries.csv', 
                           from = 'https://docs.google.com/spreadsheets/d/1VgVEg1iJrbjOU6kh_9i4ZxmINBaVkrxoF-kGl0m51WM/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_treatments <- function(file = 'treatments/treatments.csv',
                              from = 'https://docs.google.com/spreadsheets/d/1pWgFVl_rQ2HUtPSIthf0REjGvLNEUfzr7ob9rRAw2NI/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}


# ---- Strain collections / information ----
#' @rdname update_from_googlesheet
#' @export
update_strain_collections <- function() {
  # Cisplatin libraries
  update_strain_collection_cispt_sensitive_1()
  update_strain_collection_MATa_cispt_384_hb_V2()
  update_strain_collection_MATa_cispt_384_hb()
  update_strain_collection_MATa_cispt_96_hb_V2()
  update_strain_collection_MATa_cispt_96_hb()
  
  # Other libraries
  update_strain_collection_MATalpha_CTF4_RAD5_384_hb()
  update_strain_collection_MATalpha_CTF4_RAD5_96_hb()
  
  # Community libraries
  update_strain_collection_MATa_deletion_384()
  update_strain_collection_MATalpha_deletion_384_hb()
  update_strain_collection_MATalpha_deletion_96()
  update_strain_collection_MATa_ts_384_hb()
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_cispt_sensitive_1 <- function(file = 'strain-collections/cispt-sensitive-1.csv', 
                                                       from = 'https://docs.google.com/spreadsheets/d/1vTC2im5t2u-PW4QJNNia_M2fsAVXqodOj3UPmhJ7-4A/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATa_cispt_384_hb_V2 <- function(file = 'strain-collections/MATa-cispt-384-hb-V2.csv', 
                                                          from = 'https://docs.google.com/spreadsheets/d/1nHFLrPqMODgpJaBkkcDniEnp8Cq7tA8iOowvCPI5wuQ/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATa_cispt_384_hb <- function(file = 'strain-collections/MATa-cispt-384-hb.csv', 
                                                       from = 'https://docs.google.com/spreadsheets/d/1r9VKumSwjGvEcBP8XebfCsXJnKDHvh_AnCKiiBty3rk/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATa_cispt_96_hb_V2 <- function(file = 'strain-collections/MATa-cispt-96-hb-V2.csv', 
                                                         from = 'https://docs.google.com/spreadsheets/d/1eCjWKnrF-H9A0qdaufYuir3CdIy-ueMffLv3MVezjp0/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATa_cispt_96_hb <- function(file = 'strain-collections/MATa-cispt-96-hb.csv', 
                                                      from = 'https://docs.google.com/spreadsheets/d/1f9XuSwpJeN-j5_wQ3mwUGOlsGH-dixgBnHeamzl7PzU/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATa_deletion_384 <- function(file = 'strain-collections/MATa-deletion-384.csv', 
                                                       from = 'https://docs.google.com/spreadsheets/d/1dl89Nf4oI9c6lhAxAfApdM1-mz-RqkzZbrE4DOPQjxk/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATalpha_deletion_384_hb <- function(file = 'strain-collections/MATalpha-deletion-384-hb.csv', 
                                                              from = 'https://docs.google.com/spreadsheets/d/1yyBZh8OEfjmR0oMENEqgCgUj03STMELHIZ08_Kl_fHQ/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATalpha_deletion_96 <- function(file = 'strain-collections/MATalpha-deletion-96.csv', 
                                                          from = 'https://docs.google.com/spreadsheets/d/1whG-KPW9htHvIR8LutivZpW3w3pEBSdEtJpg66aOaq0/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATalpha_CTF4_RAD5_384_hb <- function(file = 'strain-collections/MATalpha-CTF4,RAD5-384-hb.csv', 
                                                               from = 'https://docs.google.com/spreadsheets/d/145J9KpYN38Ra8TAw2-dvAc_g80EJ_iAeWtqqp-NlPf8/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATalpha_CTF4_RAD5_96_hb <- function(file = 'strain-collections/MATalpha-CTF4,RAD5-96-hb.csv', 
                                                              from = 'https://docs.google.com/spreadsheets/d/19CJrfjAgspW-zE2dkRa9nvQDsrfuWYO45dc3PSxiQII/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_MATa_ts_384_hb <- function(file = 'strain-collections/MATa-ts-384-hb.csv', 
                                                    from = 'https://docs.google.com/spreadsheets/d/1zlM9wxoVHtSk5yac209veMCGZPvOnGRUdt5OJSqdGmA/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strain_collection_info <- function(file = 'strain-collection-info/strain-collection-info.csv', 
                                          from = 'https://docs.google.com/spreadsheets/d/12osNdQV8Al6lVhvssnvnVajLXxMZ7qH6Bv8bwvrUg5c/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}


# ---- Plasmids ----
#' @rdname update_from_googlesheet
#' @export
update_plasmids <- function() {
  update_plasmids_pWJ()
  update_plasmids_ZOO()
}

#' @rdname update_from_googlesheet
#' @export
update_plasmids_pWJ <- function(file = 'plasmids/pWJ.csv', 
                                from = 'https://docs.google.com/spreadsheets/d/15DXiQUBR4QQsthB7v9IU84ZcmiZB-kBKfrOjpiyisvw/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_plasmids_ZOO <- function(file = 'plasmids/ZOO.csv', 
                                from = 'https://docs.google.com/spreadsheets/d/1P-F7s-SGvey6B3_PWf-HGv07iwjkEqhwuz-muSyQuyU/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}


# ---- Strains ----
#' @rdname update_from_googlesheet
#' @export
update_strains <- function() {
  update_strains_J()   # Mutagenized
  update_strains_R()   # Foreign labs
  update_strains_rec() # Deletion consortium
  update_strains_special()
  update_strains_TS()  # Temperature sensitive
  update_strains_U()   # Transformed
  update_strains_W()   # Crosses
}

#' @rdname update_from_googlesheet
#' @export
update_strains_special <- function(file = 'strains/special.csv',
                                   from = 'https://docs.google.com/spreadsheets/d/1M-3LGm8P4X_fQ4hYO9BRha1JvwrmTHVExGRAwHXA_kw/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strains_J <- function(file = 'strains/J.csv', 
                             from = 'https://docs.google.com/spreadsheets/d/1XaqVTIWVFpV4PWi2-N8qUrZJZpxuRfDUqdB8CdsmnMw/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strains_R <- function(file = 'strains/R.csv', 
                             from = 'https://docs.google.com/spreadsheets/d/1UT1e4LyINK1F_z_koIpNC217R-om5fKL1IdwRPt0CLI/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strains_rec <- function(file = 'strains/rec.csv', 
                               from = 'https://docs.google.com/spreadsheets/d/19UL9r7UNu72PAZ6fqXIJ8TpP1YB5JfCcpcknPCrw0W8/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strains_TS <- function(file = 'strains/TS.csv', 
                              from = 'https://docs.google.com/spreadsheets/d/1J4MvBQTcJJLk2umR5jhsXPqn1SZ-SNWi763nBq7mNQY/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strains_U <- function(file = 'strains/U.csv',
                             from = 'https://docs.google.com/spreadsheets/d/1sj6cKonhGEcsE65pt_0QnIageBLEwv2MG2BW-b5NwoM/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}

#' @rdname update_from_googlesheet
#' @export
update_strains_W <- function(file = 'strains/W.csv', 
                             from = 'https://docs.google.com/spreadsheets/d/1sDqH44-LWbf4_loE_BVDxsehEoiZ44E93apYCtpgs1M/edit?usp=sharing') {
  update_from_googlesheet(file, from)
}
