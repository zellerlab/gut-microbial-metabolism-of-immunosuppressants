# load first few lines of /g/scb/zeller/karcher/PRISMA/scripts/PRISMA/2024_11_27_prepare_PRISMA_WGS.r
# then execute this 
meta  %>% 
	select(PSN, sampleID, visit, batch) %>% 
	distinct() %>%  
	right_join(read_csv('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/MPA_Tox_df_for_Nic.csv') %>% 
	rename(PSN = v65_pat_id, visit = v62_visit_number) %>% 
	select(PSN, visit)) %>% 
	mutate(
		batch_orig = case_when(
			batch == "modelling" ~ "240813_MG011_PRISMA_NovaSeq_final", 
			batch == "interim_Batch1" ~ "230111_PRISMA_Batch1_P3", 
			batch == "interim_Batch2" ~ "230727_PRISMA_Batch2_P3", 
			batch == "interim_Batch4" ~ "230804_PRISMA_Batch4_P3")
		) %>% 
	mutate(
		fastq_base_path = case_when(
			#str_detect(batch, "interim") ~ "/g/scb/zeller/karcher/PRISMA/data/WGS/",
			batch == "interim_Batch1" ~ "/g/scb/zeller/karcher/PRISMA/data/WGS/230111_PRISMA_Batch1_P3/",
			batch == "interim_Batch2" ~ "/g/scb/zeller/karcher/PRISMA/data/WGS/230727_PRISMA_Batch2_P3/",
			batch == "interim_Batch4" ~ "/g/scb/zeller/karcher/PRISMA/data/WGS/230804_PRISMA_Batch4_P3/",
			str_detect(batch, "modelling") ~ "/g/scb/mzimmerm/raw_GeneCore_data/2024/240813_MG011_PRISMA_NovaSeq/"
		)
	) %>%
	mutate(sampleID_fixed = str_replace(sampleID, "MG_", "MG")) %>%
	mutate(file_suffix = case_when(
		str_detect(batch, "interim") ~ "_R1_001.fastq.gz",
		str_detect(batch, "modelling") ~ "_R1_001.fastq.gz"
	)) %>%
	mutate(fastq_names = map2(
		fastq_base_path,
		sampleID_fixed, 
		\(ff, ss) {
			if (is.na(ff)) {
				return(NA)
			}
			if (str_detect(ff, "Batch4")) {
				ss <- str_replace(ss, "MG", "")
			}
			if (!str_detect(ff,  "Batch")) {
				tmp <- list.files(ff, full.names = TRUE)
				tmp <- tmp[str_detect(tmp, str_c(".*lane1", ss, ".*"))]
				tmp <- tmp[str_detect(tmp, "sequence.txt.gz")]
				#stopifnot(length(tmp) == 2)
				if (length(tmp) != 2) {
					browser()
				}
			} else {
				tmp <- list.files(ff, full.names = TRUE)
				tmp <- tmp[str_detect(tmp, str_c(".*lane1", ss, ".*"))]
				#stopifnot(length(tmp) == 2)
				if (length(tmp) != 1) {
					if (ss == "MG48") {
						return(NA)
					}
					tmp <- list.files("/g/scb/zeller/karcher/PRISMA/data/WGS/230804_PRISMA_Batch5_P3", full.names = TRUE)
					tmp <- tmp[str_detect(tmp, str_c(".*lane1", ss, ".*"))]
					last_field <- str_split(tmp, "/")[[1]][length(str_split(tmp, "/")[[1]])]
					tmp <- c(
							str_c(tmp, "/", last_field, "_1.fastq.gz"),
							str_c(tmp, "/", last_field, "_2.fastq.gz")
						)
					return(tmp)					
				}				
				last_field <- str_split(tmp, "/")[[1]][length(str_split(tmp, "/")[[1]])]
				tmp <- c(
						str_c(tmp, "/", last_field, "_1.fastq.gz"),
						str_c(tmp, "/", last_field, "_2.fastq.gz")
					)
			}
			return(tmp)
		}
	)) %>%
	mutate(fastq_path_NA = map_lgl(fastq_names, \(x) any(is.na(x)))) %>%
	print(n=30) %>%
	select(
		PSN, sampleID, visit, fastq_names
	) %>%
	unnest() %>%
	write_tsv('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/MPA_Tox_PRISMA_sample_list.tsv')
