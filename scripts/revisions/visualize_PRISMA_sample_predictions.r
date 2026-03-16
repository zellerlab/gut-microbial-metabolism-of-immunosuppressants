library(tidyverse)
library(ggembl)

ll <- list.files('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/tmp/predictions_for_27_samples', full.names = TRUE)
ll <- ll[str_detect(ll, "Mofetil")]

predictions <- lapply(ll, \(x) {
	return(read_tsv(x) %>%
		mutate(name_tmp = x))
})

predictions <- enframe(predictions) %>%
	unnest() %>%
	select(-name) %>%
	select(-name_tmp) %>%
	mutate(sampleID = str_replace(sampleID, ".*MG", "MG_")) %>%
	print(n=35) %>%
	# Notice here that 4 samples are not having properly cleaned up sampleID entries yet.
	# I will fix these here now
	# note also it is no mistake that I remove the first '1' after lane - this is all correct
	mutate(sampleID = str_replace(sampleID, ".*lane1", "MG_")) %>%
	group_by(sampleID, target) %>%
	summarize(
		prediction = mean(prediction, na.rm = TRUE)
	) %>%
	inner_join( read_tsv('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/MPA_Tox_PRISMA_sample_list.tsv') %>%
	select(PSN, sampleID, visit) %>% 
	distinct() %>%
	relocate(PSN, visit, sampleID)) %>%
	inner_join(
		read_csv('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/MPA_Tox_df_for_Nic.csv'),
		by = c("PSN" = "v65_pat_id", 'visit' = "v62_visit_number")
	) %>%
	rename(`predicted MMF\nmetabolism` = prediction) %>%
	mutate(`gastroint. toxicity` = case_when(
		GIT_tox ~ "Yes",
		!GIT_tox ~ "No"
	)) %>%
	mutate(
		`gastroint. toxicity` = factor(`gastroint. toxicity`, levels = c("Yes", "No"))
	) 

p <- ggplot(predictions, aes(x = `predicted MMF\nmetabolism`, y = MPA_MMF_Ratio, color = `gastroint. toxicity`)) +
	geom_point() +
	scale_color_manual(
		values = c("Yes" = "red", "No" = "darkgreen"),
		name = "Gastrointestinal\ntoxicity"
	) +
	theme_presentation()

ggsave(
	str_c("/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/tmp/",
	'pred_mmf_vs_mpa_mmf_ratio.pdf'),
	width = 6, height = 4
)


