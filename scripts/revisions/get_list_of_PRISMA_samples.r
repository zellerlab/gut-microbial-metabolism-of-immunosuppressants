# load first few lines of /g/scb/zeller/karcher/PRISMA/scripts/PRISMA/2024_11_27_prepare_PRISMA_WGS.r
# then execute this 
# meta  %>% select(PSN, sampleID, visit, batch) %>% distinct() %>%  right_join(read_csv('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/MPA_Tox_df_for_Nic.csv') %>% rename(PSN = v65_pat_id, visit = v62_visit_number) %>% select(PSN, visit)) %>% write_tsv('/g/scb/zeller/karcher/gut-microbial-metabolism-of-critical-dose-immunosuppressants/data/MPA_Tox_PRISMA_sample_list.tsv')
