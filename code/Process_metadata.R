## only working with metadata

meta_sample <- read_delim("../metadata/TCMA/MetaLight_Sample_PatAge_V2.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) 
meta_case <- read_delim("../metadata/TCMA/MetaLight_Case_SampleInfo.txt", delim = "\t", escape_double = F, trim_ws = T)
meta_case_org <- read_delim("../metadata/TCMA/metadata.TCMA.case.txt", delim = "\t", escape_double = F, trim_ws = T)

## join meta case

meta_case_2 <- left_join(x=meta_case, y=meta_case_org[,c("bcr_patient_barcode", "Stage")], by = c("bcr_patient_barcode"))
meta_sample_2 <- left_join(x=meta_sample, y=meta_case_2[, c("bcr_patient_barcode","gender", "race", "history_of_neoadjuvant_treatment", "radiation_therapy", "tss_site", "Stage")], by= "bcr_patient_barcode")

write_delim(meta_sample_2, "../metadata/TCMA/MetaLight_Sample_PatAge_V3.txt", delim = "\t", col_names = TRUE, quote = "none")
