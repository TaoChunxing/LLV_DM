library(tidyverse)
library(openxlsx)
library(readxl)
#info数据清洗，PSM不允许存在NA
#提取进行轨迹分析的ID号
info <- read_excel("info.xlsx")
gly.flw84ms.glygroup.unique <- read_excel("2.4.gly.flw84ms.glygroup.unique.xlsx")
info.gly.age.vgroup <- gly.flw84ms.glygroup.unique %>% select(1,14,16,17)

ids <- gly.flw84ms.glygroup.unique %>% pull(ID)
info <- info %>% filter(ID %in% ids)
info <- info %>% left_join(info.gly.age.vgroup, by = "ID")
info <- info %>% ungroup %>% mutate_at(vars(ends_with("Date")), as.Date) #批量转化为日期型变量
write.xlsx(info,"5.0.info.xlsx")


####1.Age####
info.1 <- info %>%
  mutate(`ART initial age (years)`= case_when(
    ARTinitialAge >= 18 & ARTinitialAge <= 34 ~ "18-34",
    ARTinitialAge >= 35 & ARTinitialAge <= 49 ~ "35-49",
    ARTinitialAge >= 50 ~ ">=50"))
table(addNA(info.1$`ART initial age (years)`))
####2.Sex####
info.1 <- info.1 %>% mutate(Sex = ifelse(Gender == 1, 'Male', 'Female'))
table(addNA(info.1$Sex))

####3.MaritalStatus为4（不详）的人众数插补为2（已婚）####
info.1 <- info.1 %>% 
  mutate('Marital status' = ifelse(MaritalStatus == 4, 2, MaritalStatus))
info.1 <- info.1 %>% 
  mutate('Marital status' = case_when(
    `Marital status`==1~'Single',
    `Marital status`==2~'Marrried/Partnered',
    `Marital status`==3~'Divorced/Widowed'
  ))
table(addNA(info.1$`Marital status`))

####4.Ethnic为NA的众数插补为Minority#####
info.1$Ethnic.a = case_when(
  info.1$Ethnic =="汉族"~"Han",
    TRUE ~ "Minority"
  )
table(addNA(info.1$Ethnic.a))

####5.Education#####
info.1 <- info.1 %>% 
  mutate(`Educational attainment` =case_when(
    EducationLevel=="小学"|EducationLevel=="文盲" ~ "Primary school and below",
    EducationLevel=="初中" ~ "Junior school",
    EducationLevel=="大专及以上"|EducationLevel=="高中或中专" ~ "High school and above",
    TRUE ~ NA_character_
  ))
info.1 <- info.1 %>% 
  mutate(`Educational attainment` = ifelse(is.na(`Educational attainment`), "Junior school", `Educational attainment`))
table(addNA(info.1$`Educational attainment`))

####6.Occupation####
info.1 <- info.1 %>% 
  mutate(Occupation.a = case_when(
    grepl("民", Occupation)~"Farmer",
    !is.na(Occupation)~"Others",
    TRUE~NA_character_
  ))
info.1 <- info.1 %>% mutate(Occupation.a = ifelse(is.na(Occupation.a), "Farmer", Occupation.a))
table(addNA(info.1$Occupation.a))

####7.Transmission####
info.1 <- info.1 %>% mutate(`Transmission`= ifelse(InfectedRout == 2, "Heterosexual contact", "Other or unkonwn"))
table(addNA(info.1$`Transmission`))

###8.FFXNM####
info.1$`FFXNM use history` <- ifelse(info.1$UsedFFXNM==0,"No" , "Yes")
info.1 <- info.1 %>% mutate(`FFXNM use history` = ifelse(is.na(`FFXNM use history`), "No", `FFXNM use history`))
table(addNA(info.1$`FFXNM use history`))

####9.WHOStage为NA的人众数插补为1####
info.1 <- info.1 %>% 
  mutate(`Baseline WHO HIV stage` = ifelse(is.na(WHOStage), "01", WHOStage))
info.1 <- info.1 %>%
  mutate(`Baseline WHO HIV stage` = case_when(
    `Baseline WHO HIV stage` == "01" ~ "Stage I",
    `Baseline WHO HIV stage` == "02" ~ "Stage II",
    `Baseline WHO HIV stage` == "03" ~ "Stage III",
    `Baseline WHO HIV stage` == "04" ~ "Stage IV"))
table(addNA(info.1$`Baseline WHO HIV stage`))

####10.Regimen####
info.1$InitialRegimen <- gsub("双汰芝","AZT+3TC", info.1$InitialRegimen)
info.1$InitialRegimen <- gsub("克力芝","LPV/r", info.1$InitialRegimen)
info.1$InitialRegimen <- gsub("/", "+", info.1$InitialRegimen)
NRTIs <- c("3TC","D4T", "TDF", "TAF", "FTC", "AZT", "DDI")
NNRTIs <- c("EFV", "RPV", "DOR","NVP")
PIs <- c("LPV", "r", "DRV", "c", "NFV", "ANV")
INSTIs <- c("DTG", "RAL", "BIC", "EVG")

f.Regimenbased <- function(regimen) {
  drugs <- unlist(str_split(regimen, "\\+"))
  if (any(drugs %in% NNRTIs)) {
    return("NNRTIs-based")
  } else if (any(drugs %in% PIs)) {
    return("PIs-based")
  } else if (any(drugs %in% INSTIs)) {
    return("INSTIs-based")
  } else {
    return("Other/combinations")
  }
}

info.1 <- info.1 %>% mutate(`ART initial regimen` = sapply(InitialRegimen, f.Regimenbased))
table(addNA(info.1$`ART initial regimen`))

info.1 <- info.1 %>%
  mutate(
    `ART initial regimen.a` = case_when(
      str_detect(InitialRegimen, "EFV") ~ "EFV-based",
      str_detect(InitialRegimen, "RPV") ~ "RPV-based",
      str_detect(InitialRegimen, "DOR") ~ "DOR-based",
      str_detect(InitialRegimen, "NVP") ~ "NVP-based",
      TRUE ~ `ART initial regimen`
    )
  )
table(addNA(info.1$`ART initial regimen.a`))


####11.CD4####
info.1 <- info.1 %>% mutate_at(vars(ends_with("CD4")),as.numeric)#批量转化为数值型变量

info.1 <- info.1 %>%
  mutate(`Baseline CD4 count (cells/µL)` = case_when(
    InitialCD4 < 200 ~ "<200",
    InitialCD4 >= 200 & InitialCD4 < 350 ~ "200-350",
    InitialCD4 >= 350 & InitialCD4 < 500 ~ "350-500",
    InitialCD4 >= 500 ~ ">=500",
    TRUE ~ NA_character_  
  )) %>%
  mutate(`Baseline CD4 count (cells/µL)`=ifelse(
    is.na(`Baseline CD4 count (cells/µL)`), "<200", `Baseline CD4 count (cells/µL)`)) 
table(addNA(info.1$`Baseline CD4 count (cells/µL)`))

info.1 <- info.1 %>%
  mutate(`Latest CD4 count (cells/µL)` = case_when(
    RecentCD4 < 200 ~ "<200",
    RecentCD4 >= 200 & RecentCD4 < 350 ~ "200-350",
    RecentCD4 >= 350 & RecentCD4 < 500 ~ "350-500",
    RecentCD4 >= 500 ~ ">=500",
    TRUE ~ NA_character_    
  )) %>%
  mutate(`Latest CD4 count (cells/µL)`=ifelse(
    is.na(`Latest CD4 count (cells/µL)`), "<200", `Latest CD4 count (cells/µL)`)) 
table(addNA(info.1$`Latest CD4 count (cells/µL)`))

####12.VL####
info.1 <- info.1 %>% mutate_at(vars(ends_with("VL")),as.numeric)#批量转化为数值型变量
info.1 <- info.1 %>%
  mutate(`Baseline viral load (copies/mL)` = case_when(
    InitialVL < 50 ~ "<50",
    InitialVL >= 50 & InitialVL < 1000 ~ "50-1000",
    InitialVL >= 1000  ~ ">=1000",
    TRUE ~ NA_character_    
  )) 
table(addNA(info.1$`Baseline viral load (copies/mL)`))

info.1 <- info.1 %>%
  mutate(`Latest viral load (copies/mL)` = case_when(
    RecentVL < 50 ~ "<50",
    RecentVL >= 50 & RecentVL < 1000 ~ "50-1000",
    RecentVL >= 1000  ~ ">=1000",
    TRUE ~ NA_character_
  ))
table(addNA(info.1$`Latest viral load (copies/mL)`))

####12.VLgroup 012####
info.1 <- info.1 %>% 
  mutate(`Enrollment VL group (VS/LV)` = case_when(
    VLgroup == "VS" ~ "VS",
    VLgroup == "Blips"| VLgroup == "LLV" ~ "LV"))
table(info.1$`Enrollment VL group (VS/LV)`)
info.1 <- info.1 %>% rename(`Enrollment VL group (VS/Blips/LLV)` = VLgroup)
table(info.1$`Enrollment VL group (VS/Blips/LLV)`)


####13.血糖####
info.1 <- info.1 %>% rename(`Blood glucose group` = Glygroup)
info.1 <- info.1 %>% mutate(`Blood glucose group`= ifelse(
  `Blood glucose group` == "Hyperglycemia", "Diabetes mellitus", `Blood glucose group`))
table(info.1$`Blood glucose group`)
write.xlsx(info.1, "5.1.info.xlsx")


####15.简化数据框####
info.2 <- info.1 %>% select(1, 20,22:38)
info.2 <-info.2 %>% select(1, 2, 19, everything())
info.2 <-info.2 %>% select(-14)
table(info.2$`Enrollment VL group (VS/LV)`, info.2$`Blood glucose group`)
write.xlsx(info.2, "5.2.info.xlsx")


#将变量名中的空格.替换为下划线
colnames(info.2) <- str_replace_all(colnames(info.2), " ", "_")
colnames(info.2) <- str_replace_all(colnames(info.2), "\\.", "_")
write.xlsx(info.2, "5.3.info.xlsx")
