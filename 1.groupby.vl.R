####1.定位到治疗半年以上的数据####
#####1.1找到治疗半年以上的病载记录，去除半年以前的记录#
check <- check %>% group_by(ID) %>% arrange(ID, TestDate) %>% 
  mutate(VLTfromART.m = round(interval(start=ARTinitialDate, end=VLDate) / months(1), digits = 2)) %>% 
  mutate(VLTfromART.6m = ifelse(!is.na(VL) & VLTfromART.m > 6, cumsum(!is.na(VL) & VLTfromART.m > 6),NA)) %>% 
  mutate(first_occurrence = if_else(any(VLTfromART.6m == 1), which(VLTfromART.6m == 1)[1], NA_integer_)) %>%
  filter(is.na(first_occurrence) | row_number() >= first_occurrence) %>% 
  select(-first_occurrence) %>% 
  filter(!all(is.na(VLTfromART.6m))) %>% 
  mutate(VLTfromART.6m = ifelse(is.na(VLTfromART.6m), -1, VLTfromART.6m)) 

#对治疗半年以上只有一次记录的ID，若不是VF则去除,因为两次记录才能判定LLV和VS
ids.1record <- check %>% group_by(ID) %>% filter(!any(VLTfromART.6m == 2))
ids.1record.removed <- ids.1record %>% group_by(ID) %>% filter(any(VLTfromART.6m == 1 & VL< 1000))
check <- anti_join(check, ids.1record.removed, by = "ID")


####2.VLgroup,分VS,LLV,VF组####
#####2.1提取出只有一条记录的ID，分组为VF#####
ids.1record.keep <- ids.1record %>% anti_join(ids.1record.removed, by = "ID") %>% 
  mutate(VLTfromART.6m = ifelse(VLTfromART.6m == -1, NA, VLTfromART.6m)) %>% 
  mutate(VLgroup = "VF")
#####2.2两条以上记录的，依据条件分VS,LLV,VF组#####
ids.more2record <- check %>% anti_join(ids.1record.keep, by = "ID")
ids.2record <- ids.more2record %>% group_by(ID) %>% filter(VLTfromART.6m %in% c(1, 2))
ids.2record <- ids.2record %>% group_by(ID) %>% 
  mutate(
    VLgroup = case_when(
      all(VL <= 50) ~ "VS",
      first(VL) >= 1000 |
        (first(VL) < 1000 & lead(VL) >= 1000) |
        (first(VL) >= 200 & lead(VL) >= 200) ~ "VF",
      TRUE ~ NA_character_
    )) %>%
  mutate(VLgroup = max(VLgroup, na.rm = TRUE))

ids.2record <- ids.2record %>% 
  mutate(VLgroup = ifelse(is.na(VLgroup)&(first(VL) > 50 & lead(VL) > 50), "LLV", 
                          ifelse(is.na(VLgroup), "Blips", VLgroup)))

table(ids.2record$VLgroup)

ids.2record.select <- ids.2record %>% select(1, 16) %>% distinct(ID, .keep_all = TRUE)
ids.1record.keep.select <- ids.1record.keep %>% select(1, 16) %>% distinct(ID, .keep_all = TRUE)
ids.vlgroup <- bind_rows(ids.2record.select, ids.1record.keep.select)
check.vlgroup <- check %>% left_join(ids.vlgroup, by = "ID")
table(check.vlgroup$VLgroup)

#####2.3计算入组时间#####
#第一次记录VL不是大于等于1000则删除，只有VF可以通过一次判断入组，其余都要观察两次，根据第二次值入组
check.vlgroup.flw <- check.vlgroup %>% group_by(ID) %>% arrange(ID, TestDate) %>%
  mutate(row_number = row_number()) %>%
  filter(!(VLgroup == "VF" & row_number == 1 & VL < 1000)) %>%
  filter(!(VLgroup == "VS" & row_number == 1)) %>%
  filter(!(VLgroup == "LLV" & row_number == 1)) %>%
  filter(!(VLgroup == "Blips" & row_number == 1)) %>% 
  select(-row_number)%>%
  mutate(StartflwDate = first(VLDate[!is.na(VLDate)]))%>% 
  mutate(Flwtime.m = round(interval(start=StartflwDate, end=TestDate) / months(1), digits = 2)) %>% 
  filter(Flwtime.m >= 0)
check.vlgroup.flw <- check.vlgroup.flw %>% select(-14,-15)
check.vlgroup.flw <- check.vlgroup.flw %>% select(1,2,3,5,15,everything())
check.vlgroup.flw <- check.vlgroup.flw %>% select(1:5,7:16,6)


######2.3.1去除入组日期晚于死亡日期的错误数据######
check.select <- check.vlgroup.flw %>% select(1,5) %>% group_by(ID) %>% distinct(ID, .keep_all = TRUE) %>% ungroup()
check.select <- check.select %>% left_join(death, by = "ID")
wrongID <- check.select %>% filter(DiedDate < StartflwDate) %>% pull(ID)
check.vlgroup.flw <- check.vlgroup.flw %>% filter(!ID %in% wrongID)
write.xlsx(check.vlgroup.flw, "1.0.check.vlgroup.flw.xlsx")

#####2.4去除入组为VF的数据#####
check.vlgroup.flw.rmVF <- check.vlgroup.flw %>% filter(VLgroup!="VF") %>% ungroup()
table(check.vlgroup.flw.rmVF$VLgroup)
write.xlsx(check.vlgroup.flw.rmVF, "1.1.check.vlgroup.flw.rmVF.xlsx")
str(check.vlgroup.flw.rmVF)

#####2.5提取每个ID的VL分组信息
vlgroup.info <- check.vlgroup.flw.rmVF %>% select(1,14) %>% distinct(ID, .keep_all = TRUE)
table(vlgroup.info$VLgroup)
