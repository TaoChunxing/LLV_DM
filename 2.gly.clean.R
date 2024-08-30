####3.入组后开始观察####
#####3.1血糖#####
#去除血糖缺失值
gly <- check.vlgroup.flw.rmVF %>% filter(!is.na(Gly))
#识别基线血糖异常的ID
ids.abnormalGLY <- gly %>% group_by(ID) %>% arrange(ID, TestDate) %>%
  slice(1) %>% filter(Gly >= 7 | Gly <= 2.8) %>% pull(ID)
#剔除基线血糖异常的ID,对随访的血糖异常值进行剔除
gly <- gly %>% filter(!ID %in% ids.abnormalGLY)
gly <- gly %>% filter(Gly>1&Gly<33)
gly <- gly %>% group_by(ID) %>% filter(n()>2)
summary(gly$Gly)
write.xlsx(gly, "2.0.gly.xlsx")
# hist(data.3.3.3$Gly, breaks = 120)
# boxplot(data.3.3.3$Gly)
# summary(data.3.3.3$Gly)
#IQR <- Q3 - Q1
#upper.IQR <- Q3 + 1.5 * IQR
#hist(data.hlme$Gly, breaks = 30)

#####3.2检测时间间隔####
#去除检测时间间隔小于等于30天的数据,去除检测时间间隔小于等于180天往后的数据
gly.testinterval.30.180 <- gly %>% group_by(ID) %>% arrange(ID, TestDate) %>%
  mutate(Testintervalday = c(0, as.numeric(diff(TestDate)))) %>% 
  filter(row_number() == 1 & Testintervalday == 0 | Testintervalday > 30) %>% filter(n()>2) %>% 
  mutate(Testintervalday = c(0, as.numeric(diff(TestDate)))) %>% 
  filter(row_number() == 1 & Testintervalday == 0 | Testintervalday > 30) %>% filter(n()>2) %>% 
  mutate(Testintervalday = c(0, as.numeric(diff(TestDate)))) %>% 
  mutate(FollowStatus = ifelse(Testintervalday > 180, "N", "Y")) %>% 
  mutate(cumsum_status = cumsum(FollowStatus == "N")) %>%
  filter(cumsum_status == 0) %>%
  select(-cumsum_status) %>% 
  filter(n()>2)
write.xlsx(gly.testinterval.30.180, "2.1.gly.testinterval.30.180.xlsx")

#####3.3观察终点#####
#统一观察7年截至，不管有没有发生血糖异常
gly.flw84ms <- gly.testinterval.30.180 %>% select(-17, -18) %>% 
  group_by(ID) %>% arrange(ID, TestDate) %>% 
  filter(Flwtime.m <= 84) %>% filter(n()>2)
write.xlsx(gly.flw84ms, "2.2.gly.flw.84ms.xlsx")

#观察7年后，连续出现2次血糖异常为1，未出现为0
gly.flw84ms.glygroup <- gly.flw84ms %>% group_by(ID) %>% 
  mutate(Glygroup = ifelse(any(Gly >= 7 & lead(Gly, 1) >= 7, na.rm = TRUE), "Hyperglycemia" , "Normal"))
gly.flw84ms.glygroup.unique <- gly.flw84ms.glygroup %>% group_by(ID) %>% distinct(ID, .keep_all = TRUE) %>% ungroup()
table(gly.flw84ms.glygroup.unique$Glygroup)
write.xlsx(gly.flw84ms.glygroup, "2.3.gly.flw84ms.glygroup.xlsx")
write.xlsx(gly.flw84ms.glygroup.unique, "2.4.gly.flw84ms.glygroup.unique.xlsx")
