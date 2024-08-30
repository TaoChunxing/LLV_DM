library(tidyverse)
library(readxl)


####1.整理变量####
surv.glygroup.last <- read_excel("4.0.surv.glygroup.last.xlsx")
info.2 <- read_excel("5.3.info.xlsx")
data.multicox <- surv.glygroup.last %>% select(1,4) %>% left_join(info.2, by = "ID")

#将自变量改为因子型变量
data.multicox <- data.multicox %>% mutate_at(vars(-ID, -Flwtime.m), as.factor)
#重新排列因子的level
data.multicox <- data.multicox %>%
  mutate(`ART_initial_age_(years)` = fct_relevel(`ART_initial_age_(years)`, '18-34', '35-49', '>=50'),
         Sex = fct_relevel(Sex, 'Male', 'Female'),
         FFXNM_use_history = fct_relevel(FFXNM_use_history, 'No', 'Yes'),
         Marital_status = fct_relevel(Marital_status, 'Marrried/Partnered', 'Divorced/Widowed', 'Single'),
         Educational_attainment = fct_relevel(Educational_attainment, 'Primary school and below', 'Junior school', 'High school and above'),
         ART_initial_regimen_a = fct_relevel(ART_initial_regimen_a, 'EFV-based','NVP-based', 'PIs-based', 'INSTIs-based', 'Other/combinations'),
         `Baseline_CD4_count_(cells/µL)` = fct_relevel(`Baseline_CD4_count_(cells/µL)`,  '<200', '200-350', '350-500', '>=500'),
         `Baseline_viral_load_(copies/mL)` = fct_relevel(`Baseline_viral_load_(copies/mL)`, '<50', '50-1000', '>=1000'),
         `Latest_CD4_count_(cells/µL)` = fct_relevel(`Latest_CD4_count_(cells/µL)`,  '<200', '200-350', '350-500', '>=500'),
         `Enrollment_VL_group_(VS/Blips/LLV)` = fct_relevel(`Enrollment_VL_group_(VS/Blips/LLV)`, 'VS','Blips', 'LLV'),
         `Latest_viral_load_(copies/mL)` = fct_relevel(`Latest_viral_load_(copies/mL)`, '<50', '50-1000', '>=1000'),
         `Enrollment_VL_group_(VS/LV)` = fct_relevel(`Enrollment_VL_group_(VS/LV)`, 'VS', 'LV'))
#cox回归的应变量需为数值型
data.multicox <- data.multicox %>% 
  mutate(Diabetes_mellitus=ifelse(Blood_glucose_group=="Diabetes mellitus", 1,0))

str(data.multicox)
table(data.multicox$ART_initial_regimen_b)
table(data.multicox$Diabetes_mellitus)
data.multicox <- data.multicox %>% rename(Time = Flwtime.m)
summary(data.multicox$Time)
write.xlsx(data.multicox, "6.0.data.multicox.xlsx")



####2.多因素cox回归####
library(survminer)
library(survival)
library(autoReg)#该函数的变量名不能有括号
str(data.multicox)
#cox回归模型构建
#婚姻状况不满足等比例风险假设
ccoxmod <- coxph(Surv(Time, Diabetes_mellitus) ~ ART_initial_age + Sex +
                  Ethnic_a + Educational_attainment + Occupation_a + Transmission +
                  Baseline_WHO_HIV_stage + ART_initial_regimen_b + FFXNM_use_history +
                  Baseline_CD4_count + Latest_CD4_count +
                  Baseline_viral_load+ Latest_viral_load +
                  Enrollment_VL_group,
                  data.multicox)
summary(ccoxmod)
# 使用 Schoenfeld 残差法检验 PHA,P>0.05满足等比例风险假设
ccox.zph <- cox.zph(ccoxmod)
print(ccox.zph)
#将回归结果生成表格
ft3<-autoReg(ccoxmod,uni=TRUE,threshold=0.1, final= TRUE)
myft(ft3)
#输出表格
library(rrtable)
table2docx(ft3)



####3.PSM匹配####
library(MatchIt)
#Nearest Neighbor Matching
matchlist <- matchit(`Enrollment_VL_group_a` ~ ART_initial_age+Sex,
                     data = data.multicox,
                     method = "nearest",
                     distance = "glm",
                     caliper = 0.01,
                     ratio = 1,
                     replace = F)
summary(matchlist)
#3.2提取匹配后的数据
matchdata <- match.data(matchlist,
                        group = "all",
                        distance = "distance",
                        weights = "weights",
                        subclass = "subclass",
                        data = NULL,
                        include.s.weights = TRUE,
                        drop.unmatched = TRUE)
str(matchdata)
write.xlsx(matchdata,"6.1.matched.age.sex.xlsx")

#检验VLgroup与Regimen之间是否还有统计学差异（无）
#chisq.test(matchdata$`Enrollment_VL_group_(VS/LV)`, matchdata$`ART_Initial_regimen_a`)
#Warning message:
#In chisq.test(table) : Chi-squared近似算法有可能不准,改用Fisher 精确检验
#fisher.test(matchdata$`Enrollment_VL_group_(VS/LV)`, matchdata$`ART_initial_regimen_a`)
# #检验VLgroup与Gender之间是否还有统计学差异（无）
# chisq.test(matchdata$Enrollment_VL_group_a, matchdata$Sex)
# #检验VLgroup与Age之间是否还有统计学差异（无）
# chisq.test(matchdata$Enrollment_VL_group_a, matchdata$ART_initial_age)
# #检验VLgroup与Glygroup之间是否有统计学差异（有）
# chisq.test(matchdata$Enrollment_VL_group_a, matchdata$Blood_glucose_group)



####4.匹配后ID的多因素cox模型####
matchID <- matchdata %>% pull(ID)
data.multicox.match <- data.multicox %>% filter(ID %in% matchID)
str(data.multicox.match)
write.xlsx(data.multicox.match, "6.2.data.multicox.match.xlsx")

#匹配后的数据多因素cox
#Baseline_viral_load不符合等比例风险假设
matcoxmod <- coxph(Surv(Time, Diabetes_mellitus) ~ ART_initial_age + Sex+
                     Ethnic_a + Educational_attainment +
                     Occupation_a + Transmission + 
                     Baseline_WHO_HIV_stage + ART_initial_regimen_b + FFXNM_use_history +
                     Baseline_CD4_count + Latest_CD4_count +
                     Latest_viral_load +
                     Enrollment_VL_group_a,
                     data.multicox.match)
summary(matcoxmod)
#使用Schoenfeld 残差法检验 PHA,P>0.05满足等比例风险假设
matcox.zph <- cox.zph(matcoxmod)
print(matcox.zph) 

#输出结果
ft3<-autoReg(matcoxmod,uni=TRUE,threshold=0.1, final= TRUE)
myft(ft3)
table2docx(ft3)




####5.基线统计表####
library(autoReg)
str(data.baseline)
data.baseline <- data.multicox %>% select(-1, -3, -5, -15,-20)
data.baseline <- data.baseline %>% select(1:11,16,everything())
write.xlsx(data.baseline, "6.3.data.baseline.xlsx")

baseline.table1 = gaze(Enrollment_VL_group_a~.,data = data.baseline)
table2docx(baseline.table1)

baseline.table2 = gaze(Sex~.,data = data.baseline)
table2docx(baseline.table2)

#匹配后的数据基线特征
str(data.multicox.match)
data.baseline.match <- data.multicox.match %>% select(-1, -3, -5, -15,-20)
data.baseline.match <- data.baseline.match %>% select(1:11,16,everything())
write.xlsx(data.baseline, "6.4.data.baseline.match.xlsx")

baseline.match.table1 = gaze(Enrollment_VL_group_a~.,data = data.baseline.match)
table2docx(baseline.match.table1)

baseline.match.table2 = gaze(Enrollment_VL_group_a+Sex~.,data = data.baseline.match)
table2docx(baseline.match.table2)
baseline.match.table2



####6.分层Cox回归####
data.male <- data.multicox %>% filter(Sex == 'Male')
ccoxmod.male <-  coxph(Surv(Time, Diabetes_mellitus) ~ ART_initial_age +
                         Ethnic_a + Educational_attainment + Occupation_a + Transmission +
                         Baseline_WHO_HIV_stage + ART_initial_regimen_b + FFXNM_use_history +
                         Baseline_CD4_count + Latest_CD4_count +
                         Baseline_viral_load+ Latest_viral_load +
                         Enrollment_VL_group,
                       data.male)
summary(ccoxmod.male)
ft3<-autoReg(ccoxmod.male,uni=TRUE,threshold=0.1, final= TRUE)
myft(ft3)
table2docx(ft3)

data.female <- data.multicox %>% filter(Sex == 'Female')
ccoxmod.female <- coxph(Surv(Time, Diabetes_mellitus) ~ ART_initial_age + 
                          Ethnic_a + Educational_attainment + Occupation_a + Transmission +
                          Baseline_WHO_HIV_stage + ART_initial_regimen_b + FFXNM_use_history +
                          Baseline_CD4_count + Latest_CD4_count +
                          Baseline_viral_load+ Latest_viral_load +
                          Enrollment_VL_group,
                        data.female)
summary(ccoxmod.female)
ft3<-autoReg(ccoxmod.female,uni=TRUE,threshold=0.1, final= TRUE)
myft(ft3)
table2docx(ft3)


data.age1 <- data.multicox %>% filter(ART_initial_age == '18-34')
data.age2 <- data.multicox %>% filter(ART_initial_age == '35-49')
data.age3 <- data.multicox %>% filter(ART_initial_age == '>=50')
coxmod.age1 <- coxph(Surv(Time, Diabetes_mellitus) ~ Sex + 
                       Ethnic_a + Educational_attainment + Occupation_a + Transmission +
                       Baseline_WHO_HIV_stage + ART_initial_regimen_b + FFXNM_use_history +
                       Baseline_CD4_count + Latest_CD4_count +
                       Baseline_viral_load+ Latest_viral_load +
                       Enrollment_VL_group,
                     data.age1)
summary(coxmod.age1)
ft3<-autoReg(coxmod.age1,uni=TRUE,threshold=0.1, final= TRUE)
myft(ft3)
table2docx(ft3)
ccox.zph <- cox.zph(coxmod)
print(ccox.zph)
coxmod.age2 <- coxph(Surv(Time, Diabetes_mellitus) ~ 
                       Enrollment_VL_group,
                     data.age2)
summary(coxmod.age2)
coxmod.age3 <- coxph(Surv(Time, Diabetes_mellitus) ~ 
                       Enrollment_VL_group,
                     data.age3)
summary(unicoxmod.age3)



data.match.male <- data.multicox.match %>% filter(Sex == 'Male')
unicoxmod.match.male <- coxph(Surv(Time, Diabetes_mellitus) ~ 
                                Enrollment_VL_group_a,
                              data.match.male )
summary(unicoxmod.match.male)

data.match.female <- data.multicox.match %>% filter(Sex == 'Female')
unicoxmod.match.female <- coxph(Surv(Time, Diabetes_mellitus) ~ 
                                  Enrollment_VL_group_a,
                                data.match.female)
summary(unicoxmod.match.female)



match.data.age1 <- data.multicox.match %>% filter(ART_initial_age == '18-34')
match.data.age2 <- data.multicox.match %>% filter(ART_initial_age == '35-49')
match.data.age3 <- data.multicox.match %>% filter(ART_initial_age == '>=50')
unicoxmod.match.age1 <- coxph(Surv(Time, Diabetes_mellitus) ~ 
                                Enrollment_VL_group_a,
                              match.data.age1)
summary(unicoxmod.match.age1)
unicoxmod.match.age2 <- coxph(Surv(Time, Diabetes_mellitus) ~ 
                                Enrollment_VL_group_a,
                              match.data.age2)
summary(unicoxmod.match.age2)
unicoxmod.match.age3 <- coxph(Surv(Time, Diabetes_mellitus) ~ 
                                Enrollment_VL_group_a,
                              match.data.age3)
summary(unicoxmod.match.age3)
