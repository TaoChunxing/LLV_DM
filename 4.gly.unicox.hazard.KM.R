library(survival)
library(survminer)

####1.对hlme截至Flwtime第二次发生血糖异常####
surv.gly0 <- data.hlme %>% filter(Glygroup=="Normal")
surv.gly1 <- data.hlme %>% filter(Glygroup=="Hyperglycemia") %>% 
  group_by(ID) %>% arrange(ID,Flwtime.m) %>% 
  mutate(flag = cumsum(lag(Gly >= 7, default = FALSE) & Gly >= 7)) %>%
  filter(row_number() <= which(flag == 1)[1]) %>% 
  select(-flag)

surv.glygroup <- bind_rows(surv.gly0, surv.gly1)
surv.glygroup.last <- surv.glygroup %>% group_by(ID) %>% arrange(ID, Flwtime.m) %>% slice_tail(n = 1) %>% ungroup()

surv.glygroup.last$VLgroup <- factor(surv.glygroup.last$VLgroup,levels = c("VS", "Blips", "LLV"))
str(surv.glygroup.last)
table(surv.glygroup.last$VLgroup)
table(surv.glygroup.last$Glygroup)

write.xlsx(surv.glygroup.last, "4.0.surv.glygroup.last.xlsx")


data.surv <- surv.glygroup.last %>% 
  mutate(VLgroup.a=case_when(
  VLgroup=="VS"~0,
  VLgroup=="Blips"~1,
  VLgroup=="LLV"~2
))

data.surv <- surv.glygroup.last %>% 
  mutate(Glygroup.a=case_when(
  Glygroup=="Normal"~0,
  Glygroup=="Hyperglycemia"~1
))
str(data.surv)
data.surv <- data.surv %>% mutate(VLgroup = fct_relevel(VLgroup, 'VS', 'Blips', 'LLV'))
write.xlsx(data.surv, "4.1.data.surv.xlsx")
# ####2.单因素cox回归获取HR值####
# table(data.surv$Glygroup)
mcoxmod <- coxph(Surv(Flwtime.m, Glygroup.a) ~ VLgroup, data = data.surv)
summary(mcoxmod)
# #计算HR(95%CI)
# summary.cox <- summary(cox)
# vs.hr <- summary.cox$coef[1, "exp(coef)"]
# vs.hrlower <- summary.cox$conf.int[1, "lower .95"]
# vs.hrupper <- summary.cox$conf.int[1, "upper .95"]
# vs.hr.95ci <- paste0("HR = ", round(vs.hr, 2), " (95% CI: ", round(vs.hrlower, 2), "-", round(vs.hrupper, 2), ")")
# vs.hr.95ci
# 
# llv.hr <- summary.cox$coef[2, "exp(coef)"]
# llv.hrlower <- summary.cox$conf.int[2, "lower .95"]
# llv.hrupper <- summary.cox$conf.int[2, "upper .95"]
# llv.hr.95ci <- paste0("HR = ", round(llv.hr, 2), " (95% CI: ", round(llv.hrlower, 2), "-", round(llv.hrupper, 2), ")")
# llv.hr.95ci
# write_rds(cox,"cox.rds")

####3.拟合血糖风险函数####
fit <- survfit(Surv(Flwtime.m,Glygroup.a) ~ VLgroup,
               data = data.surv)
write_rds(fit,"survfit.rds")

#画血糖风险函数图（KM风险曲线图）
ggsurvplot(
  fit,
  data = data.surv,
  fun = "cumhaz",#将生存曲线转化为风险曲线
  linetype = 1, # 根据分层更改线型c(0,1) or  c("solid", "dashed") or "strata"
  #surv.median.line = "hv", # 同时显示垂直和水平参考线 即增加中位生存时间 可选 "none"、"hv"、"h"、"v"
  palette = "lancet" ,#定义颜色 可选调色板有 "hue" "grey","npg","aaas","lancet","jco", "ucscgb","uchicago","simpsons"和"rickandmorty".
  xlab = "Follow Time (months)",
  ylab = "Cumulative Hazard Rate of Diabetes Mellitus",
  title = "",
  legend = "right", # 指定图例位置 "top"(默认),"left","right","none","bottom"
  legend.title = "",
  legend.labs = c("VS", "Blips", "LLV"),
  xlim = c(0,85),
  ylim = c(0,0.8),
  break.x.by = 12,
  break.y.by = .05,
  axes.offset = F, # 逻辑词，默认为TRUE。FALSE则生存曲线图的坐标轴从原点开始
  #conf.int = TRUE,#增加置信区间
  #conf.int.alpha = .3,# 数值，指定置信区间填充颜色的透明度； # 数值在0-1之间，0为完全透明，1为不透明
  #pval = TRUE, #log 秩检验
  #pval.size = 5,# 指定p值文本大小的数字，默认为 5
  pval.coord = c(36,0.05),# 长度为2的数字向量，指定p值位置x、y，如pval.coord=c(x,y)
  censor = T, # 逻辑词，默认为TRUE，在图上绘制删失点。
  censor.shape = 3, # 数值或字符，用于指定删失点的形状；默认为"+"(3), 可选"|"(124)
  censor.size = 1.5,# 指定删失点形状的大小，默认为4.5
  risk.table = "absolute", #"absolute"、"percentage"、"abs_pct"", #绝对人数、百分比和危险之中
  risk.table.col = "strata", #按组更改风险表颜色
  risk.table.y.text.col = TRUE, #颜色风险表文本注释（按层）
  risk.table.y.text = FALSE, #在风险表图例中的文本注释中显示条形而不是名称
  risk.table.height = 0.2,
  risk.table.title="Number at the time",
  fontsize=4,#风险表字体
  ggtheme = theme_bw()+ theme(
    panel.grid = element_blank(), # 去除网格线
    axis.title.x = element_text(face = "bold", size = 12), # X轴标签加粗
    axis.title.y = element_text(face = "bold", size = 12)  # Y轴标签加粗
  ))

