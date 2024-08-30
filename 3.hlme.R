
library(lcmm)
####5.匹配前的血糖轨迹——不分组####
#####5.1简化数据框(id,time,var-gly,fac-vl)，生成Subject#####
data.hlme <- gly.flw84ms.glygroup %>% select( 1, 13, 14, 15, 17) %>% 
  group_by(ID) %>% arrange(ID, Flwtime.m)%>% 
  mutate(Subject = cur_group_id())
write.xlsx(data.hlme, "3.0.data.hlme.xlsx")

#####5.2计算模型#####
model.1 <- hlme(fixed = Gly ~ 1 + Flwtime.m + I(Flwtime.m^2),
                random = ~ 1 + Flwtime.m,
                ng = 1, 
                data = data.hlme, 
                subject = "Subject")
lin <- c(unmatch.model.1$ng, unmatch.model.1$BIC)
model <- list(model.1)
for (i in 2:4) {
  mi <- hlme(fixed = Gly ~ 1 + Flwtime.m + I(Flwtime.m^2),
             mixture = ~ 1 + Flwtime.m + I(Flwtime.m^2),
             random = ~ 1 + Flwtime.m,
             ng = i,
             nwg = TRUE, 
             B = model.1,
             idiag = FALSE, 
             maxiter = 50,
             data = data.hlme, 
             subject = "Subject")
  lin <- rbind(lin, c(i, mi$BIC))
  model[[i]] <- mi
}
modelout <- knitr::kable(lin, col.names = c("k", "BIC"), row.names = FALSE, align = "c")
modelout

#####5.3绘图#####
datnew  <- data.frame(Flwtime.m = seq(0,84, length = 1000))
plotpred <- predictY(model.2, datnew, var.time ="Flwtime.m")
windows(width=10, height=8)
plot(plotpred, lty = 1, lwd = 3, marg=FALSE, shades=T,
     xlab="Follow Time (Months)", ylab="Fasting Blood Glucose Level (mmol/L)", 
     xaxt='n', # 关闭自动的 x 轴刻度
     legend.loc = "topleft", cex=0.75)
axis(1, at = seq(0, 84, by = 12), labels = seq(0, 84, by = 12))

#####5.4提取不同轨迹组的对象#####
class <- model.2$pprob[,1:2]
#####5.5合并hlme与class#####
data.hlme.class <- data.hlme %>% left_join(class, by = "Subject")
data.hlme.class.unique <- data.hlme.class %>% distinct(ID, .keep_all = TRUE)
write.xlsx(data.hlme.class, "3.1.data.hlme.class.xlsx")
write.xlsx(data.hlme.class.unique, "3.2.data.hlme.class.unique.xlsx")

#####5.6制作VLgroup与class的交叉表#####
table <- table(data.hlme.class.unique$VLgroup, data.hlme.class.unique$Class)
percent <- prop.table(table, margin = 1) * 100
count.percent <- table
for (i in 1:nrow(table)) {
  for (j in 1:ncol(table)) {
    count.percent[i, j] <- paste0(table[i, j], " (", round(percent[i, j], 2), "%)")
  }
}
count.percent <- as.data.frame.matrix(count.percent)
rowtotals <- rowSums(table)
coltotals <- colSums(table)
count.percent$Total <- rowtotals
count.percent <- rbind(count.percent, Total = c(coltotals, sum(coltotals)))
#count.percent <- count.percent %>% select(`1`, `0`,`Total`)
count.percent
chisq.test(table)
fisher.test(table)




####6.匹配前的血糖轨迹——VLgroup分组####
df.vs.hlme <- data.hlme %>% filter(VLgroup=="VS")
df.blips.hlme <- data.hlme %>% filter(VLgroup=="Blips")
df.llv.hlme <- data.hlme %>% filter(VLgroup=="LLV")
df.lv.hlme <- data.hlme %>% filter(VLgroup == "Blips" | VLgroup == "LLV")


#####6.1 llv####
model.llv.1 <- hlme(fixed = Gly ~ 1 + Flwtime.m + I(Flwtime.m^2),
                    random = ~ 1 + Flwtime.m,
                    ng = 1, 
                    data = df.llv.hlme, 
                    subject = "Subject")
lin1 <- c(model.llv.1$ng, model.llv.1$BIC)
model.llv <- list(model.llv.1)
for (i in 2:4) {
  mi <- hlme(fixed = Gly ~ 1 + Flwtime.m + I(Flwtime.m^2),
             mixture = ~ 1 + Flwtime.m + I(Flwtime.m^2),
             random = ~ 1 + Flwtime.m,
             ng = i,
             nwg = TRUE, 
             B = model.llv.1,
             idiag = FALSE, 
             maxiter = 50,
             data = df.llv.hlme, 
             subject = "Subject")
  lin1 <- rbind(lin1, c(i, mi$BIC))
  model.llv[[i]] <- mi
}
modelout <- knitr::kable(lin1, col.names = c("k", "BIC"), row.names = FALSE, align = "c")
modelout


#####6.2 blips####
model.blips.1 <- hlme(fixed = Gly ~ 1 + Flwtime.m + I(Flwtime.m^2),
                      random = ~ 1 + Flwtime.m,
                      ng = 1, 
                      data = df.blips.hlme, 
                      subject = "Subject")
lin2 <- c(model.blips.1$ng, model.blips.1$BIC)
model.blips <- list(model.blips.1)
for (i in 2:4) {
  mi <- hlme(fixed = Gly ~ 1 + Flwtime.m + I(Flwtime.m^2),
             mixture = ~ 1 + Flwtime.m + I(Flwtime.m^2),
             random = ~ 1 + Flwtime.m,
             ng = i,
             nwg = TRUE, 
             B = model.blips.1,
             idiag = FALSE, 
             maxiter = 50,
             data = df.blips.hlme, 
             subject = "Subject")
  lin2 <- rbind(lin2, c(i, mi$BIC))
  model.blips[[i]] <- mi
}
modelout <- knitr::kable(lin2, col.names = c("k", "BIC"), row.names = FALSE, align = "c")
modelout 


#####6.3 vs####
model.vs.1 <- hlme(fixed = Gly ~ 1 + Flwtime.m + I(Flwtime.m^2),
                   random = ~ 1 + Flwtime.m,
                   ng = 1, 
                   data = df.vs.hlme, 
                   subject = "Subject")
lin3 <- c(model.vs.1$ng, model.vs.1$BIC)
model.vs <- list(model.vs.1)
for (i in 2:4) {
  mi <- hlme(fixed = Gly ~ 1 + Flwtime.m + I(Flwtime.m^2),
             mixture = ~ 1 + Flwtime.m + I(Flwtime.m^2),
             random = ~ 1 + Flwtime.m,
             ng = i,
             nwg = TRUE, 
             B = model.vs.1,
             idiag = FALSE, 
             maxiter = 50,
             data = df.vs.hlme, 
             subject = "Subject")
  lin3 <- rbind(lin3, c(i, mi$BIC))
  model.vs[[i]] <- mi
}
modelout <- knitr::kable(lin3, col.names = c("k", "BIC"), row.names = FALSE, align = "c")
modelout
postprob(model.vs[[2]])

#####6.4 将预测结果转化为数据框#####
# 预测
datnew <- data.frame(Flwtime.m = seq(0, 84, length = 1000))
plotpred.vs <- predictY(mod.vs.2, datnew, var.time ="Flwtime.m")
plotpred.blips <- predictY(mod.blips.2, datnew, var.time ="Flwtime.m")
plotpred.llv <- predictY(mod.llv.2, datnew, var.time ="Flwtime.m")
# 提取预测结果为数据框的函数
extract.predictions <- function(prediction, group) {
  data.frame(Flwtime.m = prediction$times$Flwtime.m,
             Class1 = prediction$pred[, 1],
             Class2 = prediction$pred[, 2],
             Group = group)
}
# 提取并转换预测结果
df.plotpred.vs <- extract.predictions(plotpred.vs, "VS")
df.plotpred.blips <- extract.predictions(plotpred.blips, "Blips")
df.plotpred.llv <- extract.predictions(plotpred.llv, "LLV")
# 合并数据框
df.plotpred<- bind_rows(df.plotpred.vs, df.plotpred.blips, df.plotpred.llv)
# 修改数据框中的列名
  colnames(df.plotpred)[which(names(df.plotpred) == "Class1")] <- "Rapid increase"
colnames(df.plotpred)[which(names(df.plotpred) == "Class2")] <- "Stable"


#####6.5 绘图
library(ggplot2)
library(ggsci)#里面有lancet主题颜色

windows(width=10, height=8)
ggplot(df.plotpred, aes(x = Flwtime.m)) + 
  geom_line(aes(y = `Rapid increase`, color = Group, linetype = "Rapid increase"), linewidth = 1) + 
  geom_line(aes(y = `Stable`, color = Group, linetype = "Stable"), linewidth = 1) + 
  labs(x = "Follow Time (Months)", y = "Fasting Blood Glucose Level (mmol/L)", 
       title = NULL) + 
  scale_color_lancet() + #使用lancet主题颜色
  scale_linetype_manual(values = c("Rapid increase" = "solid", "Stable" = "dashed"))+
  scale_x_continuous(breaks = seq(0, 84, by = 12)) +  # 设置 X 轴刻度范围和间隔
  scale_y_continuous(breaks = seq(0, 14, by = 2))+ 
  theme_minimal() +
  theme(
    # panel.grid = element_blank(),  # 去除网格背景
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),  # 添加边框
    legend.position = "right",
    legend.box = "vertical", 
    axis.title.x = element_text(face = "bold"), # X轴标签字体加粗
    axis.title.y = element_text(face = "bold"),  # Y轴标签字体加粗
    axis.text.x = element_text(face = "bold"),   # X轴刻度字体加粗
    axis.text.y = element_text(face = "bold")    # Y轴刻度字体加粗
  ) + 
  guides(linetype = guide_legend(title = "FBG trajectory group"), color = guide_legend(title = "Enrollment VL group"))



####6. 统计每一年估计组的人数####
llv.class <- mod.llv.2$pprob[,1:2]
blips.class <- mod.blips.2$pprob[,1:2]
vs.class <- mod.vs.2$pprob[,1:2]
class <- bind_rows(llv.class, blips.class, vs.class)
hlme.class <- data.hlme %>% left_join(class, by = "Subject")

hlme.class.1 <- hlme.class %>%
  group_by(ID) %>%
  mutate(
    m12 = ifelse(any(Flwtime.m > 12), "Y", "N"),
    m24 = ifelse(any(Flwtime.m > 24), "Y", "N"),
    m36 = ifelse(any(Flwtime.m > 36), "Y", "N"),
    m48 = ifelse(any(Flwtime.m > 48), "Y", "N"),
    m60 = ifelse(any(Flwtime.m > 60), "Y", "N"),
    m72 = ifelse(any(Flwtime.m > 72), "Y", "N"),
    m84 = ifelse(any(Flwtime.m == 84), "Y", "N")
  ) %>%
  ungroup()
head(hlme.class.1)

hlme.class.1.uni <- hlme.class.1 %>% group_by(ID) %>% distinct(ID, .keep_all = TRUE)
write.xlsx(hlme.class.1.uni, "3.3.hlme.class.uni.xlsx")

#统计
library(dplyr)
library(tibble)

#定义统计函数
count_Y <- function(df, vlgroup, class_val) {
  filtered_df <- df %>% filter(VLgroup == vlgroup, class == class_val)
  results <- colSums(filtered_df[, c("m12", "m24", "m36", "m48", "m60", "m72", "m84")] == "Y")
  return(results)
}

results_df <- tibble(
  VLgroup = character(),
  Class = integer(),
  m12 = integer(),
  m24 = integer(),
  m36 = integer(),
  m48 = integer(),
  m60 = integer(),
  m72 = integer(),
  m84 = integer()
)

# 定义VLgroup和class的组合
combinations <- list(
  list(vlgroup = "VS", class_val = 1),
  list(vlgroup = "VS", class_val = 2),
  list(vlgroup = "Blips", class_val = 1),
  list(vlgroup = "Blips", class_val = 2),
  list(vlgroup = "LLV", class_val = 1),
  list(vlgroup = "LLV", class_val = 2)
)

#循环统计每个组合的结果并存储到数据框
for (comb in combinations) {
  stats <- count_Y(hlme.class.1.uni, comb$vlgroup, comb$class_val)
  results_df <- results_df %>%
    add_row(
      VLgroup = comb$vlgroup,
      Class = comb$class_val,
      m12 = stats["m12"],
      m24 = stats["m24"],
      m36 = stats["m36"],
      m48 = stats["m48"],
      m60 = stats["m60"],
      m72 = stats["m72"],
      m84 = stats["m84"]
    )
}

# 打印表格
print(results_df)
