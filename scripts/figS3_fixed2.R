setwd("path/to/ARDS_github")
suppressPackageStartupMessages({ library(pROC); library(ggplot2); library(dplyr) })

clean <- function(x) gsub("\\.", "-", x)

# ===== 1. 训练固定模型 =====
scores <- as.matrix(read.csv("data_processed/ssGSEA_32707.csv", row.names=1, check.names=FALSE))
grp <- read.csv("data_processed/GSE32707_group.csv")
lbl <- ifelse(grp$group=="Disease",1,0); names(lbl) <- grp$gsm
cm <- intersect(colnames(scores), names(lbl)); scores <- scores[,cm]; lbl <- lbl[cm]
key <- c("Immature_neutrophils","Non-classical_monocytes")
df_train <- data.frame(t(scores[key,]), label=lbl)
colnames(df_train) <- c(key, "label")
fit_fixed <- glm(label~., data=df_train, family=binomial)
cat("系数:\n"); print(round(coef(fit_fixed),4))

# ===== 2. GSE243066 =====
s243 <- as.matrix(read.csv("data_processed/ssGSEA_243066.csv", row.names=1, check.names=FALSE))
l243 <- ifelse(grepl("^A", colnames(s243)), 1, 0)
df243 <- data.frame(t(s243[key,]))
colnames(df243) <- key
pred243 <- predict(fit_fixed, newdata=df243, type="response")
roc243 <- roc(l243, pred243, quiet=TRUE)
cat("GSE243066 fixed AUC:", round(auc(roc243),3), "\n")

# ===== 3. GSE200847 =====
s200 <- as.matrix(read.csv("data_processed/ssGSEA_GSE200847.csv", row.names=1, check.names=FALSE))
gm <- c(TA203="Hypo",TA213="Hypo",TA225="Hypo",TA227="Hyper",TA229="Hyper",TA234="Hypo",TA235="Hypo",
  TA251="Hypo",TA257="Hypo",TA288="Hypo",TA311="Hyper",TA314="Hyper",TA315="Hypo",TA319="Control",
  TA320="Hyper",TA334="Hypo",TA335="Hyper",TA337="Hyper",TA341="Hypo",TA353="Hyper",TA358="Hypo",
  TA392="Hypo",TA400="Hypo",TA402="Hypo",TA409="Hypo",TA411="Hypo",TA413="Hypo",TA414="Hypo",
  TA458="Hypo",TA462="Control",TA463="Hypo",TA477="Control",TA479="Hypo",TA488="Hypo",TA503="Hyper",
  TA524="Hypo",TA530="Hypo",TA539="Hypo",TA540="Hypo",TA541="Hypo",TA543="Hypo",TA551="Hypo",
  TA555="Control",TA563="Control",TA565="Hyper",TA575="Hypo")
l200 <- ifelse(gm[colnames(s200)] %in% c("Hyper","Hypo"), 1, 0)
key200 <- intersect(key, rownames(s200))
df200 <- data.frame(t(s200[key200,]))
colnames(df200) <- key200
pred200 <- predict(fit_fixed, newdata=df200, type="response")
roc200 <- roc(l200, pred200, quiet=TRUE)
cat("GSE200847 fixed AUC:", round(auc(roc200),3), "\n")

# ===== 4. GSE171524 =====
s171 <- as.matrix(read.csv("data_processed/pb171524_scores.csv", row.names=1, check.names=FALSE))
l171 <- ifelse(grepl("^L", colnames(s171)), 1, 0)
key171 <- intersect(key, rownames(s171))
df171 <- data.frame(t(s171[key171,]))
colnames(df171) <- key171
pred171 <- predict(fit_fixed, newdata=df171, type="response")
roc171 <- roc(l171, pred171, quiet=TRUE)
cat("GSE171524 fixed AUC:", round(auc(roc171),3), "\n")

# ===== 5. FigS3A =====
png("data_processed/FigS3A_fixed.png", width=700, height=700, res=150)
plot(roc243, col="red", lwd=2)
plot(roc200, col="blue", lwd=2, add=TRUE)
plot(roc171, col="darkgreen", lwd=2, add=TRUE)
abline(a=1, b=1, lty=2, col="gray")
legend("bottomright",
  legend=c(paste0("GSE243066 (",round(auc(roc243),3),")"),
           paste0("GSE200847 (",round(auc(roc200),3),")"),
           paste0("GSE171524 (",round(auc(roc171),3),")")),
  col=c("red","blue","darkgreen"), lwd=2, bty="n", cex=0.9)
dev.off()

# ===== 6. FigS3B =====
fd <- data.frame(
  cohort=c("GSE243066","GSE200847","GSE171524"),
  AUC=c(as.numeric(auc(roc243)), as.numeric(auc(roc200)), as.numeric(auc(roc171))),
  lo=c(ci.auc(roc243)[1], ci.auc(roc200)[1], ci.auc(roc171)[1]),
  hi=c(ci.auc(roc243)[3], ci.auc(roc200)[3], ci.auc(roc171)[3]))
write.csv(fd, "data_processed/FigS3B_fixed_data.csv", row.names=FALSE)
p <- ggplot(fd, aes(x=AUC, y=reorder(cohort,AUC))) +
  geom_point(size=4, color="red") +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.2, color="red") +
  geom_vline(xintercept=0.5, linetype="dashed", color="gray") + xlim(0.4,1) +
  theme_minimal() + labs(x="AUC (95% CI)", y="")
ggsave("data_processed/FigS3B_fixed.png", p, width=7, height=4, dpi=150)

# ===== 7. FigS3C =====
df_box <- rbind(
  data.frame(cohort="GSE243066", group=ifelse(l243==1,"ARDS","Control"), score=pred243),
  data.frame(cohort="GSE200847", group=ifelse(l200==1,"ARDS","Control"), score=pred200),
  data.frame(cohort="GSE171524", group=ifelse(l171==1,"ARDS","Control"), score=pred171))
p <- ggplot(df_box, aes(x=cohort, y=score, fill=group)) +
  geom_boxplot(outlier.size=0.5) +
  scale_fill_manual(values=c("ARDS"="salmon","Control"="lightblue")) +
  theme_minimal() + labs(x="", y="Predicted probability", fill="")
ggsave("data_processed/FigS3C_fixed.png", p, width=8, height=5, dpi=150)

cat("=== Figure S3 完成 ===\n")
