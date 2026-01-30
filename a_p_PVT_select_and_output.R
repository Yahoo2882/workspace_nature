library(Seurat)
library(dplyr)

# sub_all即已经筛选过得e-T1和e-T2神经元集合
DefaultAssay(sub_all) <- "RNA"

# 1) aPVT/pPVT markers(根据文档)
a_markers <- c("Drd1","Parm1","Cbln2","Mylk","Pbx3","Gda","Ntk1")
p_markers <- c("Drd2","Slc30a3","Snx31","Col12a1","Col25a1")

# 2) 取交集
a_use <- intersect(a_markers, rownames(sub_all))
p_use <- intersect(p_markers, rownames(sub_all))

a_use; p_use 

# 3) 模块打分
sub_all <- AddModuleScore(sub_all, features = list(a_use), name = "aScore")
sub_all <- AddModuleScore(sub_all, features = list(p_use), name = "pScore")

# 设置合理的delta阈值
diff <- sub_all$aScore1 - sub_all$pScore1
check_delta <- function(delta){
  ambiguous_n <- sum(abs(diff) <= delta, na.rm = TRUE)
  total_n <- sum(!is.na(diff))
  data.frame(delta = delta,
             ambiguous_n = ambiguous_n,
             ambiguous_pct = round(ambiguous_n/total_n*100, 2))
}

do.call(rbind, lapply(c(0.02, 0.05, 0.08, 0.10, 0.15), check_delta))


# 4) 贴标签
delta <- 0.05
diff <- sub_all$aScore1 - sub_all$pScore1
sub_all$a_pvt_label <- ifelse(diff >  delta, "aPVT",
                              ifelse(diff < -delta, "pPVT", "ambiguous"))

# sub_all$a_pvt_label <- ifelse(sub_all$aScore1 >= sub_all$pScore1, "aPVT", "pPVT")
# 总体
table(sub_all$a_pvt_label)
# 分期
table(sub_all$age, sub_all$a_pvt_label)

# 5) 删除anbiguous细胞
sub_ap <- subset(sub_all, subset = a_pvt_label %in% c("aPVT", "pPVT"))
table(sub_ap$a_pvt_label)
table(sub_ap$age, sub_ap$a_pvt_label)
sub_ap$cell.type <- sub_ap$a_pvt_label
Idents(sub_ap) <- "cell.type"

# 6)输出
saveRDS(sub_ap, file = file.path(out_dir, "subset_aPVT_pPVT.rds"))


