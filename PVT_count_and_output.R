library(Seurat)
library(dplyr)

stopifnot(all(c("age", "cell.type") %in% colnames(obj@meta.data)))
age_levels <- c("e16","e18","p0","p4","p10","p18","p28","p65")
obj$age <- factor(obj$age, levels = age_levels)


target_types <- c("e-T1", "e-T2")

# 1) 筛出所有 e-T1/e-T2 并计数并输出到一个数据集
sub_all <- subset(obj, subset = `cell.type` %in% target_types)

# 计数
count_by_type <- as.data.frame(table(sub_all$`cell.type`))
colnames(count_by_type) <- c("cell.type", "n_cells")
print(count_by_type)

# 输出
out_dir  <- "/home/yahoo/Desktop/workspace_nature/output"
saveRDS(sub_all, file = file.path(out_dir, "subset_eT1_eT2_all.rds"))


# 2) 按age和cell.type
count_by_age_type <- as.data.frame(table(sub_all$age, sub_all$`cell.type`))
colnames(count_by_age_type) <- c("age", "cell.type", "n_cells")
count_by_age_type <- count_by_age_type %>%
  arrange(factor(age, levels = age_levels), cell.type)
print(count_by_age_type)

# 为每个时期输出一个Seurat子对象
sub_by_age <- setNames(vector("list", length(age_levels)), age_levels)

for (ag in age_levels) {
  sub_ag <- subset(sub_all, subset = age == ag)
  
  # 计数
  tab_ag <- as.data.frame(table(sub_ag$`cell.type`))
  colnames(tab_ag) <- c("cell.type", "n_cells")
  tab_ag$age <- ag
  tab_ag <- tab_ag[, c("age", "cell.type", "n_cells")]
  print(tab_ag)

  # 输出
  saveRDS(sub_ag, file = file.path(out_dir, paste0("subset_eT1_eT2_", ag, ".rds")))
  
}

