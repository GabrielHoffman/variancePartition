
# test when some expression values are missing

test_missing = function(){

  library(variancePartition)
  library(RUnit)

  set.seed(123)
  expr <- matrix(rnorm(6 * 6), nrow = 6)
  rownames(expr) <- paste0("Gene", 1:6)
  colnames(expr) <- paste0("Sample", 1:6)
   
  metadata <- data.frame(
      sample = colnames(expr),
      group = factor(rep(c("A", "B", "C"), each = 2)),
      subject = rep(c("S1", "S2"), 3)  
  )
  rownames(metadata) <- metadata$sample
   
  group_B_samples <- metadata$sample[metadata$group == "B"]
  expr["Gene1", group_B_samples] <- NA
   
  group_C_samples <- metadata$sample[metadata$group == "C"]
  expr["Gene2", group_C_samples] <- NA
   
  form <- ~ group 
  fit1 <- variancePartition::dream(expr, form, metadata)
  fit1 <- eBayes(fit1)

  form <- ~ group + (1|subject)
  fit2 <- variancePartition::dream(expr, form, metadata)
  fit2 <- eBayes(fit2)

  tab1 = topTable(fit1, coef="groupC", sort.by="none", number=Inf)
  tab2 = topTable(fit2, coef="groupC", sort.by="none", number=Inf) 
  checkEqualsNumeric( tab1$logFC, tab2$logFC, tolerance=1e-6)


  tab1 = topTable(fit1, coef="groupB", sort.by="none", number=Inf)
  tab2 = topTable(fit2, coef="groupB", sort.by="none", number=Inf) 
  checkEqualsNumeric( tab1$logFC, tab2$logFC, tolerance=1e-6)

  library(variancePartition)
  library(RUnit)

  data(varPartData)

  # simualte Group so that random effect explains zero variance
  info$Group = sample(LETTERS[1:2], nrow(info), replace=TRUE)

  form <- ~ Batch 
  fit1 <- dream(geneExpr, form, info)
  fit1 <- eBayes(fit1)

  form <- ~ Batch + (1 | Group)
  fit2 <- dream(geneExpr, form, info)
  fit2 <- eBayes(fit2)

  # check that F-statistics from fixed and random model
  # are close enough
  a = topTable(fit1, sort.by="none")$F
  b = topTable(fit2, sort.by="none")$F
  checkEquals(abs(mean(a-b)) < .05, TRUE)

  # run random effect with Group
  form <- ~ Batch + (1 | Group)
  geneExpr[1,info$Batch == 1] = NA
  fit3 <- dream(geneExpr, form, info)
  fit3 <- eBayes(fit3)

  head(coef(fit3))

  # on the second coef of the first gene should be NA
  checkEquals(which(is.na(coef(fit3))), 201)

  # check that F-statistics from fixed and random model
  # are close enough
  a = topTable(fit2, coef=3:4, sort.by="none")$F
  b = topTable(fit3, coef=3:4, sort.by="none")$F
  checkEquals(abs(mean(a-b)) < .05, TRUE)

}










