library(data.table)

ru_viral_expr_chi <- fread("output/deseq2/ruakura/ru_viral_expression_chi/ru_viral_expression_chi.csv")
chi_sq_res <- list(chisq.test(table(ru_viral_expr_chi)))
##chi-sq p value = 0.676, viral expression independent of behaviour
fwrite(chi_sq_res, "output/deseq2/ruakura/ru_viral_expression_chi/chi_sq_res.txt")