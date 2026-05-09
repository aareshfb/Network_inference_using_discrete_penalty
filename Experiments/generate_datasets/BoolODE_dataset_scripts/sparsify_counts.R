library(vroom)
library(readxl)
library(tidyverse)
library(glue)

list.files()


k <- 10
for(jj in 1:k){
  cxid <- glue('cluster_{jj}')
  message(cxid)
  expr <- vroom(file.path(cxid, 'ExpressionData.csv')) 
  expr <- column_to_rownames(expr, '...1')
  expr <- as.matrix(expr)
  print(dim(expr))
  
  ## make mask for sparsity
  mask <- matrix(data = sample(c(0,1), size = prod(dim(expr)), prob = c(.8, .2), replace = T),
                 nrow = nrow(expr))
  expr <- mask * expr
  expr <- expr %>% as.data.frame() %>% rownames_to_column('gene')
  write_excel_csv(expr, file = glue('{cxid}_sparse_expression.csv'))
}



