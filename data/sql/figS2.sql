select gene_type,testes from mammal_expression 
where gene_type = 'protein_coding' OR gene_type = 'lincRNA';
