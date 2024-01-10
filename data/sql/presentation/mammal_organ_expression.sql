select id,brain,heart,kidney,liver,testes from mammal_expression 
where gene_type = 'protein_coding' and species = 'hg19';
