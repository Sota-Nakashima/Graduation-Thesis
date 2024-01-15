SELECT id,species,brain,heart,kidney,liver,testes 
FROM mammal_expression 
WHERE gene_type = "protein_coding"
