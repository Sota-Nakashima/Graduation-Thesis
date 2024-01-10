SELECT fed.`Gene ID`, fed.Run, fed.TPM 
FROM Fruitfly_Expression_Data_ver2 AS fed 
INNER JOIN ( 
    SELECT DISTINCT gene_ID, transcript_type 
    FROM Flybase_ID_Convert 
    WHERE transcript_type = 'ncRNA' 
) AS fbc 
ON fed.`Gene ID` = fbc.gene_ID 
INNER JOIN Fruitfly_Expression_Index AS fei 
ON fed.Run = fei.Run;
