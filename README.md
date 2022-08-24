# TrabajoFinDeMaster_MC2
Uso de datos masivos de hepatocarcinoma para estratificar pacientes y predecir nuevas estrategias terapéuticas

## Objetivo:
Estudios previos basados en el análisis transcriptómicos de genes codificantes han conseguido estratificar pacientes con Hepatocarcinoma en agrupaciones con relevancia clínica. Teniendo en cuenta que los genes no codificantes son más específicos de tumor y tejido que los codificantes, pensamos que podrían permitir una clasificación más relevante clínicamente de los pacientes, en especial, con una asociación más significativa con la supervivencia. 
Para analizar si el transcriptoma no codificante permite la estratificación de los pacientes, en este proyecto se utilizarán datos de secuenciación de RNA de distintas cohortes. Se compararán diversas técnicas y se realizará una evaluación crítica de los resultados obtenidos con cada una de ellas. 

---

En esta web se encuentran los ficheros y scripts utilizados a lo largo del proyecto:

1. Muestras de pacientes:
  - 01_NASIR_binding.r
  - 01_NASIR_P.csv
  - 01_NASIR_T.csv


2,4. Normalización, Differential Gene Expression Analysis - Gene Set Enrichment Analysis (GSEA):
  - 02_NASIR_Differential_gene_expression_analysis.r
    - 02_resSig.rds
    - 03_NASIR_FINAL_TABLE_COUNTS.7z
  - 07_NASIR_Gene_Annotation_Gene_Ontology.r
  - 11_LICA_Normalization.r
   - 10_LICA_Normalized_Dataset.7z
  - 12_LICA_Gene_Annotation.r
  - 14_LICA_Differential_Expression_Gene_Set_Enrichment_Analysis.r
    - 14_1_SummarizedExperiment_quantile_58k.rds
    - 14_2_resSig_anotado_cuantiles58k.rds
    - 14_3_DESeq2_Sig_results_anotado_cuantiles58k.tab


3. Análisis de supervivencia:
  - 03_NASIR_survival.r
  - 10_LICA_Cox_Regression.py
  -  10_LICA_Cox_Regression.ipynb
    - 09_4_LICA_survival.xlsx

5. Modelos predictivos:

    5.1. Clustering Jerárquico:
    - 09_1_LICA_Agglomerative_Hierarchical_Clustering.py
      - 09_1_LICA_Agglomerative_Hierarchical_Clustering.ipynb
    - 09_2_LICA_Agglomerative_Hierarchical_Clustering.py
       - 09_2_LICA_Agglomerative_Hierarchical_Clustering.ipynb
    - 09_3_LICA_Agglomerative_Hierarchical_Clustering.py
       - 09_3_LICA_Agglomerative_Hierarchical_Clustering.ipynb
    - 09_4_LICA_Agglomerative_Hierarchical_Clustering.py
       - 09_4_LICA_Agglomerative_Hierarchical_Clustering.ipynb

    5.2. Deep learning:

      5.2.1.NASIR:
      - 04_NASIR_H2O_Models.r
      - 05_NASIR_R100_HCCDeReg.r
      - 06_NASIR_H2O_Models_R100_HCCDeReg.r
      - 08_NASIR_H2O_Curated_Data.r

      5.2.2. LICA:
      - 13_1_LICA_H2O_alcohol.r
      - 13_2_LICA_H2O_edmonson.r
      - 13_3_LICA_H2O_fibrosis.r
      - 13_4_LICA_H2O_vascular_invasion.r
      - 13_5_LICA_H2O_survival_iprot_coding.r
      - 13_6_1_LICA_H2O_survival_decil.r
      - 13_6_2_LICA_H2O_survival_tercil.r
      - 13_6_LICA_H2O_survival_pseudo_nc.r
      - 13_7_LICA_H2O_survival_binary.r
