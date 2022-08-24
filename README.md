# HCC_lncRNA
Uso de datos masivos de hepatocarcinoma para estratificar pacientes y predecir nuevas estrategias terapéuticas

## Objetivo:
Estudios previos basados en el análisis transcriptómicos de genes codificantes han conseguido estratificar pacientes con Hepatocarcinoma en agrupaciones con relevancia clínica. Teniendo en cuenta que los genes no codificantes son más específicos de tumor y tejido que los codificantes, pensamos que podrían permitir una clasificación más relevante clínicamente de los pacientes, en especial, con una asociación más significativa con la supervivencia. 

Para analizar si el transcriptoma no codificante permite la estratificación de los pacientes, en este proyecto se utilizarán datos de secuenciación de RNA de distintas cohortes. Se compararán diversas técnicas y se realizará una evaluación crítica de los resultados obtenidos con cada una de ellas. 

---

En esta web se encuentran los ficheros y scripts utilizados a lo largo del proyecto:

1. Muestras de pacientes:
 - 01_NASIR_binding.r


2,4. Normalización, Differential Gene Expression Analysis - Gene Set Enrichment Analysis (GSEA):
  - 02_NASIR_Differential_gene_expression_analysis.r
  - 04_NASIR_Gene_Annotation_Gene_Ontology.r
  - 4_LICA_Gene_Annotation.r
  - 4_LICA_Differential_Expression_Gene_Set_Enrichment_Analysis.r

3. Análisis de supervivencia:
  - 03_NASIR_survival.r
  - 03_LICA_Cox_Regression.py
  -  03_LICA_Cox_Regression.ipynb

5. Modelos predictivos:

    5.1. Clustering Jerárquico:
    - 05_1_1_LICA_Agglomerative_Hierarchical_Clustering.py
      - 05_1_1_LICA_Agglomerative_Hierarchical_Clustering.ipynb
    - 05_1_2_LICA_Agglomerative_Hierarchical_Clustering.py
       - 05_1_2_LICA_Agglomerative_Hierarchical_Clustering.ipynb
    - 05_1_3_LICA_Agglomerative_Hierarchical_Clustering.py
       - 05_1_3_LICA_Agglomerative_Hierarchical_Clustering.ipynb
    - 09_1_4_LICA_Agglomerative_Hierarchical_Clustering_Statistical_Study.r


    5.2. Deep learning:

      NASIR:
      - 05_2_NASIR_H2O_Models.r
      - 05_2_NASIR_R100_HCCDeReg.r
      - 05_2_NASIR_H2O_Models_R100_HCCDeReg.r
      - 05_2_NASIR_H2O_Curated_Data.r

      LICA:
      - 05_2_LICA_H2O_alcohol.r
      - 05_2_LICA_H2O_edmonson.r
      - 05_2_LICA_H2O_fibrosis.r
      - 05_2_LICA_H2O_vascular_invasion.r
      - 05_2_LICA_H2O_survival_iprot_coding.r
      - 05_2_1_LICA_H2O_survival_decil.r
      - 05_2_2_LICA_H2O_survival_tercil.r
      - 05_2_3_LICA_H2O_survival_pseudo_nc.r
      - 05_2_4_LICA_H2O_survival_binary.r

__*Se ha seguido la misma numeración que en la memoria.__

---

Versiones de las librerías de R y Python utilizadas durante el proyecto:
- R: r_library_version
- Python: python_library_version 
