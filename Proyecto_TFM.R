# 1. Cargar librerías y establecemos el directorio de trabajo
install.packages("bannerCommenter")
library(bannerCommenter)
library(dplyr)
library(limma)
library(tibble)
library(ggpubr)
setwd("~/Desktop/TCGA LAML")


# 2. Cargar matriz de expresión
data <- read.table("TCGA.LAML.sampleMap_HiSeqV2", header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# 3. Cargar datos clínicos
data_clinical <- read.table("c07a64a0-7588-4653-95ef-982b41a1a804/nationwidechildrens.org_clinical_patient_laml.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 4. Arreglar identificadores de pacientes para que coincidan
data_clinical$bcr_patient_barcode <- gsub("-", ".", data_clinical$bcr_patient_barcode)
colnames(data) <- gsub("\\.03$", "", colnames(data))

# 5. Definir los grupos en la base clínica
data_clinical <- data_clinical %>%
  mutate(cytogenetic_group = case_when(
    cytogenetic_abnormality_type %in% c("Normal", "Normal|del(7q) / 7q-", "Normal|t(8;21)", "Normal|+8", "Normal|inv(16)",
                                        "Normal|t(15;17) and variants", "Normal|t(9;11)", "Normal|del(7q) / 7q-|t(9;11)",
                                        "Normal|del(5q) / 5q-|del(7q) / 7q-") ~ "Normal",
    cytogenetic_abnormality_type %in% c("t(8;21)", "inv(16)", "del(5q) / 5q-", "del(7q) / 7q-", "t(15;17) and variants",
                                        "+8", "+8|t(15;17) and variants") ~ "Simple",
    cytogenetic_abnormality_type %in% c("Complex - >= 3 distinct abnormalities", "Complex - >= 3 distinct abnormalities|+8",
                                        "Complex - >= 3 distinct abnormalities|del(5q) / 5q-", 
                                        "Complex - >= 3 distinct abnormalities|t(15;17) and variants",
                                        "Complex - >= 3 distinct abnormalities|del(5q) / 5q-|+8",
                                        "Complex - >= 3 distinct abnormalities|del(5q) / 5q-|del(7q) / 7q-|+8",
                                        "Complex - >= 3 distinct abnormalities|del(5q) / 5q-|del(7q) / 7q-",
                                        "Complex - >= 3 distinct abnormalities|del(7q) / 7q-|+8",
                                        "Normal|Complex - >= 3 distinct abnormalities",
                                        "Normal|Complex - >= 3 distinct abnormalities|del(7q) / 7q-",
                                        "Normal|Complex - >= 3 distinct abnormalities|del(7q) / 7q-|+8",
                                        "Normal|Complex - >= 3 distinct abnormalities|+8") ~ "Complex",
    TRUE ~ "Other" 
  ))

# 6. Filtrar solo los pacientes que están en ambas bases
data_clinical_filtered <- data_clinical %>% filter(bcr_patient_barcode %in% colnames(data))

# 7. Ordenar la base clínica por grupo citogenético
data_clinical_ordered <- data_clinical_filtered %>%
  arrange(factor(cytogenetic_group, levels = c("Normal", "Simple", "Complex", "Other")))

# 8. Reordenar la matriz de expresión según la base clínica ordenada
data_ordered <- data[, data_clinical_ordered$bcr_patient_barcode]

# 9. Convertir en matriz y asegurarse de que los nombres coincidan
data_matrix <- as.matrix(data_ordered)
rownames(data_matrix) <- rownames(data)

# 10. Definir los grupos en base al orden ya establecido
grupo <- factor(data_clinical_ordered$cytogenetic_group, levels = c("Normal", "Simple", "Complex", "Other"))

# 11. Crear matriz de diseño para limma
experimental.design <- model.matrix(~ -1 + grupo)
colnames(experimental.design) <- levels(grupo)
rownames(experimental.design) <- data_clinical_ordered$bcr_patient_barcode  # Asegurar que las filas coincidan

# 12. Asegurar que las dimensiones coincidan
print(dim(data_matrix))         
print(dim(experimental.design))

# 13. Aplicar modelo lineal con limma
linear.fit <- lmFit(data_matrix, experimental.design)

# 14. Definir matriz de contrastes
contrast.matrix <- makeContrasts(
  Complejos_vs_Normales = Complex - Normal, 
  Simples_vs_Normales = Simple - Normal,
  levels = colnames(experimental.design)
)

# 15. Ajustar modelos con contrastes y aplicar eBayes
contrast.linear.fit <- contrasts.fit(linear.fit, contrast.matrix)
contrast.results <- eBayes(contrast.linear.fit)

# 16. Obtener la tabla de resultados con TopTable
Complex.Normal <- topTable(contrast.results, number=20530, coef="Complejos_vs_Normales", sort.by="logFC")
Simple.Normal <- topTable(contrast.results, number=20530, coef="Simples_vs_Normales", sort.by="logFC")

# 17. Definimos los umbrales
padj.cutoff = 0.05   # p-valor ajustado (significancia estadística)
lfc.cutoff = 0.58    # log2 fold-change (equivalente a FC > 1.5 o < 0.67)


# 18. Convertimos la tabla de resultados en Tibble
Complex.Normal_Tibble = Complex.Normal %>%
  data.frame() %>%
  rownames_to_column(var = "Gene") %>%  # Convertir rownames a columna "Gene"
  as_tibble()

Simple.Normal_Tibble = Simple.Normal %>%
  data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  as_tibble()

# 19. Filtramos los DEGs significativos
resLFC_Sig_Complex = Complex.Normal_Tibble %>%
  dplyr::filter(P.Value < padj.cutoff & abs(logFC) > lfc.cutoff)

resLFC_Sig_Simple = Simple.Normal_Tibble %>%
  dplyr::filter(P.Value < padj.cutoff & abs(logFC) > lfc.cutoff)

# 20. Visualizar los primeros genes obtenidos
cat("\nTotal DEGs Complex vs Normal:", nrow(resLFC_Sig_Complex), "\n\n")
cat("\nTotal DEGs Simple vs Normal:", nrow(resLFC_Sig_Simple), "\n")


# 21. Creación de tablas

# Instalar y cargar paquete necesario
install.packages("officer")
library(officer)

# Crear un documento Word
doc <- read_docx()

# Agregar la primera tabla (Normal vs Complex)
doc <- doc %>%
  body_add_par("Top DEGs: Normal vs Complex", style = "heading 1") %>%
  body_add_table(resLFC_Sig_Complex, style = "table_template")

# Agregar la segunda tabla (Normal vs Simple)
doc <- doc %>%
  body_add_par("Top DEGs: Normal vs Simple", style = "heading 1") %>%
  body_add_table(resLFC_Sig_Simple, style = "table_template")

# Guardar el documento Word
print(doc, target = "DEGs_Results_Comparison.docx")


# 22. Contamos los genes sobreexpresados y subexpresados para cada comparación

genes_sobreexp_Complex <- sum(resLFC_Sig_Complex$logFC > 0)
genes_subexp_Complex <- sum(resLFC_Sig_Complex$logFC < 0)
total_DEGs_Complex <- nrow(resLFC_Sig_Complex)

genes_sobreexp_Simple <- sum(resLFC_Sig_Simple$logFC > 0)
genes_subexp_Simple <- sum(resLFC_Sig_Simple$logFC < 0)
total_DEGs_Simple <- nrow(resLFC_Sig_Simple)

# Mostrar resultados en formato de tabla
cat("\nComparación\t\tGenes sobreexpresados\tGenes subexpresados\tTotal DEGs\n")
cat("Complex vs Normal\t", genes_sobreexp_Complex, "\t\t\t", genes_subexp_Complex, "\t\t\t", total_DEGs_Complex, "\n")
cat("Simple vs Normal\t", genes_sobreexp_Simple, "\t\t\t", genes_subexp_Simple, "\t\t\t", total_DEGs_Simple, "\n")


# 23. Realizamos una representación gráfica de nuestros resultados

# 1. Cargar librerías necesarias
install.packages("ggplot2")
install.packages("ggrepel")
library(ggplot2)
library(ggrepel)
library(dplyr)
library(limma)
library(tibble)
library(ggpubr)

# 22. Contamos los genes sobreexpresados y subexpresados para cada comparación
genes_sobreexp_Complex <- sum(resLFC_Sig_Complex$logFC > 0)
genes_subexp_Complex <- sum(resLFC_Sig_Complex$logFC < 0)
total_DEGs_Complex <- nrow(resLFC_Sig_Complex)

genes_sobreexp_Simple <- sum(resLFC_Sig_Simple$logFC > 0)
genes_subexp_Simple <- sum(resLFC_Sig_Simple$logFC < 0)
total_DEGs_Simple <- nrow(resLFC_Sig_Simple)

# 23. Crear Pie Chart para representar genes sobreexpresados y subexpresados

datos_pie_Complex <- data.frame(
  Categoria = c("Sobreexpresados", "Subexpresados"),
  Cantidad = c(genes_sobreexp_Complex, genes_subexp_Complex)
)

datos_pie_Simple <- data.frame(
  Categoria = c("Sobreexpresados", "Subexpresados"),
  Cantidad = c(genes_sobreexp_Simple, genes_subexp_Simple)
)

# Pie chart para Complex vs Normal
ggplot(datos_pie_Complex, aes(x = "", y = Cantidad, fill = Categoria)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggtitle("Genes Diferencialmente Expresados: Complex vs Normal") +
  theme_void()

# Pie chart para Simple vs Normal
ggplot(datos_pie_Simple, aes(x = "", y = Cantidad, fill = Categoria)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  ggtitle("Genes Diferencialmente Expresados: Simple vs Normal") +
  theme_void()


# 24. Crear Volcano Plots para cada comparación


# Definimos la función generar_volcano_plot con etiquetas
generar_volcano_plot <- function(data_input, title) {
  
  # Verificamos que data_input es un data.frame
  if (!inherits(data_input, "data.frame")) {
    stop("El objeto data_input no es un data frame válido")
  }
  
  # Creamos una nueva columna 'significativo' para clasificar los genes
  data_input <- data_input %>%
    mutate(significativo = case_when(
      P.Value < padj.cutoff & logFC > lfc.cutoff ~ "Sobreexpresado",
      P.Value < padj.cutoff & logFC < -lfc.cutoff ~ "Subexpresado",
      TRUE ~ "No Significativo"
    ))
  
  # Seleccionamos los genes más significativos para etiquetar
  significant_genes <- data_input %>%
    filter(significativo != "No Significativo") %>%
    top_n(10, wt = -log10(P.Value))  # Seleccionar los 10 genes más significativos
  
  # Creamos el volcano plot utilizando ggplot2
  plot <- ggplot(data_input, aes(x = logFC, y = -log10(P.Value), color = significativo)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("No Significativo" = "gray70", 
                                  "Sobreexpresado" = "#FF5733", 
                                  "Subexpresado" = "#3498DB")) +
    geom_hline(yintercept = -log10(padj.cutoff), linetype = "dashed", color = "darkblue") +
    geom_vline(xintercept = c(-lfc.cutoff, lfc.cutoff), linetype = "dashed", color = "darkblue") +
    ggtitle(title) +
    theme_minimal() +
    geom_text_repel(data = significant_genes, aes(label = rownames(significant_genes)), size = 3, box.padding = 0.3)  # Etiquetas usando rownames
  
  # Mostrar el gráfico
  print(plot)
}

# Establecer los umbrales globales para el p-valor ajustado y el log2 fold change
padj.cutoff = 0.05   # p-valor ajustado (significancia estadística)
lfc.cutoff = 0.58    # log2 fold-change (equivalente a FC > 1.5 o < 0.67)


# Generamos los Volcano Plots para los datos de cada comparación
generar_volcano_plot(Complex.Normal, "Volcano Plot (Todos los genes): Complex vs Normal")
generar_volcano_plot(Simple.Normal, "Volcano Plot (Todos los genes): Simple vs Normal")



# 23. Realizamos un enriquecimiento

# Instalar y cargar clusterProfiler y otras librerías necesarias
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
install.packages(c("clusterProfiler", "org.Hs.eg.db", "ggplot2"))
library(clusterProfiler)
library(org.Hs.eg.db)  # Base de datos de genes humanos
library(ggplot2)
library(enrichplot)

# Convertir nombres de genes a IDs de ENTREZ
genes_complex <- bitr(resLFC_Sig_Complex$Gene, 
                      fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
  na.omit()  # Eliminar genes sin ID

genes_simple <- bitr(resLFC_Sig_Simple$Gene, 
                     fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
  na.omit()  

# Enriquecimiento GO (Biological Process)
go_complex <- enrichGO(gene = genes_complex$ENTREZID, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05)

go_simple <- enrichGO(gene = genes_simple$ENTREZID, 
                      OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", 
                      ont = "BP", 
                      pAdjustMethod = "BH", 
                      pvalueCutoff = 0.05)


# Análisis de enriquecimiento KEGG
kegg_complex <- enrichKEGG(gene = genes_complex$ENTREZID, 
                           organism = "hsa", 
                           pvalueCutoff = 0.05)

kegg_simple <- enrichKEGG(gene = genes_simple$ENTREZID, 
                          organism = "hsa", 
                          pvalueCutoff = 0.05)


# Filtrar términos enriquecidos antes de graficar
go_complex_filtered <- go_complex
go_complex_filtered@result <- go_complex@result %>%
  filter(p.adjust < 0.01) %>%
  head(100)  # Limitar a los 100 términos más significativos

go_simple_filtered <- go_simple
go_simple_filtered@result <- go_simple@result %>%
  filter(p.adjust < 0.01) %>%
  head(100)

kegg_complex_filtered <- kegg_complex
kegg_complex_filtered@result <- kegg_complex@result %>%
  filter(p.adjust < 0.05) %>%
  head(50)

kegg_simple_filtered <- kegg_simple
kegg_simple_filtered@result <- kegg_simple@result %>%
  filter(p.adjust < 0.05) %>%
  head(50)

# Simplificar términos redundantes en GO
go_complex_simplified <- simplify(go_complex_filtered, cutoff = 0.7, by = "p.adjust", select_fun = min)
go_simple_simplified <- simplify(go_simple_filtered, cutoff = 0.7, by = "p.adjust", select_fun = min)


# Visualización del enriquecimiento

# a) Dotplot para términos GO enriquecidos
dotplot(go_complex, showCategory = 15, title = "GO BP Enrichment - Complex") + theme_minimal()
dotplot(go_simple, showCategory = 15, title = "GO BP Enrichment - Simple") + theme_minimal()

# b) Barplot para KEGG
barplot(kegg_complex, showCategory = 10, title = "KEGG Pathways - Complex")
barplot(kegg_simple, showCategory = 10, title = "KEGG Pathways - Simple")


# 6. Mostrar términos enriquecidos con los genes asociados
head(go_complex@result[, c("Description", "geneID")])
head(go_simple@result[, c("Description", "geneID")])
head(kegg_complex@result[, c("Description", "geneID")])
head(kegg_simple@result[, c("Description", "geneID")])

nrow(go_complex@result)
nrow(go_simple@result)
nrow(kegg_complex@result)
nrow(kegg_simple@result)



