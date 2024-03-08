# Manhattan-Plots 

A Manhattan plot is a graphical representation used in Genome-Wide Association Studies (GWAS) to visualize and identify potential genetic variants associated with a particular trait or disease. It's named "Manhattan" due to its skyline-like appearance, where each data point resembles a building on a cityscape.

Here's how a Manhattan plot is typically constructed:

1. **X-axis: Chromosomes**
   - The X-axis represents different chromosomes in the genome. Each chromosome is divided into segments, and these segments are often labeled with chromosome numbers.

2. **Y-axis: -log10(p-value)**
   - The Y-axis represents the statistical significance of the association between genetic variants and the trait or disease. The p-values from the statistical tests are often transformed by taking the negative logarithm (base 10) to make them more manageable and to emphasize highly significant results.

3. **Data Points: Genetic Variants**
   - Each point on the plot corresponds to a genetic variant (e.g., a single nucleotide polymorphism or SNP) tested in the GWAS. The position on the X-axis indicates the location of the variant on the chromosome, and the height on the Y-axis reflects the strength of the association (smaller p-values are higher).

4. **Threshold Line: Significance Threshold**
   - A horizontal line is drawn on the plot to represent the significance threshold. This threshold is often set based on a predefined significance level (e.g., 5 x 10^(-8)). Variants above this line are considered potentially significant associations.

Interpreting a Manhattan plot involves looking for "peaks" or "towers" that rise above the significance threshold. These peaks suggest regions in the genome where genetic variants may be associated with the trait or disease being studied.

In summary, the Manhattan plot is a valuable tool in GWAS for visually identifying genomic regions that may harbor genetic variants associated with a particular phenotype. Researchers use this visualization to prioritize regions for further investigation and exploration of potential causal variants.
