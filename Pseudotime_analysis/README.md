# Pseudotime Analysis with Palantir

Pseudotime was modeled using the Palantir package (https://github.com/dpeerlab/Palantir). Palantir was used onto the dorsal aorta (Arterial ECs, Pre-HE and HE celltypes) and IAHCs celltypes. 

![UMAP](../images/UMAP.png)

Pseudotime trajectory was executed with aorta ECs as the starting state (highest expression of the Bmx marker) and IAHCs as the terminal state (highest expression of Myc)

![Selected_cells](../images/Selected_Cells_Palantir.png)

As Endothelial to Hematopoietic transition is a one way event with no subsequent divergences, the umber of branching paths was set to 1.

![Branching](../images/Branching_paths_Palantir.png)

A pseudotime gradient was produced that would group the cells across the pseudotime trayectory.

![Pseudotime](../images/Palantir_Pseudotime.png)
