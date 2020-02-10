# Antibody cross-reaction dynamics
Code for project testing/simulation platform 


These files form a testing platform for simulating the antibodies (Abs) cross-reaction dynamics to an antigen, in this case Influenza.
'funct_shift' is for code configuration and functions.

1) Run 'H3_platform_1' to simulate dynamics of naive B-cells and plasma cells. The influenza strains (H3N2) are automatically loaded from 'funct_shift'. 
Output: 'Bcel_fina_1.csv'- is B-cell counts,
        'Bcell_1.csv'- is B-cell coordinates on binary-map,
        'Plas_final_1.csv'- is plasma cell counts,
        'PDFs'- Figures.

2) Run 'Abs_H3_platform_1' to simulate the Abs dynamics. 
Output: 'Abs_total_platform_1.csv'- is Abs counts,
        'PDFs'- Figure.
        
Further details: https://www.biorxiv.org/content/10.1101/2020.01.06.896274v1

Cite the corresponding tool as follows: Hernandez-Mejia G. and Hernandez-Vargas E.A., Uncovering antibody cross-reaction dynamics in influenza A infections, 2020.
