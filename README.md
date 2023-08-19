# virtual_ecoli_mutants_libr

Given the development of a model that approximates well the growth phenotype of a typical E. coli cell, this provides now a straightforward approach to producing in silico E.coli mutants that “encode” different strategies for resource allocation and demonstrate different growth/bioproduction phenotypes. This process is straightforward and can be achieved by changing the chemical reaction rates of mRNA expression of the 3 coarse grain biomolecular classes in the model, and observing the effects of those changes in terms of growth rate and bioproduction rate.

To this end, we produced 6 “mutant” models of the original virtual cell model and simulated those models to observe any changes in the phenotypes of growth rate or bioproduction. The introduced “mutations” in the alternative model are given in the table below. 

Model name	Modification(s) in the transcription rate of biomolecular classes
model_mut1	hsk overexpression (x2)
model_mut2	met overexpression (x1.5)
model_mut3	rib overexpression (x1.5)
model_mut4	met overexpression (x1.5)
rib overexpression (x1.5)
model_mut5	hsk under-expression (x0.7)
model_mut6	met overexpression (x1.25)
rib overexpression (x1.25)
hsk under-expression (x0.7)
