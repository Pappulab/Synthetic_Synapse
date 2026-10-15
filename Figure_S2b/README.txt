This notebook extracts proteins associated with cell junctions, as well as proteins with IDRs and examines their domain and IDR properties

The following was used to extract the proteins of interest:
1. Accessions were extracted from Uniprot using the GO term 0030054 for Cell junctions (data/uniprotkb_organism_id_9606_AND_go_00300_2025_09_24.list)

2. These acessions were filtered to exclude any proteins that had a transmembrane domain using KW-0812 in Uniprot (data/uniprotkb_keyword_KW_0812_AND_model_org_2025_09_24.list)

3. IDR containing proteins were extracted using the method from the GIN manuscript. Briefly, MobiDB was used to predict all IDRs of length ≥ 30. Then a filtering was performed to extract IDRs of length ≥ 70 and non-linker IDRs of length ≥ 50. This criterion helps idenitfy IDRs likely to play functional roles with the context of the full protein. (data/All_sapiens_IDRs_clustermap_using_centroids.xlsx)

4. Additionally only accessions in which the length of the full protein was consistent in the IDR dataset and domain dataset were considered. The reason for some mismatches is that Uniprot changed which isoform was the main isoform between cases. Thus, to avoid mismatch in where domains should match to we just remove all cases of mismatch upfront

Using these criteria we get the following numbers:
Number of unique proteins in cell junctions without transmembrane domains: 1510
Number of unique proteins housing IDRs ≥ 70 and non-linker IDRs ≥ 50: 6269
Number of unique cell junction without transmembrane domain proteins housing IDRs ≥ 70 & non-linker IDRs ≥ 50: 625
Number of proteins with PDZ domains (and IDR and in cell junction): 52

Long IDRs refer to IDRs ≥ 70 and non-linker IDRs ≥ 50. The length cutoff allowed us to analyze IDRs that unambiguously confer functional preferences within
the context of full-length proteins with the chosen cutoff for the number of residues being equivalent to that of folded domains.

Output data of interest:
1. output/Human_IDRs_in_cell_junction_without_transmembrane_proteins.xlsx - List of all long IDRs and their clusters for proteins in cell junctions without transmembrane domains

2. output/Accessions_for_venn_diagram.xlsx - List of accessions for unique, proteins housing IDRs ≥ 70 and non-linker IDRs ≥ 50, list of accessions for unique proteins in cell junctions without transmembrane domains, and list of accessions for unique cell junction without transmembrane domain proteins housing IDRs ≥ 70 & non-linker IDRs ≥ 50 that can be used to plot the venn diagram

2. output/Human_Domains_in_proteins_w_IDRs_and_cell_junction_without_transmembrane_proteins.xlsx - List of Domains for all proteins with long IDRs and in cell junctions without transmembrane domains

3. output/Unique_protein_counts_by_domain_for_proteins_w_IDRs_and_in_cell_junctions_without_transmembrane_domains.xlsx - Unique protein counts by protein domain for all proteins with IDRs and in cell junctions without transmembrane domains. I.e., how many unique cell junction proteins housing a long IDR contain a given domain 

4. output/Count_of_binding_domains_in_each_protein_w_PDZ_and_IDR_and_in_cell_junction_without_transmembrane.xlsx - Count of number of binding domains in proteins with at least one PDZ, at least one long IDR, and in cell junctions but no transmembrane domain. Here binding domains are not necessarily PDZ but any protein binding domain

5. output/Count_of_IDRs_in_each_protein_w_PDZ_and_IDR_and_in_cell_junction_without_transmembrane.xlsx - Count of number of long IDRs in proteins with at least one PDZ, at least one IDR, and in cell junctions but no transmembrane domain






