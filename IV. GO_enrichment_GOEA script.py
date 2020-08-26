###       1. DOWNLOAD ONTOLOGIES AND ASSOCIATONS


### 1a. Download ontologies, if necessary

# Get http://geneontology.org/ontology/go-basic.obo

from goatools.base import download_go_basic_obo
obo_fname = download_go_basic_obo()

### 1b. Download Associations, if necessary

# Get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz

from goatools.base import download_ncbi_associations
file_gene2go = download_ncbi_associations()

###      2. LOAD ONTOLOGIES, ASSOCIATIONS AND BACKGROUND GENE SET

### 2a. Load Ontologies

from goatools.obo_parser import GODag
obodag = GODag("go-basic.obo")

### 2b. Load Associations

from __future__ import print_function
from goatools.anno.genetogo_reader import Gene2GoReader

# Read NCBI's gene2go. Store annotations in a list of namedtuples
objanno = Gene2GoReader(file_gene2go, taxids=[10090])

# Get associations for each branch of the GO DAG (BP, MF, CC)
ns2assoc = objanno.get_ns2assc()
for nspc, id2gos in ns2assoc.items():
    print("{NS} {N:,} annotated mouse genes".format(NS=nspc, N=len(id2gos)))

from goatools.cli.ncbi_gene_results_to_python import NCBIgeneToPythonCli

from genes_ncbi_10090_proteincoding import GENEID2NT as GeneID2nt_mus

### followed "https://github.com/tanghaibao/goatools/blob/1e93d26e4c93cb17786ab5fe736f90dc4f79421a/notebooks/backround_genes_ncbi.ipynb"
### to download a set of background population genes from NCBI.

from genes_ncbi_10090_proteincoding import GENEID2NT as GeneID2nt_mus

###       3. INITIALIZE A GOEA object

### The GOEA object holds the Ontologies, Associations, and background.
###  Numerous studies can then be run withough needing to re-load the above items.
### In this case, we only run one GOEA.

from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

goeaobj = GOEnrichmentStudyNS(
        GeneID2nt_mus,                             # List of mouse protein-coding genes
        ns2assoc,                                  # geneid/GO associations
        obodag,                                    # Ontologies
        propagate_counts = False,
        alpha = 0.05,                              # default significance cut-off
        methods = ['fdr_bh'])                      # defult multipletest correction method

###         4. READ STUDY GENES


# Data will be stored in this variable
import os
geneid2symbol = {}
# Get xlsx filename where data is stored
din_xlsx = r"C:\Users\krishna\Downloads\padj_converted.xlsx" ###excel file containing 3 columns:
                                                        ### gene_symbols (our test data), their respective ESENMBL gene ids, and their p adj values (test_data)

# Read data

if os.path.isfile(din_xlsx):
    import xlrd
    book = xlrd.open_workbook(din_xlsx)
    pg = book.sheet_by_index(0)
    for r in range(pg.nrows):
        symbol, geneid, pval = [pg.cell_value(r, c) for c in range(pg.ncols)]
        if geneid:
            geneid2symbol[int(geneid)] = symbol
    print('READ: {XLSX}'.format(XLSX=din_xlsx))
else:
    raise RuntimeError('CANNOT READ: {XLSX}'.format(XLSX=fin_xlsx))

###         5. Run Gene Ontology Enrichment Analysis (GOEA)

# 'p_' means "pvalue". 'fdr_bh' is the multipletest method we are currently using.
geneids_study = geneid2symbol.keys()
goea_results_all = goeaobj.run_study(geneids_study)
goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

### to export our analysis results one file with only gene symbols and second file with gene ids

goeaobj.wr_xlsx("GO_symbols.xlsx", goea_results_sig, itemid2name=geneid2symbol)
goeaobj.wr_xlsx("GO_geneids.xlsx", goea_results_sig)
