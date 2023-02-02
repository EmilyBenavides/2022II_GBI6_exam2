def fasta_downloader(id_coati): 
    """
    Funcion que pide como input un archivo que contenga los ids o codigo de acceso a la base de datos del GenBank
    y como output genere dos documentos tipo gb y fasta que contengan información de los ids utilizados
    """
    from Bio import Entrez
    in_sequence = open(id_coati, "r")
    out_sequence_gb = open("Data/coati.gb", "w")
    for linea in in_sequence: 
        Entrez.email = "emily.benavides@est.ikiam.edu.ec"
        handle=Entrez.efetch(db="nucleotide" ,id=linea ,rettype="gb", retmode="text")
        data=(handle.read())
        out_sequence_gb.write(data)
    out_sequence_gb.close()    
    
    in_sequence_fasta = open(id_coati, "r")
    out_sequence_fasta = open("Data/coati.fasta", "w")
    for linea in in_sequence_fasta: 
        Entrez.email = "emily.benavides@est.ikiam.edu.ec"
        handle=Entrez.efetch(db="nucleotide" ,id=linea ,rettype="fasta", retmode="text")
        data=(handle.read())
        out_sequence_fasta.write(data)
    out_sequence_fasta.close()

def alignmet (archivo_fasta):
    """
    Función que pide como input un archivo de secuencias multiples tipo fasta y que genere dos archivos que contengan 
    información de alineación y dendograma de las secuencias utilizadas. 
    """
    clustalw_exe = r"C:\Program Files (x86)\ClustalW2\clustalw2.exe"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile = file)
    assert os.path.isfile(clustalw_exe), "Clustal_W executable is missing or not found"
    stdout, stderr = clustalw_cline()
    print(clustalw_cline)
    ClustalAlign = AlignIO.read("Data/coati.aln", "clustal")
    print(ClustalAlign)
    tree = Phylo.read("Data/coati.dnd", "newick")

def tree (alineacion): 
    """
    Funcion que pide como entrada un archivo de alineación con extensión aln y como out un archivo pdf que contiene el arbol 
    filogenetico
    """
    from Bio.Phylo.TreeConstruction import DistanceCalculator 
    from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
    import matplotlib
    import matplotlib.pyplot as plt
    with open(alineacion,"r") as aln: 
        alignment = AlignIO.read(aln,"clustal")
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)
    constructor = DistanceTreeConstructor(calculator)
    # Construir el arbol 
    align_total = constructor.build_tree(alignment)
    align_total.rooted = True
    Phylo.write(align_total, "coati.xml", "phyloxml")
    align_protein = Phylo.read(file="coati.xml", format= "phyloxml")
    # Arbol elemental en Matplotlib
    #fig = Phylo.draw(cis_tree)
    fig = plt.figure(figsize=(30, 40), dpi=100) # create figure & set the size 
    matplotlib.rc('font', size=20)              # fontsize of the leaf and node labels 
    matplotlib.rc('xtick', labelsize=20)       # fontsize of the tick labels
    matplotlib.rc('ytick', labelsize=20)       # fontsize of the tick labels
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(align_protein, axes=axes)
    coati_phylotree.pdf
    plt.savefig('Data/coati_phylotree.pdf')  
    return()

    
