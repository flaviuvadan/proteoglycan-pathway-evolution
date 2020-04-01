# a file that holds constants relevant to the alignment domain


class Genes:
    """ A namespace of gene names """
    ACAN = "Acan"
    ARSA = "Arsa"
    ARSB = "Arsb"
    ARSC = "Arsc"
    ARSD = "Arsd"
    ARSE = "Arse"
    ARSF = "Arsf"
    ARSG = "Arsg"
    ARSI = "Arsi"
    ARSJ = "Arsj"
    ARSK = "Arsk"
    B3GALT6 = "B3galt6"
    B3GAT3 = "B3gat3"
    B4GALT7 = "B4galt7"
    CHPF = "Chpf"
    CHPF2 = "Chpf2"
    CHST1 = "Chst1"
    CHST11 = "Chst11"
    CHST12 = "Chst12"
    CHST13 = "Chst13"
    CHST14 = "Chst14"
    CHST15 = "Chst15"
    CHST2 = "Chst2"
    CHST3 = "Chst3"
    CHST4 = "Chst4"
    CHST7 = "Chst7"
    CHST9 = "Chst9"
    CHSY1 = "Chsy1"
    CHSY3 = "Chsy3"
    CSGALNACT1 = "Csgalnact1"
    CSGALNACT2 = "Csgalnact2"
    DSE = "Dse"
    EXT1 = "Ext1"
    EXT2 = "Ext2"
    FAM20A = "Fam20a"
    FAM20B = "Fam20b"
    FAM20C = "Fam20c"
    GALNS = "Galns"
    GLB1 = "Glb1"
    GNS = "Gns"
    GUSB = "Gusb"
    IDS = "Ids"
    PRG4 = "Prg4"
    SGSH = "Sgsh"
    SULF1 = "Sulf1"
    SULF2 = "Sulf2"
    SUMF1 = "Sumf1"
    SUMF2 = "Sumf2"
    UST = "Ust"
    XYLT1 = "Xylt1"
    XYLT2 = "Xylt2"


class GeneFunctions:
    """ A namespace of gene functions """
    CORE_PROTEIN = "CoreProtein"
    XYLOSYLTRANSFERASE = "Xylosyltransferase"
    GALACTOSYLTRANSFERASE = "Galactosyltransferase"
    GLUCURONYLTRANSFERASE = "Glucuronyltransferase"
    GLYCOSYLTRANSFERASE = "Glycosyltransferase"
    SULFOTRANSFERASE = "Sulfotransferase"
    SULFATASE = "Sulfatase"
    KINASE = "Kinase"
    EPIMERASE = "Epimerase"
    GLYCOSIDASE = "Glycosidase"
    SULFOHYDROLASE = "Sulfohydrolase"


# a map of gene to gene function, used for creating subplots of genes based on their function
GENE_FUNCTIONS = {
    Genes.ACAN: GeneFunctions.CORE_PROTEIN,
    Genes.ARSA: GeneFunctions.SULFATASE,
    Genes.ARSB: GeneFunctions.SULFATASE,
    Genes.ARSC: GeneFunctions.SULFATASE,
    Genes.ARSD: GeneFunctions.SULFATASE,
    Genes.ARSE: GeneFunctions.SULFATASE,
    Genes.ARSF: GeneFunctions.SULFATASE,
    Genes.ARSG: GeneFunctions.SULFATASE,
    Genes.ARSI: GeneFunctions.SULFATASE,
    Genes.ARSJ: GeneFunctions.SULFATASE,
    Genes.ARSK: GeneFunctions.SULFATASE,
    Genes.B3GALT6: GeneFunctions.GALACTOSYLTRANSFERASE,
    Genes.B3GAT3: GeneFunctions.GLUCURONYLTRANSFERASE,
    Genes.B4GALT7: GeneFunctions.GALACTOSYLTRANSFERASE,
    Genes.CHPF: GeneFunctions.GLUCURONYLTRANSFERASE,
    Genes.CHPF2: GeneFunctions.GLUCURONYLTRANSFERASE,
    Genes.CHST1: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST11: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST12: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST13: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST14: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST15: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST2: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST3: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST4: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST7: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHST9: GeneFunctions.SULFOTRANSFERASE,
    Genes.CHSY1: GeneFunctions.GALACTOSYLTRANSFERASE,
    Genes.CHSY3: GeneFunctions.GALACTOSYLTRANSFERASE,
    Genes.CSGALNACT1: GeneFunctions.GALACTOSYLTRANSFERASE,
    Genes.CSGALNACT2: GeneFunctions.GALACTOSYLTRANSFERASE,
    Genes.DSE: GeneFunctions.EPIMERASE,
    Genes.EXT1: GeneFunctions.GLYCOSYLTRANSFERASE,
    Genes.EXT2: GeneFunctions.GLYCOSYLTRANSFERASE,
    Genes.FAM20A: GeneFunctions.KINASE,
    Genes.FAM20B: GeneFunctions.KINASE,
    Genes.FAM20C: GeneFunctions.KINASE,
    Genes.GALNS: GeneFunctions.SULFATASE,
    Genes.GLB1: GeneFunctions.GLYCOSIDASE,
    Genes.GNS: GeneFunctions.SULFATASE,
    Genes.GUSB: GeneFunctions.GLYCOSIDASE,
    Genes.IDS: GeneFunctions.SULFATASE,
    Genes.PRG4: GeneFunctions.CORE_PROTEIN,
    Genes.SGSH: GeneFunctions.SULFOHYDROLASE,
    Genes.SULF1: GeneFunctions.SULFATASE,
    Genes.SULF2: GeneFunctions.SULFATASE,
    Genes.SUMF1: GeneFunctions.SULFATASE,
    Genes.SUMF2: GeneFunctions.SULFATASE,
    Genes.UST: GeneFunctions.SULFOTRANSFERASE,
    Genes.XYLT1: GeneFunctions.XYLOSYLTRANSFERASE,
    Genes.XYLT2: GeneFunctions.XYLOSYLTRANSFERASE,
}
